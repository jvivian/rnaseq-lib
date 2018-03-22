import multiprocessing
import os
import shutil

from rnaseq_lib.docker import docker_call, docker_path
from rnaseq_lib.docker.tools import gdc_version, samtools_version, picardtools_version


def download_bam_from_gdc(gdc_id, token, work_dir='.'):
    """
    Downloads BAM file from the GDC using a GDC ID and GDC access token

    :param str gdc_id: GDC ID to download
    :param str token: File containing token credentials to GDC
    :param str work_dir: Directory being mounted into Docker
    :return: Path to BAM
    :rtype: str
    """
    # Move token to local work_dir for docker mount
    assert token, 'gdc_token is missing which is required for downloading. Check config.'
    shutil.copy(os.path.abspath(token), work_dir)

    # Define parameters and call tool
    parameters = ['download',
                  '-d', '-data',
                  '-t', '/data/{}'.format(os.path.basename(token)),
                  gdc_id]
    docker_call(tool=gdc_version, parameters=parameters, work_dir=work_dir)
    files = [x for x in os.listdir(os.path.join(work_dir, gdc_id)) if x.lower().endswith('.bam')]
    assert len(files) == 1, 'More than one BAM found from GDC URL: {}'.format(files)
    bam_path = os.path.join(work_dir, gdc_id, files[0])
    return bam_path


def assert_bam_is_paired_end(bam_path, region='chr6'):
    """
    Confirm that a BAM is paired-end and not single-end. Raises an error if not paired-end

    :param str bam_path: Path to BAM
    :param str region: Region of the genome to select
    """
    # Check if BAM index exists, otherwise index BAM
    bam_no_ext = os.path.splitext(bam_path)[0]
    if not os.path.exists(bam_no_ext + '.bai') and not os.path.exists(bam_no_ext + '.bam.bai'):
        index_bam(bam_path)

    docker_bam_path = docker_path(bam_path)
    work_dir = os.path.dirname(os.path.abspath(bam_path))

    # Check for both "chr" and no "chr" format
    results = []
    regions = [region, 'chr' + region] if 'chr' not in region else [region, region.lstrip('chr')]
    for r in regions:
        parameters = ['view', '-c', '-f', '1',
                      docker_bam_path,
                      r]  # Chr6 chosen for testing, any region with reads will work
        out = docker_call(work_dir=work_dir, parameters=parameters, tool=samtools_version, check_output=True)
        results.append(int(out.strip()))
    assert any(x for x in results if x != 0), 'BAM is not paired-end, aborting run.'


def index_bam(bam_path):
    """
    Creates a BAM index (.bai) in the same directory as the BAM
    Indexing is necessary for viewing slices of the BAM

    :param str bam_path: Path to BAM
    """
    work_dir = os.path.dirname(os.path.abspath(bam_path))
    parameters = ['index', docker_path(bam_path)]
    docker_call(work_dir=work_dir, parameters=parameters, tool=samtools_version)


def convert_bam_to_fastq(bam_path, check_paired=True, ignore_validation_errors=True):
    """
    Converts BAM to a pair of FASTQ files

    :param str bam_path: Path to BAM
    :param bool check_paired: If True, checks whether BAM is paired-end
    :param bool ignore_validation_errors: If True, ignores validation errors from picardTools
    :return: Paths to R1 and R2
    :rtype: tuple(str, str)
    """
    if check_paired:
        assert_bam_is_paired_end(bam_path)

    work_dir = os.path.dirname(os.path.abspath(bam_path))
    parameters = ['SamToFastq', 'I={}'.format(docker_path(bam_path)), 'F=/data/R1.fq', 'F2=/data/R2.fq']
    if ignore_validation_errors:
        parameters.append('VALIDATION_STRINGENCY=SILENT')
    docker_call(work_dir=work_dir, parameters=parameters, tool=picardtools_version)
    return os.path.join(work_dir, 'R1.fq'), os.path.join(work_dir, 'R2.fq')


def sort_bam(bam_path, cores=None):
    """
    Sorts STAR's output BAM using samtools

    :param str bam_path: Path to BAM file
    :param int cores: If provided, uses that number of cores. Otherwise uses all available cores
    """
    # Unpack bam directory and name
    work_dir, bam_basename = os.path.split(bam_path)
    bam_name, ext = os.path.splitext(bam_basename)

    if not cores:
        cores = multiprocessing.cpu_count()

    output_name = '{}.sorted.bam'.format(bam_name)
    parameters = ['sort',
                  '-o', output_name,
                  '-O', os.path.join('/data/', output_name),
                  '-T', 'temp',
                  '-@', str(cores),
                  os.path.join('/data/', bam_basename)]

    docker_call(tool=samtools_version, parameters=parameters, work_dir=work_dir)
    return os.path.join(work_dir, output_name)
