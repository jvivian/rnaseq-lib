import os
from subprocess import call

import logging
logging.basicConfig(level=logging.INFO)


def filter_aligned_reads(bam_path, output_name=None, paired=True):
    """
    Filters bams for aligned reads

    :param str bam_path: Path to bam file
    :param str output_name: Defaults to input bam with '.filtered' inserted before the .bam extension
    :param bool paired: Only keep paired reads
    :return: Path to filtered bam
    :rtype: str
    """
    if not output_name:
        output_name = os.path.basename(os.path.splitext(bam_path)[0]) + '.filtered.bam'

    # Define parameters
    work_dir = os.path.dirname(os.path.abspath(bam_path))
    parameters = ['docker', 'run',
                  '-v', '{}:/data'.format(work_dir),
                  'quay.io/ucsc_cgl/samtools',
                  'view',
                  '-h',
                  '-o', '/data/{}'.format(output_name),
                  '-F', '0x04 ']
    if paired:
        parameters.extend(['-f', '0x02'])
    parameters.append(os.path.join('/data', os.path.basename(bam_path)))

    # Call tool
    output_path = os.path.join(work_dir, output_name)
    if not os.path.exists(output_path):
        logging.info('Filtering bam')
        call(parameters)
    else:
        logging.info('Skipping. Filtered bam exists: {}'.format(output_path))
    return output_path


def make_test_bam(bam_path, region='chr6', output_name='chr6.test.bam', clean_input_bam=False):
    """
    Makes a smaller bam based on the region argument passed

    :param str bam_path: Path to bam to use to make test bam
    :param str region: Region of the genome to subset from
    :param str output_name: Output file name test bam
    :param bool clean_input_bam: Cleans the input bam before starting
    :return: Path to test bam
    :rtype: str
    """
    if clean_input_bam:
        bam_path = filter_aligned_reads(bam_path)

    work_dir = os.path.dirname(os.path.abspath(bam_path))
    docker_parameters = ['docker', 'run',
                         '-v', '{}:/data'.format(work_dir),
                         'quay.io/ucsc_cgl/samtools']

    bam_no_ext = bam_path.split('.bam')[0]
    if not os.path.exists(bam_no_ext + '.bai') and not os.path.exists(bam_no_ext + '.bam.bai'):
        index = docker_parameters + ['index', os.path.join('/data', os.path.basename(bam_path))]
        call(index)

    parameters = docker_parameters + ['view',
                                      '-b',
                                      '-h',
                                      '-o', os.path.join('/data', output_name),
                                      os.path.join('/data', os.path.basename(bam_path)),
                                      region]
    call(parameters)
    return os.path.join(work_dir, output_name)  #filter_aligned_reads()
