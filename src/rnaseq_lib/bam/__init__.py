import os
import shutil

from rnaseq_lib.docker import docker_call


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
    docker_call(tool='jvivian/gdc-client:1.0', parameters=parameters, work_dir=work_dir)
    files = [x for x in os.listdir(os.path.join(work_dir, gdc_id)) if x.lower().endswith('.bam')]
    assert len(files) == 1, 'More than one BAM found from GDC URL: {}'.format(files)
    bam_path = os.path.join(work_dir, gdc_id, files[0])
    return bam_path
