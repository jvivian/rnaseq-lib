import os
from subprocess import call


def base_docker_call(mount):
    """
    Returns the boilerplate array used for Docker calls

    :param str mount: Directory to mount
    :return: Docker run parameters
    :rtype: list(str)
    """
    return ['docker', 'run', '--rm', '-v', '{}:/data'.format(os.path.abspath(mount))]


def fix_directory_ownership(output_dir, tool):
    """
    Uses a Docker container to change ownership recursively of a directory

    :param str output_dir: Directory to change ownership of
    :param str tool: Docker tool to use
    """
    stat = os.stat(output_dir)
    call(['docker', 'run', '--rm', '--entrypoint=chown', '-v', '{}:/data'.format(output_dir),
          tool, '-R', '{}:{}'.format(stat.st_uid, stat.st_gid), '/data'])
