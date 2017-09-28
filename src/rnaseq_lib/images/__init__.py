import os
from subprocess import call


def resize_image(img_path, resize='28x28', grayscale=True, output_prefix='resized.'):
    """
    Resizes and image using ImageMagick

    :param str img_path: Path to image
    :param str resize: Resize w x h in pixels
    :param bool grayscale: If True, image will be greyscale
    :param str output_prefix: Prefix for output image
    :return: Output path
    :rtype: str
    """
    dirname, name = os.path.split(img_path)

    params = ['convert', img_path,
              '-resize', resize + '^',
              '-gravity', 'center',
              '-extent', resize]

    if grayscale:
        params.extend(['-type', 'Grayscale'])

    output_path = os.path.join(dirname, output_prefix + name)
    call(params + [output_path])

    return output_path
