from setuptools import setup, find_packages

setup(name='rnaseq-lib',
      version='1.0a5',
      description='',
      url='http://github.com/jvivian/rnaseq-lib',
      author='John Vivian',
      author_email='jtvivian@gmail.com',
      license='MIT',
      package_dir={'': 'src'},
      packages=find_packages('src'),
      install_requires=['pandas',
                        'numpy',
                        'seaborn'])