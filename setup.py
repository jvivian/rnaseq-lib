from setuptools import setup, find_packages

setup(name='rnaseq-lib',
      version='1.0a1',
      description='',
      url='http://github.com/jvivian/rnaseq-lib',
      author='John Vivian',
      author_email='jtvivian@gmail.com',
      license='MIT',
      packages=find_packages(),
      install_requires=['pandas',
                        'numpy',
                        'seaborn'])