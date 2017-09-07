from setuptools import setup, find_packages

setup(name='rnaseq-lib',
      version='1.0a12',
      description='Library of convenience functions related to current research',
      url='http://github.com/jvivian/rnaseq-lib',
      author='John Vivian',
      author_email='jtvivian@gmail.com',
      license='MIT',
      package_dir={'': 'src'},
      packages=find_packages('src'),
      package_data={'rnaseq_lib.utils': ['data/*']},
      install_requires=['pandas',
                        'numpy',
                        'seaborn',
                        'requests',
                        'sklearn'])
