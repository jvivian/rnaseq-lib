from setuptools import setup, find_packages

setup(name='rnaseq-lib',
      version='1.0a8',
      description='',
      url='http://github.com/jvivian/rnaseq-lib',
      author='John Vivian',
      author_email='jtvivian@gmail.com',
      license='MIT',
      package_dir={'': 'src'},
      packages=find_packages('src'),
      package_data={'rnaseq_lib.utils.data': ['samples.pickle',
                                              'UCSF-genes.csv',
                                              'cancer-MAB-gene-targets.txt',
                                              'gene_map.pickle']},
      install_requires=['pandas',
                        'numpy',
                        'seaborn',
                        'requests'])
