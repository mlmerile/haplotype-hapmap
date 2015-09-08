from setuptools import setup, find_packages

setup(name='pyhapmap',
      version='0.1',
      description='Convert Hapmap data to numpy array / pandas dataframe',
      author='Thomas Dias-Alves',
      packages=find_packages(),
      test_suite='pyhapmap.tests',
      zip_safe=False)
