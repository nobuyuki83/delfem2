from setuptools import setup, find_packages

setup(name='delfem2',
      version='0.0.0',
      description='Handy toolset for implementing geometry processing and finite element simulation',
      author='Nobuyuki Umetani',
      author_email='n.umetani@gmail.com',
      url='https://github.com/nobuyuki83/delfem2', 
      package_dir={'dfm2': 'python/dfm2'},
      packages=['dfm2'])