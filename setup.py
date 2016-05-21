from setuptools import setup, find_packages
from biocodices import __version__


setup(name='biocodices',
      version=__version__,
      description='Nac & Pop Next-Generation Sequencing',
      url='http://github.com/biocodices/biocodices',
      author='Juan Manuel Berros',
      license='MIT',
      packages=find_packages(),
      zip_safe=False)
