from setuptools import setup, find_packages
from biocodices import __version__


dependencies = [
    'termcolor',
    'pandas',
    'numpy',
    'matplotlib',
    'seaborn',
    'sqlalchemy',
    'myvariant',
    'inflect',
    'docopt',
    'pymysql',
    'python-ternary',
    'pyvcf',
]

setup(name='biocodices',
      version=__version__,
      description='Nac & Pop Next-Generation Sequencing',
      url='http://github.com/biocodices/biocodices',
      author='Juan Manuel Berros',
      author_email='juanma.berros@gmail.com',
      install_requires=dependencies,
      license='MIT',
      packages=find_packages(),
      entry_points={
          'console_scripts': [
            'bioco = biocodices.__main__:cli'
          ]
      },
      zip_safe=False)
