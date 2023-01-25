from setuptools import setup, find_packages

setup(
    entry_points={
          'console_scripts': [
              'scRNA2BoNI=scRNA2BoNI.scRNA2BoNI:run',
          ]
      },
)