from setuptools import setup, find_packages

setup(
    entry_points={
          'console_scripts': [
              'pipeline=pipeline.pipeline:run',
          ]
      },
)