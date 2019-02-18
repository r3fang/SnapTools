from setuptools import setup

setup(name='snaptools',
      version='1.0',
      description='A module for working with snap files in Python',
      url='https://github.com/r3fang/snaptools.git',
      author='Rongxin Fang',
      author_email='r4fang@gmail.com',
      license='MIT',
      packages=['snaptools'],
      install_requires=[
          "pysam",
          "h5py",
          "numpy",
          "pybedtools"
      ],
      zip_safe=False)
