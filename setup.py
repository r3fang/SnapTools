from setuptools import setup

snaptools_version = '1.2.9'

setup(
      name='snaptools',
      version=snaptools_version,
      author='Rongxin Fang',
      author_email='r4fang@gmail.com',
      license='LICENSE',
      packages=['snaptools'],
      description='A module for working with snap files in Python',
      url='https://github.com/r3fang/SnapTools.git',
      python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, !=3.3.*, <3.0',
      
      install_requires=[
          "pysam",
          "h5py",
          "numpy",
          "pybedtools==0.7.10",
          "python-louvain"
      ],
      keywords = ["Bioinformatics pipeline",
                  "Single cell analysis",
                  "Epigenomics",
                  "Epigenetics",
                  "ATAC-seq",
                  "Chromatin Accessibility",
                  "Functional genomics"],
      scripts = ["bin/snaptools"],      
      zip_safe=False)

if __name__ == '__main__':
    f = open("snaptools/__init__.py",'w')
    f.write("__version__ = \'"+snaptools_version+"\'"+"\n")
    f.close()
