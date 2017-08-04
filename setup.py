from setuptools import setup, find_packages

setup(name = "focal_mech",
      version = "0.1",
      author = "Ben Lasscock",
      author_email = "blasscoc@gmail.com",
      description = ("An Earthquake classifier."),
      packages=find_packages(),
      install_requires = ['numpy', 'scipy', 'matplotlib',
                          'mplstereonet',
                          'obspy'])


