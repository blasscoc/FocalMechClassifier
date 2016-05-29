from setuptools import setup

setup(name = "focal_mech",
      version = "0.1",
      author = "Ben Lasscock",
      author_email = "blasscoc@gmail.com",
      description = ("An Earthquake classifier."),
      packages=['focal_mech'],
      install_requires = ['numpy', 'scipy', 'matplotlib',
                          'mplstereonet',
                          'obspy'])


