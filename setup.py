#!/usr/bin/python
# #
# # setup(name='yUMItools',
# #       description='UMI tools for ERR-Seq Project',
# #       version='1.0.0',
# #       author_email=anitido@g.harvard.edu,
# #       url='https://github.com/adamn102/yUMItools/',
# #       packages=['yUMItools'])
# #
#
# setuptools.setup(name='yUMItools')

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="yUMItools",
    version="1.0.3",
    author="Adam Nitido",
    author_email="anitido@g.harvard.edu",
    description="UMI tools for ERR-Seq Project",
    long_description=long_description,
    url="https://github.com/adamn102/yUMItools/",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)