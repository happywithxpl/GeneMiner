#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2022/7/13 22:58
# @Author  : xiepulin
# @File    : setup.py
# @Software: PyCharm



"""Setup the GeneMiner environment."""

# Test pip
# 1) Clean the /dist directory
# 2) python3 setup.py sdist bdist_wheel
# 3) pip install --index-url https://test.pypi.org/simple/
#    --extra-index-url https://pypi.org/simple atram
# 4) twine upload --repository-url https://test.pypi.org/legacy/ dist/*


import platform
import sys
import setuptools
from setuptools import setup

# python libs
install_dependencies = []
try:
    import Bio
    print(Bio.__version__)
except ImportError:
    install_dependencies.append("biopython>=1.70")
else:
    sys.stdout.write("Existed module numpy " + str(Bio.__version__) + "\n")



sys.stdout.write("Python " + str(sys.version).replace("\n", " ") + "\n")
sys.stdout.write("PLATFORM: " + " ".join(platform.uname()) + "\n")
sys.stdout.write("Using setuptools " + str(setuptools.__version__) + "\n")

scripts_to_install = ["geneminer.py"]

setup(
    name="GeneMiner",
    version="1.0.0",
    author="xiepulin",
    author_email="xiepulin@stu.scu.edu.cn , happywithxpl@126.com",
    description="""GeneMiner : a software for extracting phylogenetic markers from next generation sequencing data""",
    license="GPL-3.0 license",
    # Project home
    url="https://github.com/happywithxpl/Geneminer",
    python_requires='>=3.7',
    install_requires=install_dependencies,
    scripts=scripts_to_install,
    packages=["lib"],
    zip_safe=False
)