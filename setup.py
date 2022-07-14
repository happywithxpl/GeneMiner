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


'''compare version'''
def compare_version(a: str, b: str):
    lena = len(a.split('.'))  # Get the components of the version string
    lenb = len(b.split('.'))
    a2 = a + '.0' * (lenb-lena)  # Completing a when b is longer than a
    b2 = b + '.0' * (lena-lenb)

    for i in range(max(lena, lenb)):  # To compare each part, you need to convert to integers for comparison
        if int(a2.split('.')[i]) > int(b2.split('.')[i]):
            return 1
        elif int(a2.split('.')[i]) < int(b2.split('.')[i]):
            return 0
        else:						# If they are equal at the end of the comparison, the first version is returned
            if i == max(lena, lenb)-1:
                return 1



# python libs
install_dependencies = []
try:
    import Bio
    biopython_current=Bio.__version__
    biopython_requirement="1.70"
    flag=compare_version(biopython_current, biopython_requirement)
    if not flag:
        install_dependencies.append("biopython>=1.70")     
except ImportError:
    install_dependencies.append("biopython>=1.70")
else:
    sys.stdout.write("Existed module biopython " + str(Bio.__version__) + "\n")



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