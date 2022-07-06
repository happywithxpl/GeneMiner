#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2022/7/6 16:59
# @Author  : xiepulin
# @File    : install.py
# @Software: PyCharm

import sys
import subprocess

#检测python版本
def check_python_version():
    this_python = sys.version_info[:2]
    min_version = (3, 6)
    if this_python < min_version:
        message_parts = [
            "GeneMiner does not work on Python {}.{}.".format(*this_python),            #*可以用来接受元组
            "The minimum supported Python version is {}.{}.".format(*min_version),
            "Please download the eligible version from https://www.python.org/.".format(*this_python)]
        print("ERROR: " + " ".join(message_parts))
        sys.exit(1)

#检测pip是否存在
def check_pip():
    try:
        import  pip
    except ImportError:
        message_parts=["pip does not detect in the system.",
                       "You can download the appropriate version from https://pip.pypa.io/en/stable/installation/ ."
                       ]
        print("ERROR: " + " ".join(message_parts))
        sys.exit(1)

#检测依赖并安装
def check_biopython():
    mirror_sources = ["https://mirrors.aliyun.com/pypi/simple/","https://pypi.douban.com/simple/", "https://pypi.mirrors.ustc.edu.cn/simple/"]

    biopython_required=[]
    try:
        import  Bio
        biopython_version=float(Bio.__version__)
        if biopython_version < 1.79:
            biopython_required.append("biopython==1.79")
    except ImportError:
        biopython_required.append("biopython==1.79")


    flag=0

    if biopython_required==[]:
        print("Requirements already satisfied")
    else:
        for i in mirror_sources:
            flag=install_biopython(biopython_required[0],i)
            if flag==1:
                break

        if flag==0:
            print("Try a common installation method ...")
            flag_v1=install_biopython_normal(biopython_required[0])
            if flag_v1==0:
                message_parts = ["Failed to install biopython , please check pip tool or install manually"
                                 ]
                print("ERROR: " + " ".join(message_parts))
                sys.exit(1)
            else:
                print("Requirements already satisfied")

        else:
            print("Requirements already satisfied")



def install_biopython(package,mirror):
    cmd="pip3 install {0} -i  {1} --default-timeout=100 --user ".format(package,mirror)
    subprocess.call(cmd,shell=True)
    try:
        import Bio
        flag=1
    except ImportError:
        flag=0
    return flag

def install_biopython_normal(package):
    cmd = "pip3 install {0} --default-timeout=100 --user ".format(package)
    subprocess.call(cmd, shell=True)
    try:
        import Bio
        flag = 1
    except ImportError:
        flag = 0
    return flag




if __name__ == '__main__':
    check_python_version()
    check_pip()
    check_biopython()






