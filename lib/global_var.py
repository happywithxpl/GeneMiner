#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2021/12/17 16:58
# @Author  : xiepulin
# @File    : global_var.py
# @Software: PyCharm

def _init_():#初始化
    global _global_dict
    _global_dict = {}

def set_value(key,value):
    """ 定义一个全局变量 """
    global _global_dict
    _global_dict[key] = value

def get_value(key,defValue=None):
    """ 获得一个全局变量,不存在则返回默认值 """
    global _global_dict
    try:
        return _global_dict[key]
    except KeyError:
        return defValue