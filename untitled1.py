# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 15:38:39 2020

@author: Lara
"""

b = 3


def func(par):
    a = 2.0
    
    c = a + par
    
    def func2(foo):
        print(foo)
        d = 10
        return d
    
    
    func2(par)
    
    return 1.


func(b)