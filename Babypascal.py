# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 16:54:38 2021

@author: Lara
"""
def Pascal_Mol(P, V=0.0001087, R = 8.3145, T = 277.15) :
    
    n = (P*V) / (R*T)
    
    print("For",P, "Pa, we need", n, "mol H2")
    return n


def Pascal_Pa( n, T = 277.15,  V=0.0001087, R = 8.3145) :
    
    P = (n*R*T) / V
    print("At",n, "mole H2, we have a Pressure of", P, "Pa")
    return P

