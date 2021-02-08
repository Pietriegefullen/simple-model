# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 16:54:38 2021

@author: Lara
"""
from scipy.constants import convert_temperature

def Pa_nach_Mol(P, V=0.0001087, R = 8.3145, T = 277.15) :
    
    n = (P*V) / (R*T)
    
    print("For",P, "Pa, we need", n, "mol H2")
    return n


def Mol_nach_Pa( n, T = 277.15,  V=0.0001057, R = 8.3145) :
    
    P = (n*R*T) / V
    print("At",n, "mole H2, we have a Pressure of", P, "Pa")
    return P

def g_nach_mol( g, gproMol, Substance) :
    
    MolKonz = g/gproMol
    KnFlasche =(0.0133 * MolKonz *10**6)  # 0,0133 ist das Volumen von Boden und wasser in unserer flasche
    
    
    print("At",g, "g", Substance, " pro L we have", MolKonz , "mol pro L", Substance)
    print("This would be ", KnFlasche, "Mikromol in our study")

# g pro mol : 
Acetate = 59.04

g_nach_mol(g=5, gproMol = Acetate, Substance = "Acetate")


