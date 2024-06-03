import numpy as np
import matplotlib.pyplot as plt
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 09:54:21 2024

@author: Christophe
"""

def example_function(x, y):
    
    f=3*(1-x)**2*np.exp(-x**2-(y+1)**2) -10*(x/5-x**2-y**5)*np.exp(-x**2-y**2)-1/3*np.exp(-(x+1)**2-y**2)
    f1=np.cos(10*y)+np.sin(10*(x-y))
    f2=np.exp(-(5-10*x)**2/2)+0.75*np.exp(-(5-10*y)**2/2)+0.75*np.exp(-(5-10*x)**2/2)*np.exp(-(5-10*y)**2/2)
    f3=np.sin(2*np.pi*y)*np.sin(np.pi*x)
    f4=0.75*np.exp(-((9*x-2)**2+(9*y-2)**2)/4)+0.75*np.exp(-(9*x+1)**2/49-(9*y+1)/10)+0.5*np.exp(-((9*x-7)**2+(9*y-3)**2)/4)+0.2*np.exp(-((9*x-4)**2-(9*y-7)**2))  
    f5=1/9*(np.tanh(9*y-9*x)+1)
    f6=(1.25+np.cos(5.4*y))/(6*(1+(3*x-1)**2))
    f7=1/3*np.exp(-(81/16)*((x-0.5)**2+(y-0.5)**2) )
    f8=1/3*np.exp(-(81/4)*((x-0.5)**2+(y-0.5)**2) )
    f9=3+2*np.sin(2*np.pi*x)*np.sin(2*np.pi*y+13)
    
    
    return f1
