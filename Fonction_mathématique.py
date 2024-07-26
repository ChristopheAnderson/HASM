import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import random
import pandas as pd


def f(x, y):
    return 3+2*np.sin(2*np.pi*x)*np.sin(2*np.pi*y) + 13
    return np.sin(np.pi*x)*np.sin(2*np.pi*y)
    return 1/3*np.exp(-(81/4)*((x-0.5)**2+(y-0.5)**2) )
    return np.cos(10*y) + np.sin(10*(x-y))
    return x**2 - y**2
    return np.exp((81/(x-0.5)**2)+(y-0.5)**2)
    return x*y