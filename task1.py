"""
author: Geert Schulpen
Task 1; random disorder
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import random as random

g_AxColor = "lightgoldenrodyellow"

random.seed(1248)

def disorder(size,scale):
    """
    A function that generates a disorder-potential
    :param size: how big the matrix has to be
    :param scale: how big the disorder has to be
    :return: A square array of size size, which contains random entries on the main diagonal
    """
    deltaHhatDiagonaal=[]
    deltaHhatDiagonaal.append([random.random() for i in range(size)])*scale 
    deltaHhat=np.diagflat(deltaHhatDiagonaal) #aanmaken extra matrix
    return(deltaHhat)
    
def task1():
    gridSize = 400
    t=1
    box=1
    big=1
    medium=0.075
    small=0.00125
    
    standardHamilton=getHamiltonMatrix(gridSize,t)
    deltaHamiltonSmall=disorder(gridSize,small)
    deltaHamiltonMedium=disorder(gridSize,medium)
    deltaHamiltonBig=disorder(gridSize,big)
    
    valuesSmall,vectorsSmall=calculateEigenstate(standardHamilton+deltaHamiltonSmall, box)
    valuesMedium,vectorsMedium=calculateEigenstate(standardHamilton+deltaHamiltonMedium, box)
    valuesBig,vectorsBig=calculateEigenstate(standardHamilton+deltaHamiltonBig, box)
    
    plotEigenstate(valuesSmall,vectorsSmall,box,'Task1, small disorder')
    plotEigenstate(valuesMedium,vectorsMedium,box,'Task1, medium disorder')
    plotEigenstate(valuesBig,vectorsBig,box,'Task1, big disorder')


task1()    
