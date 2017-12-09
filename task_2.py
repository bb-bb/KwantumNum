"""
Author Frank Otto
Task 2 the delta potential.
"""

from quantummodule import *
import matplotlib.pyplot as plt
import numpy as np


gridSize = 400
boxSize = 1
t = 1
deltapotential = np.zeros(gridSize)

for V0 in [-5,-2,0,2,5]:
    deltapotential[int(gridSize / 2)] = t * V0
    Hamiltonmatrix = getHamiltonMatrix(gridSize, t) - deltapotential
    eigenvalues, eigenvectors = calculateEigenstate(Hamiltonmatrix, boxSize)
    plotEigenstate(eigenvalues, eigenvectors, boxSize, 'Task 2.0, V0 = {:4.3f}'.format(V0))

    plt.show()
