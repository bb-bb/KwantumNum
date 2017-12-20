"""
Author Frank Otto
Editted by Bouke Jansen, Ravi de Berg

In this document the schrodinger equation is solved on a discrete space, and the eigenValues for the energie are calculated.

WARNING: Dragging sliders might cause figures to crash!!
"""

from quantummodule import plotEigenstate, calculateEigenstate, getHamiltonMatrix
import matplotlib.pyplot as plt
import numpy as np

def theory(n, x):
	return np.sqrt(2) * np.sin((n + 1) * x * np.pi)
values, vectors = calculateEigenstate(getHamiltonMatrix(100, 1), 1)
plotEigenstate(values, vectors, 1, 'Task 0.0', theory)

plt.show()

