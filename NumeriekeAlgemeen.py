"""
Author Frank Otto
Editted by Bouke Jansen, Ravi de Berg

In this document the schrodinger equation is solved on a discrete space, and the eigenValues for the energie are calculated.

WARNING: Dragging sliders might cause figures to crash!!
"""

from quantummodule import plotEigenstate, calculateEigenstate, getHamiltonMatrix
import matplotlib.pyplot as plt

values, vectors = calculateEigenstate(getHamiltonMatrix(100, 1), 1)
plotEigenstate(values, vectors, 1, 'Task 0.0')
plt.show()

