"""
Author Frank Otto
Editted by Bouke Jansen

In this document the schrodinger equation is solved on a discrete space, and the eigenvalues for the energie are calculated.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

g_AxColor = "lightgoldenrodyellow"

# Construct base Hamilton matrix
def getHamiltonMatrix(gridSize, t):
	return np.diag(-t * np.ones(gridSize - 1), k = -1) \
		   + np.diag(2 * t * np.ones(gridSize)) \
		   + np.diag(-t * np.ones(gridSize - 1), k = 1)

def calculateAndPlot(hamiltonMatrix, gridSize, spacing):
	fig, (ax1, ax2, ax3) = plt.subplots(3, 1)

	# Calculate the eigenvalues
	eigenValues, eigenVectors = np.linalg.eigh(hamiltonMatrix)

	for i in range(gridSize):
		eigenVectors[:, i] /= np.sqrt(np.sum(eigenVectors[:, i]**2 / gridSize))

	ax1.plot(eigenValues, "o", markersize = 3)
	ax1.set_title("Energy")
	ax1.set_ylabel("E (A.U.)")
	ax1.set_xlabel("n")

	n0 = 0

	eigenVector = eigenVectors[:, n0]
	l, = ax2.plot(np.linspace(0, spacing, num = gridSize), eigenVector ** 2)
	ax2.set_ylabel(r"$|\psi|^2$")
	ax2.set_xlabel("x (A.U)")
	ax2.set_title("Probability density")

	l2, = ax3.plot(np.linspace(0, spacing, num = gridSize), eigenVector)
	ax3.set_ylim([np.min(eigenVectors), np.max(eigenVectors)])
	ax3.set_ylabel(r"$\psi$")
	ax3.set_xlabel("x (A.U)")
	ax3.set_title("Wave function")

	plt.tight_layout()

	sliderAx = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor = g_AxColor)
	slider = Slider(sliderAx, "n", 0, gridSize - 1, valinit = n0, valfmt = "%d")
	def update(val):
		n = int(slider.val)
		y = eigenVectors[:, n]
		l.set_ydata(y ** 2)
		l2.set_ydata(y)
		ax2.set_title("Probability distribution for eigenvalue = {:4.3f}".format(eigenValues[n]))
		fig.canvas.draw_idle()
	slider.on_changed(update)

	resetAx = plt.axes([0.75, 0.025, 0.1, 0.04])
	resetButton = Button(resetAx, "Reset", color = g_AxColor, hovercolor="0.975")
	def reset(event):
		slider.reset()
	resetButton.on_clicked(reset)

	plt.subplots_adjust(bottom = 0.25)

	closeAx = plt.axes([0.85, 0.025, 0.1, 0.04])
	closeButton = Button(closeAx, "Close", color = g_AxColor, hovercolor="0.975")
	def close(event):
		plt.close()
	closeButton.on_clicked(close)


#calculateAndPlot(getHamiltonMatrix(100, 1), 100, 1)

def task3():
	gridSize = 400
	wellDepth = 1
	wellWidth = 10
	t = 1
	spacing = 1

	wellPosition = int(gridSize/2 - wellWidth/2)
	hamiltonMatrix = getHamiltonMatrix(gridSize, t)
	for i in range(wellWidth):
		hamiltonMatrix[wellPosition + i][wellPosition + i] -= wellDepth
	calculateAndPlot(hamiltonMatrix, gridSize, spacing)
	plt.show()

def task5():
	for i in range(10):
		task5help(2, 50**-1, 4, i, 101)
	plt.show()


def task5help(n1, a1, n2, a2, gridSize):
	def potentiaal(n, a, gridSize):
		matrix = np.diag(np.ones(gridSize)) 
		for i in range(gridSize):
			for j in range(gridSize):
				matrix[i][j] = matrix[i][j] * a * np.abs(0.5 * gridSize - j - 0.5)**n
		return matrix
	potentiaalMatrix = potentiaal(n1, a1, gridSize) + potentiaal(n1, a2, gridSize) + getHamiltonMatrix(gridSize, 1)
	calculateAndPlot(potentiaalMatrix, gridSize, 1)
#task3()
task5()
