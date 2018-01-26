"""
Authors Frank Otto, Bouke Jansen, Ravi de Berg
General functions that are relevant for all tasks will be created in this file.
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider, Button

g_AxColor = "lightgoldenrodyellow"

def getHamiltonMatrix(gridSize, t):
	"""
	:param gridSize: The size of the Hamiltonmatrix you wish to construct.
	:param t: The unit energy in hbar/8ma.
	:return: A numpy array of size (gridSize,gridSize) representing the Hamiltonmatrix.
	"""
	return np.diag(-t * np.ones(gridSize - 1), k = -1) \
		   + np.diag(2 * t * np.ones(gridSize)) \
		   + np.diag(-t * np.ones(gridSize - 1), k = 1)


def calculateEigenstate(hamiltonMatrix, boxSize):
	"""
	A function that calculates the eigenvalues and eigenvectors of a matrix.

	:param hamiltonMatrix: A numpy array representing the hamilton matrix for which you wish to know the eigenvalues.
	:param boxSize: Length of the box
	:return: A size (2,) tuple with the first element eigenValues and second element normalized eigenVectors.
	"""
	eigenValues, eigenVectors = np.linalg.eigh(hamiltonMatrix)
	gridSize = len(eigenValues)

	# Normalization procedure
	for i in range(gridSize):
		eigenVectors[:, i] /= np.sqrt(boxSize * np.sum(eigenVectors[:, i] ** 2 / gridSize))

	return eigenValues, eigenVectors


def plotEigenstate(eigenValues, eigenVectors, boxSize, windowTitle, waveFunc = False):
	"""
	A function plotting the eigenValues, an eigenStates and eigenState squared,
	with some slider to scroll through the eigenStates and eigenState squared.

	:param eigenValues: Sorted set of eigenvalues.
	:param eigenVectors: Eigenvectors in the order to their corresponding eigenvalues.
	:param boxSize: Length of the box.
	:waveFunc: Optional argument, a function to plot through the eigenstates. It should take two arguments:
	n the eigenstate number and x the position.

	WARNING: Dragging sliders might cause figures to crash!!
	"""
	gridSize = len(eigenValues)

	fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
	fig.canvas.set_window_title(windowTitle)
	ax1.plot(eigenValues, "o", markersize = 3)
	ax1.set_title("Energy")
	ax1.set_ylabel("E")
	ax1.set_xlabel("n")

	n0 = 0

	x = np.linspace(0, boxSize, num=gridSize)
	eigenVector = eigenVectors[:, n0]
	l, = ax2.plot(x, eigenVector ** 2)
	ax2.set_ylim([0, 1.1 * np.max(eigenVectors) ** 2])
	ax2.set_ylabel(r"$|\psi|^2$")
	ax2.set_xlabel("x")
	ax2.set_title("Probability density")

	l2, = ax3.plot(x, eigenVector)
	ax3.set_ylim([1.1 * np.min(eigenVectors), 1.1 * np.max(eigenVectors)])
	ax3.set_ylabel(r"$\psi$")
	ax3.set_xlabel("x")
	ax3.set_title("Wave function")

	if waveFunc != False:
		waveFuncData = waveFunc(0, x)
		lb, = ax2.plot(x, waveFuncData ** 2)
		ax2.legend([l, lb], ["Distcrete", "Continuous"], loc = "right", framealpha = 0.5)
		l2b, = ax3.plot(x, waveFuncData)
		ax3.legend([l2, l2b], ["Distcrete", "Continuous"], loc = "right", framealpha = 0.5)

	plt.tight_layout()

	sliderAx = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor = g_AxColor)
	slider = Slider(sliderAx, "n", 0, gridSize - 1, valinit = n0, valfmt="%d")

	def update(val):
		n = int(slider.val)
		y = eigenVectors[:, n]
		l.set_ydata(y ** 2)
		l2.set_ydata(y)
		if waveFunc != False:
			yb = waveFunc(n, np.linspace(0, boxSize, num = gridSize))
			lb.set_ydata(yb**2)
			l2b.set_ydata(yb)
		ax2.set_title("Probability distribution for eigenvalue = {:4.3f}".format(eigenValues[n]))
		fig.canvas.draw_idle()

	slider.on_changed(update)

	plt.subplots_adjust(bottom=0.25)

	closeAx = plt.axes([0.85, 0.025, 0.1, 0.04])
	closeButton = Button(closeAx, "Close", color = g_AxColor, hovercolor="0.975")

	def close(event):
		plt.close()

	resetAx = plt.axes([0.65, 0.025, 0.1, 0.04])
	resetButton = Button(resetAx, "Reset", color = g_AxColor, hovercolor="0.975")

	def reset(event):
		slider.reset()

	resetButton.on_clicked(reset)

	closeButton.on_clicked(close)

	nextAx = plt.axes([0.55, 0.025, 0.1, 0.04])
	nextButton = Button(nextAx, "Next", color = g_AxColor, hovercolor="0.975")

	def next(event):
		slider.set_val((slider.val + 1) % gridSize)

	nextButton.on_clicked(next)

	prevAx = plt.axes([0.45, 0.025, 0.1, 0.04])
	prevButton = Button(prevAx, "Previous", color = g_AxColor, hovercolor="0.975")

	def previous(event):
		slider.set_val((slider.val - 1) % gridSize)

	prevButton.on_clicked(previous)
