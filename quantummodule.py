"""
Author Frank Otto, Bouke Jansen, Ravi de Berg
General functions that are relevant for all tasks will be created in this file.
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider, Button

g_AxColor = "lightgoldenrodyellow"


def getHamiltonMatrix(gridSize, t):
    """
    :param gridSize: The size of the Hamiltonmatrix you wish to construct.
    :param t: The unit energy in hbar/8ma
    :return: A numpy array of size (gridSize,gridSize) representing the Hamiltonmatrix
    """
    return np.diag(-t * np.ones(gridSize - 1), k=-1) \
           + np.diag(2 * t * np.ones(gridSize)) \
           + np.diag(-t * np.ones(gridSize - 1), k=1)


def calculateEigenstate(hamiltonMatrix, boxSize):
    """
    A function that calculates the eigenvalues and eigenvectors of a matrix.

    :param hamiltonMatrix: A numpy array representing the hamilton matrix for which you wish to know the eigenvalues.
    :param gridSize: Eigenlijk hoeft dit niet als parameter meegegeven want kan je berekenen uit hamiltonMatrix.
    :param boxSize: Length of the box
    :return: A size (2,) tuple with the first element eigenValues and second element normalized eigenVectors.
    """
    eigenValues, eigenVectors = np.linalg.eigh(hamiltonMatrix)
    gridSize = len(eigenValues)

    # Normalization procedure
    for i in range(gridSize):
        eigenVectors[:, i] /= np.sqrt(boxSize * np.sum(eigenVectors[:, i] ** 2 / gridSize))

    return eigenValues, eigenVectors


def plotEigenstate(eigenValues, eigenVectors, boxSize, windowtitle):
    gridSize = len(eigenValues)

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
    fig.canvas.set_window_title(windowtitle)
    ax1.plot(eigenValues, "o", markersize=3)
    ax1.set_title("Energy")
    ax1.set_ylabel("E (a.u.)")
    ax1.set_xlabel("n")

    n0 = 0

    eigenVector = eigenVectors[:, n0]
    l, = ax2.plot(np.linspace(0, boxSize, num=gridSize), eigenVector ** 2)
    ax3.set_ylim([0, 1.1 * np.max(eigenVectors) ** 2])
    ax2.set_ylabel(r"$|\psi|^2$")
    ax2.set_xlabel("x (A.U)")
    ax2.set_title("Probability density")

    l2, = ax3.plot(np.linspace(0, boxSize, num=gridSize), eigenVector)
    ax3.set_ylim([1.1 * np.min(eigenVectors), 1.1 * np.max(eigenVectors)])
    ax3.set_ylabel(r"$\psi$")
    ax3.set_xlabel("x (A.U)")
    ax3.set_title("Wave function")

    plt.tight_layout()

    sliderAx = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=g_AxColor)
    slider = Slider(sliderAx, "n", 0, gridSize - 1, valinit=n0, valfmt="%d")

    def update(val):
        n = int(slider.val)
        y = eigenVectors[:, n]
        l.set_ydata(y ** 2)
        l2.set_ydata(y)
        ax2.set_title("Probability distribution for eigenvalue = {:4.3f}".format(eigenValues[n]))
        fig.canvas.draw_idle()

    slider.on_changed(update)

    resetAx = plt.axes([0.75, 0.025, 0.1, 0.04])
    resetButton = Button(resetAx, "Reset", color=g_AxColor, hovercolor="0.975")

    def reset(event):
        slider.reset()

    resetButton.on_clicked(reset)

    plt.subplots_adjust(bottom=0.25)

    closeAx = plt.axes([0.85, 0.025, 0.1, 0.04])
    closeButton = Button(closeAx, "Close", color=g_AxColor, hovercolor="0.975")

    def close(event):
        plt.close()

    closeButton.on_clicked(close)
