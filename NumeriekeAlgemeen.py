"""
Author Frank Otto
Editted by Bouke Jansen, Ravi de Berg

In this document the schrodinger equation is solved on a discrete space, and the eigenValues for the energie are calculated.

WARNING: Dragging sliders might cause figures to crash!!
"""

import numpy as np
import matplotlib.pyplot as plt
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


def calculateEigenstate(hamiltonMatrix, gridSize, boxSize):
    """
    A function that calculates the eigenvalues and eigenvectors of a matrix.
    
    :param hamiltonMatrix: A numpy array representing the hamilton matrix for which you wish to know the eigenvalues.
    :param gridSize: Eigenlijk hoeft dit niet als parameter meegegeven want kan je berekenen uit hamiltonMatrix.
    :param boxSize: Length of the box
    :return: A size (2,) tuple with the first element eigenValues and second element normalized eigenVectors.
    """
    eigenValues, eigenVectors = np.linalg.eigh(hamiltonMatrix)

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


def potMaker(n1, a1, n2, a2, gridSize):
    """
    gooi dit ff lekker in je eigen task
    
    :param n1: 
    :param a1: 
    :param n2: 
    :param a2: 
    :param gridSize: 
    :return: 
    """
    def potentiaal(n, a, gridSize):
        matrix = np.diag(np.ones(gridSize))
        for i in range(gridSize):
            for j in range(gridSize):
                matrix[i][j] = matrix[i][j] * a * np.abs(0.5 * gridSize - j - 0.5) ** n
        return matrix

    potentiaalMatrix = potentiaal(n1, a1, gridSize) + potentiaal(n2, a2, gridSize) + getHamiltonMatrix(gridSize, 1)
    values, vectors = calculateEigenstate(potentiaalMatrix, gridSize, 1)
    return (values, vectors, 1)

def task2():
    gridSize = 400
    boxSize = 1
    t = 1
    
    V0 = 2
    deltapotential = np.zeros(gridSize)
    deltapotential[int(gridSize / 2)] = t * V0
    
    Hamiltonmatrix = getHamiltonMatrix(gridSize,t) - deltapotential
    eigenvalues,eigenvectors = calculateEigenstate(Hamiltonmatrix,gridSize,boxSize)
    plotEigenstate(eigenvalues,eigenvectors,boxSize,'Task 2.0')
    plt.show()

def task3():
    gridSize = 400
    wellDepth = 1
    wellWidth = 10
    t = 1
    boxSize = 1

    wellPosition = int(gridSize / 2 - wellWidth / 2)
    hamiltonMatrix = getHamiltonMatrix(gridSize, t)
    for i in range(wellWidth):
        hamiltonMatrix[wellPosition + i][wellPosition + i] -= wellDepth

    values, vectors = calculateEigenstate(hamiltonMatrix, gridSize, boxSize)
    plotEigenstate(values, vectors, boxSize, 'Task 3.0')
    plt.show()


def task5():
    plotEigenstate(potMaker(2, 10 ** -4, 0, 0, 101)[0], potMaker(2, 10 ** -4, 0, 0, 101)[1], 1, 'Task 5.0')
    for i in np.arange(0, 10 ** -5, 10 ** -6):
        hulpje = potMaker(2, 10 ** -4, 4, i, 101) # Slechte naam voor een variabele
        # plotEigenstate(hulpje[0],hulpje[1],hulpje[2])

        if i == 0:  # Aangezien deze statement altijd True geeft bij de eerste stap in je for kan je wat erin staat
            # misschien beter voor je for loop uitvoeren
            aMatrix = hulpje[0]
        if i != 0:
            aiMatrix = hulpje[0] - aMatrix
            plt.figure('Task 5.2')
            plt.plot(aiMatrix, 'o', label="a = " + str(i), markersize=2.5)
            plt.legend()
            plt.xlabel('n')
            plt.ylabel('$\Delta$E (A.U)')
            plt.title('Energieverschillen van deltapotential(x) = ax\u2074 + x\u00B2 en deltapotential(x) = x\u00B2')
        plt.figure('Task 5.1')
        plt.plot(hulpje[0], 'o', label="a = " + str(i), markersize=2.5)
        plt.legend()
        plt.xlabel('n')
        plt.ylabel('E (A.U)')
        plt.title('Energie met deltapotential(x) = ax\u2074 + x\u00B2')
    plt.show()


values, vectors = calculateEigenstate(getHamiltonMatrix(100, 1), 100, 1)
plotEigenstate(values, vectors, 1, 'Task 0.0')
plt.show()
task2()
task3()
task5()
