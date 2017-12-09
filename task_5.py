from quantummodule import *
import numpy as np
import matplotlib.pyplot as plt


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
    values, vectors = calculateEigenstate(potentiaalMatrix, 1)
    return (values, vectors, 1)


def task5():
    plotEigenstate(potMaker(2, 10 ** -4, 0, 0, 101)[0], potMaker(2, 10 ** -4, 0, 0, 101)[1], 1, 'Task 5.0')
    for i in np.arange(0, 10 ** -5, 10 ** -6):
        hulpje = potMaker(2, 10 ** -4, 4, i, 101)  # Slechte naam voor een variabele
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
