from quantummodule import *
import numpy as np
import matplotlib.pyplot as plt


def potMaker(n1, a1, n2, a2, gridSize):
    """
	A function that creates a potential V(x) = a1 x^n1 + a2 x^n2 which can be added to the Hamilton matrix.
    :param gridSize: Length of the box.
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
    """
	 This plots the eigenfunctions and values for the harmonic oscillator. It also plots the energy levels of the 
    anharmonic oscillator with potential x^2 + x^4 and the difference in energy levels between the harmonic 
    oscillator and the anharmonic oscillator.
    """
    plotEigenstate(potMaker(2, 1/4, 0, 0, 101)[0], potMaker(2, 1/4, 0, 0, 101)[1], 1, 'Task 5.0')
    for i in np.arange(0, 1, 0.1):
        HVMatrix = potMaker(2, 1/4, 4, i, 101)  
        if i == 0: 
            energieMatrix = HVMatrix[0]
            """
            This is the matrix with only x^2 potential.
            """
        if i != 0:
            verschilMatrix = HVMatrix[0] - energieMatrix 
            """
            This is the difference between a1 x^2 and a1 x^2 + i x^4.
            """
            plt.figure('Task 5.2')
            plt.plot(verschilMatrix, 'o', label="a = " + str(i), markersize=2.5)
            plt.legend()
            plt.xlabel('n')
            plt.ylabel('$\Delta$E (A.U)')
            plt.title('Energieverschillen van V(x) = ax\u2074 + 1/4 x\u00B2 en V(x) = 1/4 x\u00B2')
        plt.figure('Task 5.1')
        plt.plot(HVMatrix[0], 'o', label="a = " + str(i), markersize=2.5)
        plt.legend()
        plt.xlabel('n')
        plt.ylabel('E (A.U)')
        plt.title('Energie met potentiaal V(x) = ax\u2074 + x\u00B2')
    plt.show()
	
task5()
