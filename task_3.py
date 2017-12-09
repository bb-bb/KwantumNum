from quantummodule import *
import matplotlib.pyplot as plt

gridSize = 400
wellDepth = 1
wellWidth = 10
t = 1
boxSize = 1

wellPosition = int(gridSize / 2 - wellWidth / 2)
hamiltonMatrix = getHamiltonMatrix(gridSize, t)
for i in range(wellWidth):
    hamiltonMatrix[wellPosition + i][wellPosition + i] -= wellDepth

values, vectors = calculateEigenstate(hamiltonMatrix, boxSize)
plotEigenstate(values, vectors, boxSize, 'Task 3.0')
plt.show()
