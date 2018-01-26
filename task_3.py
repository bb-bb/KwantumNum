"""
Author Bouke Jansen
Task 3 Scattering from a potential barrier
"""
from quantummodule import *
import matplotlib.pyplot as plt

def calculateTask3(pictureNumber, **data):
	if data["gridSize"] < 3 and periodicBoundaries:
		print("The gridSize is too small for periodic boundaries: " + str(data["gridSize"]) + " results might be corrupted.")
	if data["barrierWidth"] > data["gridSize"]:
		print("The barrierWidth is bigger than the gridSize:" + str(data["barierWidth"]) + ">" + str(data["gridSize"]) + " addapting the barrierWidth...")
		data["barrierWidth"] = data["gridSize"]

	# Place barier in middle of the chain (notice this does not matter when using periodic boundaries)
	# Do not care about odd and even numbers since the gridsize should be big enough
	barrierPosition = int(data["gridSize"] / 2 - data["barrierWidth"] / 2)
	hamiltonMatrix = getHamiltonMatrix(data["gridSize"], data["t"])
	# Make the barrier, a negative barrier is a well
	for i in range(data["barrierWidth"]):
	    hamiltonMatrix[barrierPosition + i][barrierPosition + i] += data["barrierHeight"]

	# For periodic boundaries make the first and the last depend on eachother if it were modulus calculation.
	if data["periodicBoundaries"] == True:
		hamiltonMatrix[0][-1] = -data["t"]
		hamiltonMatrix[-1][0] = -data["t"]

	return calculateEigenstate(hamiltonMatrix, data["boxSize"])

eigens = calculateTask3(0, gridSize = 400, barrierHeight = 2, barrierWidth = 10, periodicBoundaries = False, t = 1, boxSize = 1)
plotEigenstate(eigens[0], eigens[1], 1, "Task 3.1")
plt.show()
eigens = calculateTask3(0, gridSize = 400, barrierHeight = 4, barrierWidth = 10, periodicBoundaries = False, t = 1, boxSize = 1)
plotEigenstate(eigens[0], eigens[1], 1, "Task 3.2")
plt.show()
eigens = calculateTask3(0, gridSize = 400, barrierHeight = 2, barrierWidth = 10, periodicBoundaries = True, t = 1, boxSize = 1)
plotEigenstate(eigens[0], eigens[1], 1, "Task 3.3")
plt.show()
	
