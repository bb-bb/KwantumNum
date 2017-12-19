from quantummodule import *
import matplotlib.pyplot as plt

gridSize = 400
wellDepth = 100
wellWidth = 10
periodicBoundaries = True
t = 1
boxSize = 1
def performTask3(pictureNumber, **data):
	if data["gridSize"] < 3:
		print("The gridSize is too small: " + str(data["gridSize"]) + " in picture " + str(pictureNumber) + " results might be corrupted.")
	if data["wellWidth"] > data["gridSize"]:
		print("The wellWidth is bigger than the gridSize:" + str(data["wellWidth"]) + ">" + str(data["gridSize"]) + " in picture " + str(pictureNumber) + " addapting the wellWidth...")
		data["wellWidth"] = data["gridSize"]

	# Place well in middle of the chain (notice this does not matter when using periodic boundaries)
	# Do not care about odd and even numbers since the gridsize should be big enough
	wellPosition = int(data["gridSize"] / 2 - data["wellWidth"] / 2)
	hamiltonMatrix = getHamiltonMatrix(data["gridSize"], data["t"])
	# Make the well, a negative well is a bump
	for i in range(data["wellWidth"]):
	    hamiltonMatrix[wellPosition + i][wellPosition + i] -= data["wellDepth"]

	# For periodic boundaries make the first and the last depend on eachother if it were modulus calculation.
	if data["periodicBoundaries"]:
		hamiltonMatrix[0][-1] = -t
		hamiltonMatrix[-1][0] = -t

	values, vectors = calculateEigenstate(hamiltonMatrix, data["boxSize"])
	plotEigenstate(values, vectors, data["boxSize"], "Task 3." + str(pictureNumber))
	plt.show()

performTask3(0, gridSize = 400, wellDepth = 100, wellWidth = 10, periodicBoundaries = True, t = 1, boxSize = 1)
