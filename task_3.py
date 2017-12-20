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
	barriePosition = int(data["gridSize"] / 2 - data["barrierWidth"] / 2)
	hamiltonMatrix = getHamiltonMatrix(data["gridSize"], data["t"])
	# Make the barrie, a negative barrier is a well
	for i in range(data["barrierWidth"]):
	    hamiltonMatrix[barriePosition + i][barriePosition + i] += data["barrierHeight"]

	# For periodic boundaries make the first and the last depend on eachother if it were modulus calculation.
	if data["periodicBoundaries"] == True:
		hamiltonMatrix[0][-1] = -data["t"]
		hamiltonMatrix[-1][0] = -data["t"]

	return calculateEigenstate(hamiltonMatrix, data["boxSize"])

scatterResults = []
for i in range(0, 400):
	eigens = calculateTask3(0, gridSize = 400, barrierHeight = i, barrierWidth = 10, periodicBoundaries = False, t = 1, boxSize = 1)
	for k in range(400):
		stop = False
		vector = eigens[1][:, k]
		for j in range(1, 190):
			# Use some epsilons here due to precisions
			if np.abs(vector[j]) > 10 ** -18 and np.abs(vector[-j]) > 10 ** -18:
				scatterResults.append([i, k - 1])
				stop = True
				break
		if stop:
			break
fig = plt.figure()
plt.plot([arr[0] for arr in scatterResults], [arr[1] for arr in scatterResults], "o", markersize = 2)
fig.canvas.set_window_title("Minimum scatter energy state distribution")
plt.xlabel("Potential height (a.u.)")
plt.ylabel("n")
plt.show()
eigens = calculateTask3(0, gridSize = 400, barrierHeight = 115, barrierWidth = 10, periodicBoundaries = False, t = 1, boxSize = 1)
plotEigenstate(eigens[0], eigens[1], 1, "Task 3.")
	
plt.show()
