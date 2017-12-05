import numpy as np
import matplotlib.pyplot as plt

g_N =99
arr = []
for i in range(g_N):
    arr.append([0 for i in range(g_N)])
    for j in range(g_N): 
        if j == i:
            arr[i][j]=2
        elif (j == i+1 or j == i-1) and i+1 <= g_N and i -1 >= -1 :
            arr[i][j] = -1
            
g_Hamilton = np.matrix(arr)
eigenvalues = np.linalg.eig(g_Hamilton)
sortedEigenvalues = np.sort(eigenvalues[0])
plt.plot(sortedEigenvalues, '.')
plt.ylabel('Eigenvalue')
plt.xlabel('n')
plt.show()
