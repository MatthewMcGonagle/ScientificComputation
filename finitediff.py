import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

N = 50 # Mesh size, should be greater than 2
xa = -0.5 # x lower bound
xb = 0.5 # x upper bound
ya = -0.5 # y lower bound
yb = 0.5 # y upper bound
meshdeltax = (xb - xa) / (N-1)
meshdeltay = (yb - ya) / (N-1)

# Make mesh for x, y coordinates

x0 = np.arange(N)
y0 = np.arange(N)
ymesh, xmesh = np.meshgrid(x0, y0)

# Function g(x,y) for boundary values
def g(x,y):
	return x**5 - 10.0*x**3*y**2+5.0*x*y**4 

# Function f(x,y) for source term
def f(x,y):
	return 0.0

# coordinate functions of array indices
def xcoord(xi):
	return xa + meshdeltax * xi 

def ycoord(yi):
	return yb - meshdeltay * yi 

# Initialize Mesh for u(x,y). Set the boundary conditions.
umesh = [ [0 for cols in range(N)] for rows in range(N)]
for i in range(N):
	umesh[i][0] = g(xcoord(i), ya)
	umesh[i][N-1] = g(xcoord(i), yb)
	umesh[0][i] = g(xa, ycoord(i))
	umesh[N-1][i] = g(xb, ycoord(i))

# Initialize Mesh for f(x, y) values
fmesh = [ [meshdeltax**2*f(xcoord(xi), ycoord(yi)) for yi in range(N)] for xi in range(N)]

# Initialize Mesh for Equation Coefficients. 2D array of 1D array.
# Really only need to keep track of coefficients for at most 3 cols at a time.
# equation[i][j] = Linear Equation Coefficients for xij 
# [xij, x(i+1)j, xi(j+1), x(i-1)j, xi(j-1)]
# There are no linear equations for the boundary, so we have (N-2)x(N-2) equations.

eqvert = [[[0.0 for eqi in range(N)] for yi in range(N)] for xi in range(N)]
eqright = [[[0.0 for eqi in range(N)] for yi in range(N)] for xi in range(N)]
eqleftboundary = [[[0.0 for eqi in range(N)] for yi in range(N)] for xi in range(N)]
eqbottomboundary = [[[0.0 for eqi in range(N)] for yi in range(N)] for xi in range(N)]
eqtopboundary = [[[0.0 for eqi in range(N)] for yi in range(N)] for xi in range(N)]

# Properly initialize some values in equation arrays

for xi in range(1,N-1):
	for yi in range(1,N-1):
		eqright[xi][yi][yi] = 1.0 
		eqleftboundary[1][yi][yi] = 1.0
		if yi>1: # Leave bottom boundary point empty	
			eqvert[xi][yi][yi-1] = 1.0
		eqvert[xi][yi][yi] = -4.0
		if yi < N-1: # Leave top boundary point empty	
			eqvert[xi][yi][yi+1] = 1.0
	eqtopboundary[xi][N-2][xi] = 1.0
	eqbottomboundary[xi][1][xi] = 1.0

for yi in range(1,N-1):
	eqleftboundary[1][yi][yi] = 1.0

# Do Reductions
# For continuous memory considerations, we do reduction by columns

for xi in range(1,N-1):
	for yi in range(1,N-2):

		# Normalize with respect to coefficient for xij

		scaling = 1.0/eqvert[xi][yi][yi]
		for eqi in range(yi, N):
			eqvert[xi][yi][eqi] *= scaling
		for eqi in range(N):
			eqright[xi][yi][eqi] *= scaling 
			eqtopboundary[xi][yi][eqi] *= scaling 
			eqbottomboundary[xi][yi][eqi] *= scaling
			eqleftboundary[xi][yi][eqi] *=scaling
		fmesh[xi][yi] *= scaling 

		# ! Need to use xij to reduce all xik for k>j, not just j+1 !
		# Reduce the equation for xi(j+1)

		for eqyi in range(yi+1, N-1):	
			scaling = eqvert[xi][eqyi][yi]	
			for eqi in range(yi, N):
				eqvert[xi][eqyi][eqi] -= eqvert[xi][yi][eqi]*scaling
		
			for eqi in range(N):
				eqright[xi][eqyi][eqi] -= eqright[xi][yi][eqi]*scaling 
				eqleftboundary[xi][eqyi][eqi] -= eqleftboundary[xi][yi][eqi]*scaling
				eqbottomboundary[xi][eqyi][eqi] -= eqbottomboundary[xi][yi][eqi]*scaling
				eqtopboundary[xi][eqyi][eqi] -= eqtopboundary[xi][yi][eqi]*scaling
			fmesh[xi][eqyi] -= fmesh[xi][yi]*scaling

	for yi in reversed(range(2,N-1)):

		# Normalize the equation for xij
		
		scaling = 1.0/eqvert[xi][yi][yi]
		for eqi in range(yi, N):
			eqvert[xi][yi][eqi] *= scaling
		for eqi in range(N):
			eqright[xi][yi][eqi] *= scaling 
			eqtopboundary[xi][yi][eqi] *= scaling 
			eqbottomboundary[xi][yi][eqi] *= scaling
			eqleftboundary[xi][yi][eqi] *= scaling
		fmesh[xi][yi] *= scaling

		# Reduce the equation for xik for k<j
		
		for eqyi in range(1, yi):

			scaling = eqvert[xi][eqyi][yi]

			eqvert[xi][eqyi][yi] -= eqvert[xi][yi][yi]*scaling	
			for eqi in range(N):
				eqright[xi][eqyi][eqi] -= eqright[xi][yi][eqi]*scaling
				eqleftboundary[xi][eqyi][eqi] -= eqleftboundary[xi][yi][eqi]*scaling
				eqbottomboundary[xi][eqyi][eqi] -= eqbottomboundary[xi][yi][eqi]*scaling
				eqtopboundary[xi][eqyi][eqi] -= eqtopboundary[xi][yi][eqi]*scaling

			fmesh[xi][eqyi] -= fmesh[xi][yi]*scaling

	for yi in range(1,N-1):

		# Reduce the equation for x(i+1)j
		scaling = 1.0

		# Don't need to explicitly reduce coefficient for xij in x(i+1)j. We don't need to keep track of it.

		for eqi in range(1,N-1):
			eqvert[xi+1][yi][eqi] -= eqright[xi][yi][eqi]*scaling
		for eqi in range(N):
			eqtopboundary[xi+1][yi][eqi] -= eqtopboundary[xi][yi][eqi]*scaling
			eqbottomboundary[xi+1][yi][eqi] -= eqbottomboundary[xi][yi][eqi]*scaling
			eqleftboundary[xi+1][yi][eqi] -= eqleftboundary[xi][yi][eqi]*scaling
		fmesh[xi+1][yi] -= fmesh[xi][yi]*scaling

# Now compute umesh[xi][xj]
for xi in reversed(range(1,N-1)):
	for yi in range(1,N-1):
		umesh[xi][yi] = fmesh[xi][yi]
		for eqi in range(N):
			umesh[xi][yi] -= eqright[xi][yi][eqi]*umesh[xi+1][eqi]
			umesh[xi][yi] -= eqbottomboundary[xi][yi][eqi]*umesh[eqi][0]
			umesh[xi][yi] -= eqleftboundary[xi][yi][eqi]*umesh[0][eqi]
			umesh[xi][yi] -= eqtopboundary[xi][yi][eqi]*umesh[eqi][N-1]

# test against known function
error = 0.0
for xi in range(N):
	for yi in range(N):
		error += abs(umesh[xi][yi] - g(xcoord(xi), ycoord(yi)))
error /= N**2

print(error)

# Plot the computed function

for xi in range(N):
	for yi in range(N):
		dummy = 1
	#	xmesh[xi][yi] = xcoord(xi)

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_wireframe(xmesh, ymesh, umesh)
plt.show()

