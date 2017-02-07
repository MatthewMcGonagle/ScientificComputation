from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt

print( "Initializing Set Up" )
varbounds = np.array( [ [-0.5, 0.5]
                      , [-0.5, 0.5] ] )
nvarpoints = np.array([10, 30], dtype = int)
varpointsexcess = [ np.linspace(varbounds[i][0], varbounds[i][1], nvarpoints[i]+2) 
                        for i in range(len(nvarpoints))
                  ] 
varpoints = [varpointsexcess[i][1:-1:1] for i in range(len(varpointsexcess))]
dvar = [varpointsexcess[i][1] - varpointsexcess[i][0] for i in range(len(varpointsexcess))]

print(dvar)

# Boundary values function

def f(x,y):
	return (x**2 - y**2)*(-np.log(x**2+y**2))**0.5
        # return x**3 + y**3 

# Source function 

def g0(x,y):
	return -np.log(x**2 + y**2)
def g1(x,y):
	return 2*x*g0(x,y)**-0.5 * -1 / (x**2+y**2) * 2*x 
def g2(x,y):
	# dx = (-np.log(x**2 + y**2)**-0.5 * -x / (x**2 + y**2)
	temp = g0(x,y)**-1.5 * x**2 / (x**2 + y**2)**2 - g0(x,y)**-0.5 / (x**2 + y**2) + g0(x,y)**-0.5 * 2*x**2 / (x**2 + y**2)**2 
	return (x**2 - y**2) * temp
def g(x,y):
	return g1(x,y) - g1(y,x) + g2(x,y) + g2(y,x) 
	# return 6*x + 6*y 

print( "Initializing Data Arrays" )

varcoeff = np.zeros( ( nvarpoints[0]
		     , nvarpoints[1]
		     , nvarpoints[1]
		     ) )
scaling = -1.0 / (dvar[0]**-2 + dvar[1]**-2) / 2.0
diffcoeffsize = np.array( [dvar[0]**-2, dvar[1]**-2] ) * scaling
for i in range(nvarpoints[0]):
        for j in range(nvarpoints[1]):
                if j > 0:
                        varcoeff[i][j][j-1] = diffcoeffsize[1] 
                if j < nvarpoints[1] - 1 :
                        varcoeff[i][j][j+1] = diffcoeffsize[1] 
                varcoeff[i][j][j] = 1.0 

topbottom = [ [ f(varpoints[0][j], varbounds[1][i])
	for j in range(nvarpoints[0]) ]
	for i in range(2) ]

leftright = [ [ f(varbounds[0][i], varpoints[1][j])
	for j in range(nvarpoints[1])]
	for i in range(2) ]

rhsterms = [ [g(varpoints[0][i], varpoints[1][j]) * scaling
	for j in range(nvarpoints[1])]
	for i in range(nvarpoints[0])]
rhsterms = np.array( rhsterms)

for i in range(nvarpoints[0]) :
        rhsterms[i][0] -= topbottom[0][i] * diffcoeffsize[1] 
        rhsterms[i][-1] -= topbottom[1][i] * diffcoeffsize[1] 

for i in range(nvarpoints[1]):
        rhsterms[0][i] -= leftright[0][i] * diffcoeffsize[0] 
        rhsterms[-1][i] -= leftright[1][i] * diffcoeffsize[0] 

# Reduce
# After completion, varpoints[xi][j][k] = coefficient for ( xi + 1, k ) when k <= j
#                                       = coefficient for ( xi, k ) when k > j
def reduceupycoord(xi):
	for thiseqy in range(nvarpoints[1]):

		# varpoints[xi][j][k] holds information on the coefficients in the equation for the
		# points (xi, j). Some coefficients are held implicitly and some coefficients are held
		# explicitly in varpoints[xi][j][[:].
		#
		# At the beginning of loop on thiseqy:
		#
		# Should have varpoints[xi][j][k] = coeff of ( xi + 1 , k) for all j and all k < thiseqy,
		#                                 = coeff of (xi, k) for all k >= thiseqy.
		#
		# Also, implictly the coefficient for (xi + 1, thiseqy) is diffcoeffsize[0].
		#
		# Finally, implicitly the equation for any point (xi , k) with k < thiseqy 
		# has the coefficient 1.0 for itself, (xi, k).

		scaling = 1.0 / varcoeff[xi][thiseqy][thiseqy]
		rightcoeff = scaling * diffcoeffsize[0]

		for yk in range(nvarpoints[1]):
			varcoeff[xi][thiseqy][yk] *= scaling
		rhsterms[xi][thiseqy] *= scaling

		# Store the coeff of ( xi + 1, thiseqy ) at varcoeff[xi][thiseqy][thiseqy]
		varcoeff[xi][thiseqy][thiseqy] = scaling * diffcoeffsize[0] 

		# Reduce for 0 to thiseqy - 1
		for toreduce in range(0, thiseqy):
			scaling = -varcoeff[xi][toreduce][thiseqy]
			for yl in range(nvarpoints[1]):
				varcoeff[xi][toreduce][yl] += scaling * varcoeff[xi][thiseqy][yl]
			varcoeff[xi][toreduce][thiseqy] = scaling * varcoeff[xi][thiseqy][thiseqy]
			rhsterms[xi][toreduce] += scaling * rhsterms[xi][thiseqy]

		for toreduce in range(thiseqy + 1, nvarpoints[1]):
			scaling = -varcoeff[xi][toreduce][thiseqy]	
			for yl in range(nvarpoints[1]):
				varcoeff[xi][toreduce][yl] += scaling * varcoeff[xi][thiseqy][yl]
			varcoeff[xi][toreduce][thiseqy] = scaling * varcoeff[xi][thiseqy][thiseqy]
			rhsterms[xi][toreduce] += scaling * rhsterms[xi][thiseqy]

def reducenextcol(xi):
	for thiseqy in range(nvarpoints[1]):
		scaling = -diffcoeffsize[0]
		for yk in range(nvarpoints[1]):
			varcoeff[xi+1][thiseqy][yk] += scaling * varcoeff[xi][thiseqy][yk]
		rhsterms[xi+1][thiseqy] += scaling * rhsterms[xi][thiseqy]
def backsolve(xi):
	for thiseqy in range(nvarpoints[1]):
		for yk in range(nvarpoints[1]):
			rhsterms[xi][thiseqy] -= varcoeff[xi][thiseqy][yk] * rhsterms[xi+1][yk]

for xi in range(nvarpoints[0]-1):
	reduceupycoord(xi)
	reducenextcol(xi)
reduceupycoord(nvarpoints[0]-1)

for xi in range(nvarpoints[0]-2, -1, -1):
	backsolve(xi)

Y, X = np.meshgrid(varpoints[1], varpoints[0])

print (Y.shape)
print (X.shape)
print (rhsterms.shape)

fig = plt.figure()
ax = fig.gca(projection = '3d')
surf = ax.plot_wireframe(X, Y, rhsterms)

truefunc = [[ f(varpoints[0][i], varpoints[1][j])
	for j in range(nvarpoints[1])]
	for i in range(nvarpoints[0])]
truefunc = np.array(truefunc)

surf2 = ax.plot_wireframe(X, Y, truefunc, color='red')
plt.show()

error = rhsterms - truefunc
print(np.amax(error))
print(np.amin(error))
