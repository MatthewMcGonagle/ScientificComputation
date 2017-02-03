import numpy as np

print( "Initializing Set Up" )
varbounds = np.array( [ [0., 1.]
                      , [0., 2.] ] )
nvarpoints = np.array([5, 3], dtype = int)
varpointsexcess = [ np.linspace(varbounds[i][0], varbounds[i][1], nvarpoints[i]+2) 
                        for i in range(len(nvarpoints))
                  ] 
varpoints = [varpointsexcess[i][1:-1:1] for i in range(len(varpointsexcess))]
dvar = [varpointsexcess[i][1] - varpointsexcess[i][0] for i in range(len(varpointsexcess))]


# Boundary values function

def f(x,y):
        return x**2 - y**2

# Source function 

def g(x,y):
        return 1.0+x

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

topbottom = np.fromfunction( lambda i, j : f(varpoints[0][j], varbounds[1][i]) 
                               , (2, nvarpoints[0])
                               , dtype = int
                               )


leftright = np.fromfunction( lambda i, j: f(varbounds[0][i], varpoints[1][j]) 
                                 , (2, nvarpoints[1])
                                 , dtype = int
                                 )

rhsterms = np.fromfunction( lambda i, j : g(varpoints[0][i], varpoints[1][j]) * scaling
                         , (nvarpoints[0], nvarpoints[1]) 
                         , dtype = int
                         )

for i in range(nvarpoints[0]) :
        rhsterms[i][0] -= topbottom[0][i] 
        rhsterms[i][-1] -= topbottom[1][i]

for i in range(nvarpoints[1]):
        rhsterms[0][i] -= leftright[0][i]
        rhsterms[-1][i] -= leftright[1][i]

# Reduce
# After completion, varpoints[xi][j][k] = coefficient for ( xi + 1, k ) when k <= j
#                                       = coefficient for ( xi, k ) when k > j
def reduceupycoord_notlastx(xi):
	for yj in range(nvarpoints[1]):

		# varpoints[xi][j][k] holds information on the coefficients in the equation for the
		# points (xi, j). Some coefficients are held implicitly and some coefficients are held
		# explicitly in varpoints[xi][j][[:].
		#
		# At the beginning of loop on yj:
		#
		# Should have varpoints[xi][j][k] = coeff of ( xi + 1 , k) for all j and all k < yj,
		#                                 = coeff of (xi, k) for all yj <= k.
		#
		# Also, implictly the coefficient for (xi + 1, yj) is diffcoeffsize[0].
		#
		# Finally, implicitly the equation for any point (xi , k) with k < yj 
		# has the coefficient 1.0 for itself, (xi, k).

		scaling = 1.0 / varpoints[xi][yj][yj]

		for yk in range(yj):
			varpoints[xi][yj][yk] *= scaling

		fixedcoeff = diffcoeffsize[0] * scaling 
		varpoints[xi][yj][yj] = fixedcoeff

		for yk in range(yj+1, nvarpoints[1]):
			varpoints[xi][yj][yk] *= scaling

		for yk in range(yj + 1, nvarpoints[1]):
			varpoints[xi][yk][yj] = -* fixedcoeff

