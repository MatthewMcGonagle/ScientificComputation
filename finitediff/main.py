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

varcoeff = np.zeros( (nvarpoints[0], nvarpoints[1], 3, nvarpoints[1]) )
for i in range(nvarpoints[0]):
        for j in range(nvarpoints[1]):
                if i > 0:
                        varcoeff[i][j][0][j] = dvar[0]**-2
                if i < nvarpoints[0] - 1 :
                        varcoeff[i][j][2][j] = dvar[0]**-2
                if j > 0:
                        varcoeff[i][j][1][j-1] = dvar[1]**-2
                if j < nvarpoints[1] - 1 :
                        varcoeff[i][j][1][j+1] = dvar[1]**-2
                varcoeff[i][j][1][j] = dvar[0]**-2 + dvar[1]**-2
                varcoeff[i][j][1][j] *= -2.0

topbottom = np.fromfunction( lambda i, j : f(varpoints[0][j], varbounds[1][i])
                               , (2, nvarpoints[0])
                               , dtype = int
                               )

leftright = np.fromfunction( lambda i, j: f(varbounds[0][i], varpoints[1][j])
                                 , (2, nvarpoints[1])
                                 , dtype = int
                                 )

rhsterms = np.fromfunction( lambda i, j : g(varpoints[0][i], varpoints[1][j])
                         , (nvarpoints[0], nvarpoints[1]) 
                         , dtype = int
                         )

for i in range(nvarpoints[0]) :
        rhsterms[i][0] -= topbottom[0][i] * dvar[1]**-2
        rhsterms[i][-1] -= topbottom[1][i] * dvar[1]**-2

for i in range(nvarpoints[1]):
        rhsterms[0][i] -= leftright[0][i] * dvar[0]**-2
        rhsterms[-1][i] -= leftright[1][i] * dvar[0]**-2
