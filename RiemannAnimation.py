import numpy as np
import matplotlib.pyplot as plt
import pylab

t1 = np.arange(0.0, 3.0, 0.0001)
numframes = 45
def f(t, n):
	result = 0
	for i in range(1, n+1):
		 result += np.sin(i**2*t)/i**2 
	return result 

y1 = f(t1, 1)

fig, ax = plt.subplots(1, 1)
for i in range(1, numframes+1):	
	y1 = f(t1, i)
	fig.suptitle('$f_n(t) = \sum_{k=1}^{'+str(i)+'}\sin(k^2 t)/k^2$', y = 0.99)
	ax.cla()
	ax.set_ylim([0,1.4])	
	ax.plot(t1, y1) 
	#ax.set_title('n = 1')
	
	#fig.tight_layout()
	fig.set_size_inches(6,6)
	if i<10:
		pylab.savefig('RiemannExample0'+str(i)+'.png', dpi=150)
	else:
		pylab.savefig('RiemannExample'+str(i)+'.png', dpi=150)
plt.show()
