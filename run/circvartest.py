import numpy as np

def CircVar(v,u):
	# Circular variance
	c  = np.hypot(v,u)
	vv = np.sum(v/c)/v.size   # equivalent to (SUM cos(*theta*) ) / n 
	uu = np.sum(u/c)/u.size   # equivalent to (SUM sin(*theta*) ) / n 
	r  = np.hypot(uu,vv)	
	dispersion = 1. - r  # should be between zero and one 
#	return dispersion
	return np.sqrt(-2*np.log(r))


#angles = [110.0,60.0,90.,120.0,150., 20., 20., 210., 240., 270., 300., 330., 360., 0.]
angles  = [0, 30, 90, 30, 0]
angs   = np.sin(angles) 
angc   = np.cos(angles)

a = sum(np.sin(angles))/len(angles)
b = sum(np.cos(angles))/len(angles)

r = np.sqrt(-2*np.log(np.hypot(a,b)))


print 'check: ', r

print 'CircVar: ', CircVar(angs,angc)

