import numpy
import matplotlib.pyplot as pyplot

epsilon = 1.0
lamba = 0.05
delta = 0.002
omega = -1.0 

f = (lamba**2 + (omega + 2.0*epsilon)**2/4.0)
t = range(1,101)
T = [2.0*numpy.pi * i/(100.0*f) for i in t]
K = [i * numpy.sqrt(f) for i in T]
coeff = lamba**2 * numpy.sin( K )**2/f;
print len(coeff)
print len(T)
fi = numpy.loadtxt("numeric")
pyplot.plot( T, coeff )
#pyplot.plot( T, fi[:,0] )
pyplot.plot( T, fi[:,1] )

pyplot.savefig( "exact" )
