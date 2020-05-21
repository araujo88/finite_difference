import numpy as np
import matplotlib.pyplot as plt
from derivatives import Diff1, Diff2


N=50 # Number of points
x=np.linspace(0,1,N) # X=points
dx=x[1]-x[0] # delta_x

y=np.exp(x**2) # function to be differentiated
dydx=2*np.exp(x**2)*x # first derivative
d2ydx2=2*np.exp(x**2)*(2*x**2+1) # second derivative

# First derivatives

d2=Diff1(N,2)/dx # second-order diff matrix
d4=Diff1(N,4)/dx # fourth-order diff matrix
d6=Diff1(N,6)/dx # sixth-order diff matrix

dydx_2=d2 @ y # second-order 1st. derivative approx.
dydx_4=d4 @ y # fourth-order 1st. derivative approx.
dydx_6=d6 @ y # sixth-order 1st. derivative approx.

# Computes relative error for first-derivatives
plt.semilogy(x, abs(dydx_2-dydx)/dydx, '-rs', label='Second order')
plt.semilogy(x, abs(dydx_4-dydx)/dydx, '--bo', label='Fourth order')
plt.semilogy( x, abs(dydx_6-dydx)/dydx, 'g*-', label='Sixth order')
plt.xlabel('x')
plt.ylabel('Relative Error')
plt.title('Relative error for the first derivative')
plt.aspect='equal'
plt.legend() 
plt.grid()
plt.show()

# Second derivatives

dd2_2=Diff2(N,2)/(dx**2) # second-order diff matrix
dd2_4=Diff2(N,4)/(dx**2) # fourth-order diff matrix
dd2_6=Diff2(N,6)/(dx**2) # sixth-order diff matrix

d2ydx2_2=dd2_2 @ y # second-order 2nd. derivative approx.
d2ydx2_4=dd2_4 @ y # fourth-order 2nd. derivative approx.
d2ydx2_6=dd2_6 @ y # sixth-order 2nd. derivative approx.

# Computes relative error for second-derivatives
plt.semilogy(x, abs(d2ydx2_2-d2ydx2)/d2ydx2, '-rs', label='Second order')
plt.semilogy(x, abs(d2ydx2_4-d2ydx2)/d2ydx2, '--bo', label='Fourth order')
plt.semilogy( x, abs(d2ydx2_6-d2ydx2)/d2ydx2, 'g*-', label='Sixth order')
plt.xlabel('x')
plt.ylabel('Relative Error')
plt.title('Relative error for the second derivative')
plt.aspect='equal'
plt.legend() 
plt.grid()
plt.show()
