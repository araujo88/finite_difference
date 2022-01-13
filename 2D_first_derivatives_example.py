import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from finite_difference.finite_difference import Diff1, PadeD1_4, LeleD1_6, Diff2, PadeD2_4, LeleD2_6


def contour(X, Y, Z):
    plt.contourf(X, Y, Z, 40, cmap='gist_rainbow_r')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.colorbar()
    plt.show()


nx = 50  # Number of points in x direction
ny = 50  # Number of points in y direction
x = np.linspace(0, 2*np.pi, nx, endpoint=True)  # X-points
y = np.linspace(0, 2*np.pi, ny, endpoint=True)  # Y-points
dx = x[1]-x[0]  # delta_x
dy = y[1]-y[0]  # delta_y
X, Y = np.meshgrid(x, y, indexing='ij')  # generate mesh
Z = np.sin(X)*np.sin(Y)  # function to be differentiated
dzdx = np.cos(X)*np.sin(Y)  # partial derivative (x)
dzdy = np.sin(X)*np.cos(Y)  # partial derivative (y)
# contour(X,Y,Z)

# First derivatives

I = np.eye(nx, ny)  # identity matrix

d_x2 = Diff1(nx, 2)/dx  # second-order diff coeffs matrix (x)
d_x4 = Diff1(nx, 4)/dx  # fourth-order diff coeffs matrix (x)
d_x6 = Diff1(nx, 6)/dx  # sixth-order diff coeffs matrix (x)

d_y2 = Diff1(nx, 2)/dy  # second-order diff coeffs matrix (y)
d_y4 = Diff1(nx, 4)/dy  # fourth-order diff coeffs matrix (y)
d_y6 = Diff1(nx, 6)/dy  # sixth-order diff coeffs matrix (y)

# Approximation of partial derivatives (x)
DX2 = sparse.kron(d_x2, I)
dzdx2 = DX2 @ np.reshape(Z, (nx*ny, 1))
dzdx2 = np.reshape(dzdx2, (nx, ny))

DX4 = sparse.kron(d_x4, I)
dzdx4 = DX4 @ np.reshape(Z, (nx*ny, 1))
dzdx4 = np.reshape(dzdx4, (nx, ny))

DX6 = sparse.kron(d_x6, I)
dzdx6 = DX6 @ np.reshape(Z, (nx*ny, 1))
dzdx6 = np.reshape(dzdx6, (nx, ny))

# Approximation of partial derivatives (x)
DY2 = sparse.kron(I, d_y2)
dzdy2 = DY2 @ np.reshape(Z, (nx*ny, 1))
dzdy2 = np.reshape(dzdy2, (nx, ny))

DY4 = sparse.kron(I, d_y4)
dzdy4 = DY4 @ np.reshape(Z, (nx*ny, 1))
dzdy4 = np.reshape(dzdy4, (nx, ny))

DY6 = sparse.kron(I, d_y6)
dzdy6 = DY6 @ np.reshape(Z, (nx*ny, 1))
dzdy6 = np.reshape(dzdy6, (nx, ny))

dzdx_pade = np.transpose(np.reshape(PadeD1_4(np.reshape(Z, (nx*ny, 1)), dx), (nx, ny)))
dzdy_pade = np.transpose(np.reshape(PadeD1_4(np.reshape(Z, (nx*ny, 1)), dy), (nx, ny)))

dzdx_lele = np.transpose(np.reshape(LeleD1_6(np.reshape(Z, (nx*ny, 1)), dx), (nx, ny)))
dzdy_lele = np.transpose(np.reshape(LeleD1_6(np.reshape(Z, (nx*ny, 1)), dy), (nx, ny)))

# Plots absolute error (x) - log scale
plt.title('Absolute error (log scale) - x-derivative (2nd order)')
contour(X, Y, np.log10(abs(dzdx2-dzdx)))
plt.title('Absolute error (log scale) - x-derivative (4th order)')
contour(X, Y, np.log10(abs(dzdx4-dzdx)))
plt.title('Absolute error (log scale) - x-derivative (6th order)')
contour(X, Y, np.log10(abs(dzdx6-dzdx)))
plt.title('Absolute error (log scale) - x-derivative (Padé 4th order)')
contour(X, Y, np.log10(abs(dzdx_pade-dzdx)))
plt.title('Absolute error (log scale) - x-derivative (Lele 6th order)')
contour(X, Y, np.log10(abs(dzdx_lele-dzdx)))

# Plots absolute error (y) - log scale
plt.title('Absolute error (log scale) - y-derivative (2nd order)')
contour(X, Y, np.log10(abs(dzdy2-dzdy)))
plt.title('Absolute error (log scale) - y-derivative (4th order)')
contour(X, Y, np.log10(abs(dzdy4-dzdy)))
plt.title('Absolute error (log scale) - y-derivative (6th order)')
contour(X, Y, np.log10(abs(dzdy6-dzdy)))
plt.title('Absolute error (log scale) - y-derivative (Padé 4th order)')
contour(X, Y, np.log10(abs(dzdy_pade-dzdx)))
plt.title('Absolute error (log scale) - y-derivative (Lele 6th order)')
contour(X, Y, np.log10(abs(dzdy_lele-dzdx)))
