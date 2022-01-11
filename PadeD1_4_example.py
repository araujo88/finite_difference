from finite_difference.finite_difference import PadeD1_4
import numpy as np
import matplotlib.pyplot as plt


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

Z = np.reshape(Z, (nx*ny, 1))

# First derivatives
# Approximation of partial derivative (x)
dzdx2 = PadeD1_4(Z, dx)

# Approximation of partial derivative (y)
dzdy2 = PadeD1_4(Z, dy)

# Plots absolute error (x) - log scale
plt.title('Abs error (log scale) - x-derivative (Padé 4th-order compact scheme)')
contour(X, Y, np.log10(abs(np.transpose(np.reshape(dzdx2, (nx, ny)))-dzdx)))

# Plots absolute error (y) - log scale
plt.title('Abs error (log scale) - y-derivative (Padé 4th-order compact scheme)')
contour(X, Y, np.log10(abs(np.reshape(dzdy2, (nx, ny))-dzdy)))
