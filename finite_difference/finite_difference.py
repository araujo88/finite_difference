import numpy as np


def Diff1(n, o=2):
    # Computes finite-difference matrices for the first derivative
    # Default order is set to second-order
    D = np.zeros((n, n))
    if o == 2:  # Second-order
        D[0, 0] = -1
        D[0, 1] = 1
        for i in range(1, (n-1)):
            D[i, i-1] = -0.5
            D[i, i] = 0
            D[i, i+1] = 0.5
        D[-1, -1] = D[0, 1]
        D[-1, -2] = D[0, 0]
        return D
    elif o == 4:  # Fourth-order
        D[0, 0] = -1
        D[0, 1] = 1
        D[1, 0] = -0.5
        D[1, 1] = 0
        D[1, 2] = 0.5
        for i in range(2, (n-2)):
            D[i, i-2] = 1/12
            D[i, i-1] = -2/3
            D[i, i] = 0
            D[i, i+1] = 2/3
            D[i, i+2] = -1/12
        D[-1, -1] = D[0, 1]
        D[-1, -2] = D[0, 0]
        D[-2, -1] = D[1, 2]
        D[-2, -2] = D[1, 1]
        D[-2, -3] = D[1, 0]
        return D
    elif o == 6:  # Sixth-order
        D[0, 0] = -1
        D[0, 1] = 1
        D[1, 0] = -0.5
        D[1, 1] = 0
        D[1, 2] = 0.5
        D[2, 0] = 1/12
        D[2, 1] = -2/3
        D[2, 2] = 0
        D[2, 3] = 2/3
        D[2, 4] = -1/12
        for i in range(3, (n-3)):
            D[i, i-3] = -1/60
            D[i, i-2] = 3/20
            D[i, i-1] = -3/4
            D[i, i] = 0
            D[i, i+1] = 3/4
            D[i, i+2] = -3/20
            D[i, i+3] = 1/60
        D[-1, -1] = D[0, 1]
        D[-1, -2] = D[0, 0]
        D[-2, -1] = D[1, 2]
        D[-2, -2] = D[1, 1]
        D[-2, -3] = D[1, 0]
        D[-3, -1] = D[2, 4]
        D[-3, -2] = D[2, 3]
        D[-3, -3] = D[2, 2]
        D[-3, -4] = D[2, 1]
        D[-3, -5] = D[2, 0]
        return D


def Diff2(n, o=2):
    # Computes finite-difference matrices for the second derivative
    # Default order is set to second-order
    D = np.zeros((n, n))
    if o == 2:  # Second-order
        D[0, 0] = 2  # Forward scheme (second-order)
        D[0, 1] = -5
        D[0, 2] = 4
        D[0, 3] = -1
        for i in range(1, (n-1)):
            D[i, i-1] = 1
            D[i, i] = -2
            D[i, i+1] = 1
        D[-1, -1] = D[0, 0]
        D[-1, -2] = D[0, 1]
        D[-1, -3] = D[0, 2]
        D[-1, -4] = D[0, 3]
        return D
    elif o == 4:  # Fourth-order
        D[0, 0] = 2  # Forward scheme (second-order)
        D[0, 1] = -5
        D[0, 2] = 4
        D[0, 3] = -1
        D[1, 0] = 1  # Central scheme (second-order)
        D[1, 1] = -2
        D[1, 2] = 1
        for i in range(2, (n-2)):
            D[i, i-2] = -1/12
            D[i, i-1] = 4/3
            D[i, i] = -5/2
            D[i, i+1] = 4/3
            D[i, i+2] = -1/12
        D[-1, -1] = D[0, 0]
        D[-1, -2] = D[0, 1]
        D[-1, -3] = D[0, 2]
        D[-1, -4] = D[0, 3]
        D[-2, -1] = D[1, 0]
        D[-2, -2] = D[1, 1]
        D[-2, -3] = D[1, 2]
        return D
    elif o == 6:  # Sixth-order
        D[0, 0] = 2  # Forward-scheme (second-order)
        D[0, 1] = -5
        D[0, 2] = 4
        D[0, 3] = -1
        D[1, 0] = 1  # Central-scheme (second-order)
        D[1, 1] = -2
        D[1, 2] = 1
        D[2, 0] = -1/12  # Central-scheme (fourth-order)
        D[2, 1] = 4/3
        D[2, 2] = -5/2
        D[2, 3] = 4/3
        D[2, 4] = -1/12
        for i in range(3, (n-3)):
            D[i, i-3] = 1/90
            D[i, i-2] = -3/20
            D[i, i-1] = 3/2
            D[i, i] = -49/18
            D[i, i+1] = 3/2
            D[i, i+2] = -3/20
            D[i, i+3] = 1/90
        D[-1, -1] = D[0, 0]
        D[-1, -2] = D[0, 1]
        D[-1, -3] = D[0, 2]
        D[-1, -4] = D[0, 3]
        D[-2, -1] = D[1, 0]
        D[-2, -2] = D[1, 1]
        D[-2, -3] = D[1, 2]
        D[-3, -1] = D[2, 0]
        D[-3, -2] = D[2, 1]
        D[-3, -3] = D[2, 2]
        D[-3, -4] = D[2, 3]
        D[-3, -5] = D[2, 4]
        return D


def PadeD1_4(f, h=1):

    N = int((np.size(f)))

    A = np.zeros((N, N))
    b = np.zeros(N)

    # Generate A matrix (f_i')

    # Boundary and adjacent points:

    A[0, 0] = 1
    A[0, 1] = 2

    A[N-1, N-1] = A[0, 0]
    A[N-1, N-2] = A[0, 1]

    # Interior points:

    for i in range(1, N-1):
        A[i, i-1] = 1
        A[i, i] = 4
        A[i, i+1] = 1

    # Generate b vector (f_i)

    # Boundary and adjacent points:

    b[0] = (-5/2)*f[0]+2*f[1]+0.5*f[2]
    b[N-1] = (5/2)*f[N-1]-2*f[N-2]-0.5*f[N-3]

    # Interior points:

    for i in range(1, N-2):
        b[i] = 3*(f[i+1]-f[i-1])/h

    # Compute the derivative:

    df = np.linalg.solve(A, b)

    return df


def LeleD1_6(f, h=1):
    # This is an implementation of a 6th order non-dissipative compact difference as shown in:
    # S. K. Lele. Compact finite difference schemes with spectral-like resolution.
    # Journal of Computational Physics, 103(1):16–42, November 1992.
    # While functional, this implementation is designed for demonstration
    # and so clarity was pursued at the expense of a more efficient implementation
    # Using a tridiagonal solver should allow for a much more efficient solution to the problem.
    # What we have is a system of the form Ax'=Bx=b
    N = int(np.sqrt(np.size(f)))
    # Initialize variable
    A = np.zeros((N, N))
    B = np.zeros((N, N))
    # Populate the A-matrix left boundary
    A[0, 0] = 1
    A[0, 1] = 2
    # 2nd point from left
    A[1, 0] = 1/4
    A[1, 1] = 1
    A[1, 2] = 1/4
    # All interior points
    for kk in range(2, N-2):
        A[kk, kk-1] = 1/3
        A[kk, kk] = 1
        A[kk, kk+1] = 1/3
    # 2nd point from right
    A[N-2, N-3] = 1/4
    A[N-2, N-2] = 1
    A[N-2, N-1] = 1/4
    # right boundary
    A[N-1, N-2] = 2
    A[N-1, N-1] = 1
    # Populate the B-matrix
    alf = 7/9/h
    bet = 1/36/h
    # left boundary
    B[0, 0] = -5/2/h
    B[0, 1] = 2/h
    B[0, 2] = 1/2/h
    # 2nd point from left
    B[1, 0] = -3/4/h
    B[1, 2] = 3/4/h
    # interior points
    for kk in range(2, N-2):
        B[kk, kk-2] = -bet
        B[kk, kk-1] = -alf
        B[kk, kk+1] = alf
        B[kk, kk+2] = bet
    # 2nd point from right
    B[N-2, N-3] = -3/4/h
    B[N-2, N-1] = 3/4/h
    # right boundary
    B[N-1, N-1] = 5/2/h
    B[N-1, N-2] = -2/h
    B[N-1, N-3] = -1/2/h
    b = np.matmul(B, f)
    outarray = np.linalg.solve(A, b)
    return outarray
