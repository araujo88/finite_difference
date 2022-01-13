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

    for i in range(1, N-1):
        b[i] = 3*(f[i+1]-f[i-1])

    b = b/h

    # Compute the derivative:

    df = np.linalg.solve(A, b)

    return df


def PadeD2_4(f, h=1):

    N = int((np.size(f)))

    A = np.zeros((N, N))
    b = np.zeros(N)

    # Generate A matrix (f_i')

    # Boundary and adjacent points:
    A[0, 0] = 1

    A[N-1, N - 1] = A[0, 0]

    # Interior points:

    for i in range(1, N-1):
        A[i, i-1] = 1/10
        A[i, i] = 1
        A[i, i+1] = 1/10

    # Generate b vector (f_i)

    # Boundary and adjacent points:

    b[0] = (2)*f[0]+(-5)*f[1]+(4)*f[2]-f[3]
    b[N-1] = (2)*f[-1]+(-5)*f[-2]+(4)*f[-3]-f[-4]

    # Interior points:

    for i in range(1, N-1):
        b[i] = (6/5)*(f[i+1]-2*f[i]+f[i-1])

    b = b/(h**2)

    # Compute the derivative:

    df = np.linalg.solve(A, b)

    return df


def LeleD1_6(f, h=1):

    N = int((np.size(f)))

    A = np.zeros((N, N))
    b = np.zeros(N)

    # Generate A matrix (f_i')

    # Boundary and adjacent points:

    A[0, 0] = 1
    A[0, 1] = 2
    A[1, 0] = 1/4
    A[1, 1] = 1
    A[1, 2] = 1/4

    A[N-1, N-1] = A[0, 0]
    A[N-1, N-2] = A[0, 1]
    A[N-2, N-1] = A[1, 0]
    A[N-2, N-2] = A[1, 1]
    A[N-2, N-3] = A[1, 2]

    # Interior points:

    for i in range(2, N-2):
        A[i, i-1] = 1/3
        A[i, i] = 1
        A[i, i+1] = 1/3

    # Generate b vector (f_i)

    # Boundary and adjacent points:

    b[0] = (-5*f[0]+4*f[1]+f[2])/(2*h)
    b[1] = (3/2)*(f[2]-f[0])/(2*h)
    
    b[N-1] = -(-5*f[N-1]+4*f[N-2]+f[N-3])/(2*h)
    b[N-2] = -(3/2)*(f[-3]-f[-1])/(2*h)

    # Interior points:

    for i in range(2, N-2):
        b[i] = (14/9)*(f[i+1]-f[i-1])/(2*h) + (1/9)*(f[i+2]-f[i-2])/(4*h)

    # Compute the derivative:

    df = np.linalg.solve(A, b)

    return df


def LeleD2_6(f, h=1):

    N = int((np.size(f)))

    A = np.zeros((N, N))
    b = np.zeros(N)

    # Generate A matrix (f_i')

    # Boundary and adjacent points:
    A[0, 0] = 1
    A[0, 1] = 11
    A[1, 0] = 1/10
    A[1, 1] = 1
    A[1, 2] = 1/10
    A[-1, -1] = A[0, 0]
    A[-1, -2] = A[0, 1]
    A[-2, -1] = A[1, 0]
    A[-2, -2] = A[1, 1]
    A[-2, -3] = A[1, 2]

    # Interior points:

    for i in range(2, N-2):
        A[i, i-1] = 2/11
        A[i, i] = 1
        A[i, i+1] = 2/11

    # Generate b vector (f_i)

    # Boundary and adjacent points:

    b[0] = (13*f[0]-27*f[1]+15*f[2]-f[3])/(h**2)
    b[1] = (6/5)*(f[0]-2*f[1]+f[2])/(h**2)
    b[-1] = (13*f[-1]-27*f[-2]+15*f[-3]-f[-4])/(h**2)
    b[-2] = (6/5)*(f[-1]-2*f[-2]+f[-3])/(h**2)

    # Interior points:

    for i in range(2, N-2):
        b[i] = (3/11)*(f[i+2]-2*f[i]+f[i-2])/(4*h**2) + \
            (12/11)*(f[i+1]-2*f[i]+f[i-1])/(h**2)

    # Compute the derivative:

    df = np.linalg.solve(A, b)

    return df
