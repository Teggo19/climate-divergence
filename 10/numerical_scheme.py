import numpy as np
import matplotlib.pyplot as plt

def create_diagonal(N, dx, x, D, B_out):
    """
    Create the diagonal of the matrix for the numerical scheme.
    """
    diag = np.empty(N+1)
    diag[0] = B_out
    diag[-1] = 0
    diag[1:-1] = D * (2 - (x[1:-1]-0.5*dx)**2 - (x[1:-1]+0.5*dx)**2) / dx**2 + B_out
    return diag


def create_lower(N, dx, x, D):
    """
    Create the lower diagonal of the matrix for the numerical scheme.
    """
    lower = np.empty(N)
    lower[:-1] = -D * (1 - (x[1:-1]-0.5*dx)**2) / dx**2
    lower[-1] = D * (2 - x[-1]**2 - x[-2]**2) / dx**2
    return lower


def create_upper(N, dx, x, D):
    """
    Create the upper diagonal of the matrix for the numerical scheme.
    """
    upper = np.empty(N)
    upper[0] = D * (1 - x[1]**2) / dx**2
    upper[1:] = -D * (1 - (x[1:-1]+0.5*dx)**2) / dx**2
    return upper


def create_matrix(N, dx, x, D, B_out):
    """
    Create a matrix for the numerical scheme.
    """
    lower = create_lower(N, dx, x, D)
    upper = create_upper(N, dx, x, D)
    diag = create_diagonal(N, dx, x, D, B_out)
    A = np.diag(diag) + np.diag(lower, -1) + np.diag(upper, 1)
    A[0,2] = -D * (1 - x[1]**2) / dx**2
    A[-1,-3] = -D * (1 - x[-2]**2) / dx**2

    return A


def create_rhs(N, dx, x, D, A_out, B_out, Q, S, a, T_s):
    """
    Create the right-hand side of the numerical scheme.
    """
    rhs = -A_out + Q * S(x) * a
    # rhs = A_out + Q * S(x) * a # Flipping the sign in front of A_out giver reasonable results, but must be wrong
    rhs[-1] += (D * (1 - x[-1]**2) / dx**2 - B_out) * T_s
    return rhs



def main():
    N = 50
    x_s = 0.6
    x, dx = np.linspace(0, x_s, N+1, retstep=True)
    D = 0.3
    a = 0.68
    A_out = 201.4
    B_out = 1.45
    T_s = 273.15
    Q = 340
    def S(x, S_2=-0.477):
        return 1 + 0.5*S_2 * (3*x**2 - 1)

    f = create_rhs(N, dx, x, D, A_out, B_out, Q, S, a, T_s)
    A = create_matrix(N, dx, x, D, B_out)
    T = np.linalg.solve(A, f)

    print(T)
    
    plt.plot(x, T)
    plt.show()

if __name__ == '__main__':
    main()
