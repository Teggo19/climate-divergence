import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim

def compute_diagonal(N, dx, x, D, B_out):
    """
    compute the diagonal of the matrix for the numerical scheme.
    """
    diag = np.zeros(N+1)
    diag[1:-1] = D * (2 - (x[1:-1]-0.5*dx)**2 - (x[1:-1]+0.5*dx)**2) / dx**2 + B_out
    return diag


def compute_lower(N, dx, x, D):
    """
    compute the lower diagonal of the matrix for the numerical scheme.
    """
    lower = np.zeros(N)
    lower[:-1] = -D * (1 - (x[1:-1]-0.5*dx)**2) / dx**2
    return lower


def compute_upper(N, dx, x, D):
    """
    compute the upper diagonal of the matrix for the numerical scheme.
    """
    upper = np.zeros(N)
    upper[1:] = -D * (1 - (x[1:-1]+0.5*dx)**2) / dx**2
    return upper


def compute_south_matrix(N, dx, x, D, B_out):
    """
    compute a matrix for the numerical scheme.
    """
    lower = compute_lower(N, dx, x, D)
    upper = compute_upper(N, dx, x, D)
    diag = compute_diagonal(N, dx, x, D, B_out)
    A = np.diag(diag) + np.diag(lower, -1) + np.diag(upper, 1)
    enforce_boundary_condition_south_matrix(A)
    return A

def enforce_boundary_condition_south_matrix(A):
    A[0, 0:3] = [3, -4, 1] # Neumann condition, second order
    # A[0, 0:2] = [1, -1] # Neumann condition, first order
    A[-1,-1] = 1            # Dirichlet condition

def enforce_boundary_condition_south_rhs(f, T_s):
    f[0] = 0  # Neumann condition
    f[-1] = T_s # Dirichlet condition

def compute_north_matrix(N, dx, x, D, B_out):
    return np.flip(compute_south_matrix(N, dx, x, D, B_out))


def compute_south_rhs(N, x, A_out, Q, S, a, T_s):
    """
    compute the right-hand side of the numerical scheme.
    """
    rhs = np.empty(N+1)
    rhs[1:-1] = -A_out + Q * S(x[1:-1]) * a
    enforce_boundary_condition_south_rhs(rhs, T_s)
    return rhs

def compute_north_rhs(N, x, A_out, Q, S, a, T_s):
    return np.flip(compute_south_rhs(N, x, A_out, Q, S, a, T_s))


def S(x, S_2=-0.477):
    return 1 + 0.5*S_2 * (3*x**2 - 1)

def compute_temperature(N, x_s):
    N_south = int(N * x_s)
    N_north = (N - N_south)*10
    x_south, dx_south = np.linspace(0, x_s, N_south+1, retstep=True)
    x_north, dx_north = np.linspace(x_s, 1, N_north+1, retstep=True)
    D = 0.3
    a_south = 0.68
    a_north = 0.38
    A_out = 201.4
    B_out = 1.45
    T_s = 0
    Q = 340

    f_south = compute_south_rhs(N_south, x_south, A_out, Q, S, a_south, T_s)
    A_south = compute_south_matrix(N_south, dx_south, x_south, D, B_out)
    T_south = np.linalg.solve(A_south, f_south)

    f_north = compute_north_rhs(N_north, x_north, A_out, Q, S, a_north, T_s)
    A_north = compute_north_matrix(N_north, dx_north, x_north, D, B_out)
    T_north = np.linalg.solve(A_north, f_north)

    return T_south, x_south, T_north, x_north


def main():
    N = 200
    x_s = 0.9

    fig, ax = plt.subplots()
    T_south, x_south, T_north, x_north = compute_temperature(N, x_s)
    
    line_south, = ax.plot(np.asin(x_south), T_south, label='South')
    line_north, = ax.plot(np.asin(x_north), T_north, label='North')
    ice_line = ax.axvline(np.arcsin(x_s), color='black', linestyle='--', label='Ice edge')
    zero_line = ax.axhline(0, color='black', linestyle='--')
    
    ax.set_xlabel('$\\varphi/\\text{rad}$')
    ax.set_ylabel('$T/C^{\circ}$')
    ax.legend()

    # plt.show()

    frames = 75
    x_low = 0.9
    x_high = 0.999

    def update(frame):
        # x_s_frame = x_s + (1+np.cos(2*np.pi*frame/frames))**2 * (1-2*x_s)
        phi_s_frame = np.asin(x_low) + (np.asin(x_high) - np.asin(x_low)) * (1 - abs(frame/frames - 1)**2)
        x_s_frame = np.sin(phi_s_frame)

        T_south, x_south, T_north, x_north = compute_temperature(N, x_s_frame)
        
        line_south.set_data(np.asin(x_south), T_south)
        line_north.set_data(np.asin(x_north), T_north)
        ice_line.set_xdata(np.arcsin(x_north[[0, 0]]))

        return line_south, line_north, ice_line

    ani = anim.FuncAnimation(fig, update, frames=2*frames, blit=True, repeat=True, interval=20)
    # writer = anim.PillowWriter(fps=30,
    #                            bitrate=1800)
    # ani.save('animation_70deg.gif', writer=writer)
    plt.show()

if __name__ == '__main__':
    main()
