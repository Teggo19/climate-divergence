import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim

from trygve import T as T_analytic

def compute_diagonal(N, dx, x, D, B_out):
    """
    compute the diagonal of the matrix for the numerical scheme.
    """
    return D * (2 - (x-0.5*dx)**2 - (x+0.5*dx)**2) / dx**2 + B_out


def compute_lower(N, dx, x, D):
    """
    compute the lower diagonal of the matrix for the numerical scheme.
    """
    return -D * (1 - (x[1:]-0.5*dx)**2) / dx**2


def compute_upper(N, dx, x, D):
    """
    compute the upper diagonal of the matrix for the numerical scheme.
    """
    return -D * (1 - (x[:-1]+0.5*dx)**2) / dx**2


def compute_central_difference_matrix(N, dx, x, D, B_out):
    lower = compute_lower(N, dx, x, D)
    upper = compute_upper(N, dx, x, D)
    diag = compute_diagonal(N, dx, x, D, B_out)
    return np.diag(diag) + np.diag(lower, -1) + np.diag(upper, 1)


def compute_south_matrix(N, dx, x, D, B_out):
    """
    compute a matrix for the numerical scheme.
    """
    A = compute_central_difference_matrix(N, dx, x, D, B_out)
    enforce_boundary_condition_south_matrix(A)
    return A

def enforce_boundary_condition_south_matrix(A):
    A[0, 0:3] = [3, -4, 1] # Neumann condition, second order
    # A[0, 0:2] = [1, -1] # Neumann condition, first order
    A[-1, -3:] = [0, 0, 1] # Dirichlet condition

def enforce_boundary_condition_south_rhs(f, T_s):
    f[0] = 0  # Neumann condition
    f[-1] = T_s # Dirichlet condition

def enforce_boundary_condition_north_rhs(f, T_s, dT_s):
    f[0] = T_s  #Dirichlet condition
    f[1] = dT_s  # Neumann condition

def compute_north_matrix(N, dx, x, D, B_out):
    A = compute_central_difference_matrix(N, dx, x, D, B_out)
    enforce_boundary_condition_north_matrix(A, dx, x, D, B_out)
    return A


def enforce_boundary_condition_north_matrix(A, dx, x, D, B_out):
    A[0,0:2] = [1, 0]  # Dirichlet condition
    A[1, :3] = [-1, 1, 0]  # Neutral condition

    # Backward difference:
    A[-1,:] = 0
    A[-1, -2] = -2*D*(1-(x[-1]-0.5*dx)**2)/dx**2
    A[-1, -1] = 2*D*(1-(x[-1]-0.5*dx)**2)/dx**2 + B_out


def compute_rhs(N, x, A_out, Q, S, a, T_s):
    rhs = np.empty(N+1)
    rhs[:] = -A_out + Q * S(x) * a
    return rhs


def compute_south_rhs(N, x, A_out, Q, S, a, T_s):
    """
    compute the right-hand side of the numerical scheme.
    """
    rhs = compute_rhs(N, x, A_out, Q, S, a, T_s)
    enforce_boundary_condition_south_rhs(rhs, T_s)
    return rhs

def compute_north_rhs(N, x, A_out, Q, S, a, T_s, dT_s):
    rhs = compute_rhs(N, x, A_out, Q, S, a, T_s)
    enforce_boundary_condition_north_rhs(rhs, T_s, dT_s)
    return rhs


def S(x, S_2=-0.477):
    return 1 + 0.5*S_2 * (3*x**2 - 1)

def compute_temperature(N, x_s):
    N_south = int(N * x_s)
    N_north = N - N_south
    x_south, dx_south = np.linspace(0, x_s, N_south+1, retstep=True)
    x_north, dx_north = np.linspace(x_s, 1, N_north+1, retstep=True)
    D = 1e-2
    a_south = 0.68
    a_north = 0.38
    A_out = 201.4
    B_out = 1.45
    T_s = 0
    Q = 340

    f_south = compute_south_rhs(N_south, x_south, A_out, Q, S, a_south, T_s)
    A_south = compute_south_matrix(N_south, dx_south, x_south, D, B_out)
    T_south = np.linalg.solve(A_south, f_south)

    f_north = compute_north_rhs(N_north, x_north, A_out, Q, S, a_north, T_s, (T_south[-1] - T_south[-2])/dx_south*dx_north)
    A_north = compute_north_matrix(N_north, dx_north, x_north, D, B_out)
    T_north = np.linalg.solve(A_north, f_north)

    return T_south, x_south, T_north, x_north


a_l = 0.68
a_u = 0.38
S_2 = -0.477
D = 1e-2
A = 201.4
B = 1.45
T_s = 0

N = 1000
beta_0 = 1
beta_1 = B/(2*D)
beta = [beta_0, beta_1]
for i in range(N):
    m = len(beta)
    beta.append(beta[-1]*(B+m*(m-1)*D)/(2*(m)**2*D))

alpha_0 = 1
alpha_1 = 0
alpha = [alpha_0, alpha_1]
for i in range(N):
    m = len(alpha)
    alpha.append(alpha[-2]*(B+D*(m-1)*(m-2))/(D*(m-1)*m))

def z(u):
    return sum([b*(1-u)**i for i, b in enumerate(beta)])

def y(x):
    return sum([a*(x)**i for i, a in enumerate(alpha)])

def z_prime(u):
    return sum([(-1)*(i+1)*b*(1-u)**(i) for i, b in enumerate(beta[1:])])

def y_prime(x):
    return sum([(i+1)*a*(x)**(i) for i, a in enumerate(alpha[1:])])


def Q(x_s):
    numerator = (A/B + T_s)*(z_prime(x_s)/z(x_s) - y_prime(x_s)/y(x_s))
    den_1 = 3*S_2*x_s/(6*D+B)*(a_l - a_u)
    den_2 = - y_prime(x_s)/y(x_s)*(a_l/B + a_l*S_2/(6*D+B)*0.5*(3*x_s**2 -1))
    den_3 = z_prime(x_s)/z(x_s)*(a_u/B + a_u*S_2/(6*D+B)*0.5*(3*x_s**2 -1))
    return numerator/(den_1 + den_2 + den_3)

def Q_new(x_s):
    return (T_s + A/B)/(a_m/B + a_m*S_2/(6*D+B)*0.5*(3*x_s**2 -1))

Q_real = 1360/4

x = np.linspace(0.1, 0.8, 1000)
Q_arr = [Q(xi) for xi in x]


fig, (ax0, ax,ax2) = plt.subplots(1, 3, figsize=(12, 4))

print(Q_real)
print(x[np.where(np.array(Q_arr) > Q_real)[0][0]])
print(Q(x[np.where(np.array(Q_arr) > Q_real)[0][0]-1]))
ax0.set_title(r"Intersection between $Q(x_s)$ and $G_{sc}/4$")
ax0.plot(x, Q_arr, label=r"$Q(x_s)$")
ax0.axhline(y=Q_real, color='r', linestyle='--', label=r"340 W/m$^2$")
ax0.axvline(x=x[np.where(np.array(Q_arr) > Q_real)[0][0]], color='g', linestyle='--', label=rf"$x_s = 0.44$")
ax0.legend()
ax0.set_xlabel(r"$x_s$")
ax0.set_ylabel(r"Q [W/m$^2$]")


N = 1000
x_s = 0.44

T_south, x_south, T_north, x_north = compute_temperature(N, x_s)

ax.plot(x_south, T_south, label='South')
ax.plot(x_north, T_north, label='North')
ax.axvline(x=x[np.where(np.array(Q_arr) > Q_real)[0][0]], color='g', linestyle='--', label=rf"$x_s = 0.44$")
ax.axhline(0, color='black', linestyle='--')

ax2.plot(x_south, T_south - np.array([T_analytic(x) for x in x_south]))
ax2.plot(x_north, T_north - np.array([T_analytic(x) for x in x_north]))

ax.set_title('Temperature profile')
ax2.set_title('Difference from analytical solution')
ax.set_xlabel('$x$')
ax2.set_xlabel('$x$')
ax.set_ylabel('$T$ [C$^{\circ}$]')
ax.legend()

plt.tight_layout()

plt.savefig("done.png")
plt.show()