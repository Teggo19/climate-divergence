import numpy as np
import matplotlib.pyplot as plt

a_l = 0.68
a_u = 0.38
S_2 = -0.477
D = 1e-2
A = 201.4
B = 1.45

Q = 340


T_s = 0

x_s = 0.41

beta_0 = 1
beta_1 = B/(2*D)
beta = [beta_0, beta_1]

N = 1000

for i in range(N):
    m = len(beta)
    beta.append(beta[-1]*(B+m*(m-1)*D)/(2*(m)**2*D))

alpha_0 = 1
alpha_1 = 0
alpha = [alpha_0, alpha_1]

for i in range(N):
    m = len(alpha)
    alpha.append(alpha[-2]*(B+D*(m-1)*(m-2))/(D*(m-1)*m))
print(alpha[:20])

def z(u):
    return sum([b*(1-u)**i for i, b in enumerate(beta)])

def y(x):
    return sum([a*(x)**i for i, a in enumerate(alpha)])

def z_prime(u):
    return sum([(-1)*(i+1)*b*(1-u)**(i) for i, b in enumerate(beta[1:])])

def y_prime(x):
    return sum([(i+1)*a*(x)**(i) for i, a in enumerate(alpha[1:])])

def a(x):
    if x < x_s:
        return a_l
    else:
        return a_u
    
def T_p(x):
    return (Q*a(x) - A)/B + (Q*a(x)*S_2)/(6*D + B)*0.5*(3*x**2 - 1)

delta = 1e-5

C = (T_s - T_p(x_s - delta))/y(x_s)
E = (T_s - T_p(x_s + delta))/z(x_s)

def T(x):
    if x <= x_s:
        return T_p(x) + C*y(x)
    else:
        return T_p(x) + E*z(x)


x = np.linspace(0, 1, 300)
T_arr = [T(xi) for xi in x]
y_arr = [y(xi) for xi in x]
T_p_arr = [T_p(xi) for xi in x]
plt.plot(x, T_arr, label="T(x)")
#plt.plot(x, y_arr, label="y(x)")
#plt.plot(x, T_p_arr, label="T_p(x)")
plt.axhline(y=T_s, color='r', linestyle='--', label=rf"T_s = {T_s}$\degree$C")
plt.axvline(x=x_s, color='g', linestyle='--', label=f"x_s = {x_s}")

plt.xlabel("x")
plt.ylabel(r"T ($\degree$C)")
plt.legend()

# plt.savefig("T_plot_analytic.png", dpi=300)
plt.show()