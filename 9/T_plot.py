import numpy as np
import matplotlib.pyplot as plt

a_l = 0.68
a_u = 0.38
S_2 = -0.477
D = 2
A = 201.4
B = 1.45

Q = 340

a_l = (a_l + a_u)/2
a_u = a_l

T_s = -20

x_s = 0.95

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
    elif x == x_s:
        return (a_l + a_u)/2
    else:
        return a_u
    
def T_p(x):
    return (Q*a(x) - A)/B + (Q*a(x)*S_2)/(6*D + B)*0.5*(3*x**2 - 1)

C = (T_s- T_p(x_s))/y(x_s)
E = (T_s - T_p(x_s))/z(x_s)

def T(x):
    if x < x_s:
        return T_p(x) + C*y(x)
    else:
        return T_p(x) + E*z(x)

print(x_s)

x = np.linspace(0, 1, 100)
T_arr = [T(xi) for xi in x]
y_arr = [y(xi) for xi in x]
T_p_arr = [T_p(xi) for xi in x]
plt.plot(x, T_arr, label="T(x)")
#plt.plot(x, y_arr, label="y(x)")
#plt.plot(x, T_p_arr, label="T_p(x)")
plt.legend()
plt.show()