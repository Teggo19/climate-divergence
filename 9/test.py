import numpy as np
import matplotlib.pyplot as plt

a_l = 0.68
a_u = 0.38
S_2 = -0.477
D = 1e-2
A = 201.4
B = 1.45

a_m = (a_l + a_u)/2

T_s = 0

x_s = 0.9

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

#print(np.where(np.array(Q_arr) > Q_real))
print(Q_real)
print(x[np.where(np.array(Q_arr) > Q_real)[0][0]])
print(Q(x[np.where(np.array(Q_arr) > Q_real)[0][0]-1]))
#print(Q_arr[37])
#print(Q_arr[38])
plt.figure(figsize=(6, 8))
plt.title(r"Intersection between $Q(x_s)$ and $G_{sc}/4$")
plt.plot(x, Q_arr, label=r"$Q(x_s)$")
# plt.plot(x, [Q_new(xi) for xi in x], label="Q(x_s)")
plt.axhline(y=Q_real, color='r', linestyle='--', label=r"340 W/m$^2$")
plt.axvline(x=x[np.where(np.array(Q_arr) > Q_real)[0][0]], color='g', linestyle='--', label=rf"$x_s = 0.44$")
plt.legend()
plt.xlabel(r"$x_s$")
plt.ylabel(r"Q (W/m$^2$)")
plt.savefig("Q_plot.png", dpi=300)
plt.show()

