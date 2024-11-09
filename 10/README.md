# Question 10

We want to use finite differences to solve the equation

$$
    -D \frac{\partial}{\partial x} \left( (1-x^2) \frac{\partial T}{\partial x} \right) = -I(x) + Q S(x) a(x, x_s)
$$

in the intervals $[0, x_s]$ and $[x_s, 1]$, with the intensity being approximated by

$$I(x) = B_\text{out} T(x) + A_\text{out}.$$

## Finite difference operators

Define the central difference operator

$$
    \delta_x f_j = \frac{f_{j+1/2} - f_{j-1/2}}{\Delta x},
$$

the forward difference operator

$$
    \Delta_{x} f_j = \frac{f_{j+1} - f_j}{\Delta x},
$$

and the backward difference operator

$$
    \nabla_{x} f_j = \frac{f_j - f_{j-1}}{\Delta x}.
$$

## Discretization

We can now use the finite difference operators to discretize the equation in the interior:

$$-D \delta_x \left[ (1-x^2) \delta_x T \right]_j = -I_j + Q S_j a_j$$

This results in the equations

$$-D \left(\frac{1-x_{j+1/2}^2}{\Delta x} \frac{T_{j+1}-T_j}{\Delta x} - \frac{1-x_{j-1/2}^2}{\Delta x} \frac{T_j-T_{j-1}}{\Delta x}\right) = -I_j + Q S_j a_j,$$

for $j=1, 2, \ldots, N-1$. On the boundaries, we could use the forward and backward difference operators and implicitly embed the boundary conditions in the equations. However, it is easier to just explicitly enforce the boundary conditions.


## Boundary conditions

At the ice cap $x_s$, we have $T=T_s=0^\circ C$ and at the pole and equator, we use Neumann boundary conditions $\frac{\partial T}{\partial x} = 0$. This can be discretized as $T_0 = T_1$ in the first domain and $T_N = T_{N-1}$ in the second domain.

### Higher order boundary conditions

The central difference operator is of second dergree, but that helps little, when the boundary conditions are of first degree. A second order approximation of the Neumann boundary conditions would be

$$3T_0 - 4T_1 + T_2 = 0$$

and

$$3T_N - 4T_{N-1} + T_{N-2} = 0.$$


## Matrix Reformulation

We first consider the domain $[0,x_s]$. In the interior, the equations can be rewritten as

$$\alpha_j T_{j-1} + \beta_j T_j + \gamma_j T_{j+1} = f_j,$$

where

$$
\begin{aligned}
    \alpha_j &= - D \frac{1-x_{j-1/2}^2}{\Delta x^2}, \\
    \beta_j &= D \frac{2 - x_{j+1/2}^2 - x_{j-1/2}^2}{\Delta x^2} + B_\text{out} \\
    \gamma_j &= - D \frac{1-x_{j+1/2}^2}{\Delta x^2}, \\
    f_j &= -A_\text{out} + Q S_j a_j.
\end{aligned}
$$

Additionally, we apply the boundary conditions given above. We can now write this system of equations as

$$A \mathbf{T} = \mathbf{f},$$

where $\mathbf{T} = [T_0, T_1, \ldots, T_N]^T$, $\mathbf{f} = [f_0, f_1, \ldots, f_{N-1}, f_N]^T$, 

$$A = \begin{bmatrix}
    3 & -4 & 1 & 0 &\cdots & 0 & 0 & 0 & 0 \\
    \alpha_1 & \beta_1 & \gamma_1 & 0 & \cdots & 0 & 0 & 0 & 0 \\
    0 & \alpha_2 & \beta_2 & \gamma_2 & \cdots & 0 & 0 & 0 & 0 \\
    & & & &\ddots & & & \\
    0 & 0 & 0 & 0 & \cdots & 0 & \alpha_{N-1} & \beta_{N-1} & \gamma_{N-1} \\
    0 & 0 & 0 & 0 & \cdots & 0 & 0 & 0 & 1
\end{bmatrix},$$

and $f_0 = \frac{\partial T}{\partial x}=0$ and $f_N = T_s$ enforce the boundary conditions.


Above the ice cap, $[x_s, 1]$, due to symmetry, we just flip theindices of $\mathbf{f}$ and $A$:

$$
\begin{aligned}
    f_j &\gets f_{N+1-j} \\
    A_{j,k} &\gets A_{N+1-j, N+1-k}
\end{aligned}
$$