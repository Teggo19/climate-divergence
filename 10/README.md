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
    \delta_x f_j = \frac{f_{j+1/2} - f_{j-1/2}}{\Delta x}
$$

and the backward difference operator

$$
    \Delta_x f_j = \frac{f_j - f_{j-1}}{\Delta x}.
$$

## Boundary conditions

At the ice cap $x_s$, we have $T=T_s=0^\circ C$ and at equator, we use Neumann boundary conditions $\frac{\partial T}{\partial x} = 0$. This can be discretized as $T_0 = T_1$ in the first domain and $T_N = T_{N-1}$ in the second domain. At the pole, we do not have a boundary condition, so we need to an explicit integrator to determine the temperature there, say backward difference $\Delta_x$. This gives

$$\tilde\alpha_N T_{N-2} + \tilde\beta_N T_{N-1} + \tilde\gamma T_N = f_N,$$

where 

$$
\begin{aligned}
    \tilde\alpha_N &= D \frac{1-x_{N-1}^2}{\Delta x^2}, \\
    \tilde\beta_N &= - D \frac{2 - x_N^2 - x_{N-1}^2}{\Delta x^2} \\
    \tilde\gamma_N &= D \frac{1-x_N^2}{\Delta x^2} + B_{out}, \\
    f_N &= -A_\text{out} + Q S_N a_N.
\end{aligned}
$$

## Interior discretization

We can now use the finite difference operators to discretize the equation in the interior:

$$-D \delta_x \left[ (1-x^2) \delta_x T \right]_j = -I_j + Q S_j a_j$$

This results in the equations

$$-D \left(\frac{1-x_{j+1/2}^2}{\Delta x} \frac{T_{j+1}-T_j}{\Delta x} - \frac{1-x_{j-1/2}^2}{\Delta x} \frac{T_j-T_{j-1}}{\Delta x}\right) = -I_j + Q S_j a_j,$$

for $j=1, 2, \ldots, N-1$.

## Matrix Reformulation

The equations above can be rewritten as

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


In the domain $[x_s, 1]$ above the ice cap, we have almost the same system of equations

$$A = \begin{bmatrix}
    1 & 0 & 0 & 0 &\cdots & 0 & 0 & 0 & 0 \\
    \alpha_1 & \beta_1 & \gamma_1 & 0 & \cdots & 0 & 0 & 0 & 0 \\
    0 & \alpha_2 & \beta_2 & \gamma_2 & \cdots & 0 & 0 & 0 & 0 \\
    & & & &\ddots & & & \\
    0 & 0 & 0 & 0 & \cdots & 0 & \alpha_{N-1} & \beta_{N-1} & \gamma_{N-1} \\
    0 & 0 & 0 & 0 & \cdots & 0 & \tilde\alpha_N & \tilde\beta_N & \tilde\gamma_N
\end{bmatrix}$$

with $f_0 = T_s$ and $f_N = -A_\text{out} + Q S_N a_N$.