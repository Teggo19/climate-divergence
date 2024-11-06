# Question 10

We want to use finite differences to solve the equation

$$
    -D \frac{\partial}{\partial x} \left( (1-x^2) \frac{\partial T}{\partial x} \right) = -I(x) + Q S(x) a(x, x_s)
$$

in the intervals $[0, x_s]$ and $[x_s, 1]$.

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

We can now use the finite difference operators to discretize the equation:

$$-D \mathcal L_x \left[ (1-x^2) \mathcal L_x T \right]_j = -I_j + Q S_j a_j,$$

where $\mathcal L_x \in \{ \delta_x, \Delta_x, \nabla_x \}$. In the interior, we use $\mathcal L_x = \delta_x$:

$$-D \left(\frac{1-x_{j+1/2}^2}{\Delta x} \frac{T_{j+1}-T_j}{\Delta x} - \frac{1-x_{j-1/2}^2}{\Delta x} \frac{T_j-T_{j-1}}{\Delta x}\right) = -I_j + Q S_j a_j$$

On the left boundary (equator or ice cap), we use $\mathcal L_x = \Delta_x$:

$$-D\left(\frac{1-x_1^2}{\Delta x} \frac{T_2-T_1}{\Delta x} - \frac{1-x_0^2}{\Delta x} \frac{T_1-T_0}{\Delta x} \right) = -I_0 + Q S_0 a_0$$

On the right boundary (ice cap or pole), we use $\mathcal L_x = \nabla_x$:

$$-D\left(\frac{1-x_N^2}{\Delta x} \frac{T_N-T_{N-1}}{\Delta x} - \frac{1-x_{N-1}^2}{\Delta x} \frac{T_{N-1}-T_{N-2}}{\Delta x} \right) = -I_N + Q S_N a_N$$


## Boundary conditions

We consider first the region $[0, x_s]$. The following also applies symmetrically to the region $[x_s, 1]$. At the ice cap $x_s$, we have $T=T_s=273.15 K$ and at the pole and equator, we Neumann boundary conditions $\frac{\partial T}{\partial x} = 0$. 

These can be discretized as $T_0 = T_1$ and $T_N = T_{N-1}$. This simplifies the second equation to

$$-D\frac{1-x_1^2}{\Delta x} \frac{T_2-T_1}{\Delta x} = -I_0 + Q S_0 a_0$$


## Matrix Reformulation

In the interior, the equations can be rewritten as

$$\alpha_j T_{j-1} + \beta_j T_j + \gamma_j T_{j+1} = f_j,$$

where

$$
\begin{aligned}
    \alpha_j &= - D \frac{1-x_{j-1/2}^2}{\Delta x^2}, \\
    \beta_j &= D \frac{2 - x_{j+1/2}^2 - x_{j-1/2}^2}{\Delta x^2}, \\
    \gamma_j &= - D \frac{1-x_{j+1/2}^2}{\Delta x^2}, \\
    f_j &= -I_j + Q S_j a_j.
\end{aligned}
$$

At the equator, we have

$$\beta_0 T_1 + \gamma_0 T_2 = f_0,$$

with $f_0$ as above and

$$
\begin{aligned}
    \beta_1 &= D \frac{1-x_1^2}{\Delta x^2}, \\
    \gamma_1 &= - D \frac{1-x_1^2}{\Delta x^2}.
\end{aligned}
$$

At the ice cap, we have

$$\alpha_N T_{N-2} + \beta_N T_{N-1} = f_N,$$

with

$$
\begin{aligned}
    \alpha_N &= - D \frac{2-x_{N-1}^2}{\Delta x^2}, \\
    \beta_N &= D \frac{2-x_s^2-x_{N-1}^2}{\Delta x^2}, \\
    f_N &= D \frac{1-x_s^2}{\Delta x^2} T_s - I_N + Q S_N a_N.
\end{aligned}
$$

We can now write this system of equations as

$$A \mathbf{T} = \mathbf{f},$$

where $\mathbf{T} = [T_0, T_1, \ldots, T_N]^T$, $\mathbf{f} = [f_0, f_1, \ldots, f_N]^T$, and $A=\operatorname{tridiag}(\alpha, \beta, \gamma)$:

$$A = \begin{bmatrix}
    \beta_0 & \gamma_0 & 0 & 0 &\cdots & 0 & 0 & 0 \\
    \alpha_1 & \beta_1 & \gamma_1 & 0 & \cdots & 0 & 0 & 0 \\
    & & & &\ddots & & & \\
    0 & 0 & 0 & 0 & \cdots & \alpha_{N-1} & \beta_{N-1} & \gamma_{N-1} \\
    0 & 0 & 0 & 0 & \cdots & 0 & \alpha_N & \beta_N
\end{bmatrix}$$