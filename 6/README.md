TODO: tikz figur som forklarer koordinatene på jorda. Inkluderer også kontrollvolumet med volum $r\Delta\phi\cdot r\cos(\phi)\Delta\theta\cdot \Delta r$.

Using the control volume definition of divergence and the fact that $\nabla T = \frac{1}{r} \frac{\partial T}{\partial \phi}$, we obtain:

$$
\begin{align*}
    \nabla\cdot(k\nabla T)
    &= \lim_{V\rightarrow 0} \frac{\oiint_{\partial V} k\nabla T\cdot \boldsymbol{dS}}{\iiint_V dV}
    \\&=\lim_{V\rightarrow 0} \frac{
        \frac{1}{r} k(\phi+\Delta \phi) \frac{\partial T}{\partial \phi}(\phi+\Delta \phi)r\cos(\phi+\Delta\phi)\Delta\theta\Delta r
        -
        \frac{1}{r} k(\phi) \frac{\partial T}
        {\partial \phi}(\phi)r\cos(\phi)\Delta\theta\Delta r}{r^2\cos(\phi)\Delta\phi\Delta\theta\Delta r}
    \\&=\frac{1}{r^2\cos(\phi)}
    \lim_{\Delta\phi\rightarrow 0} \frac{
        k(\phi+\Delta \phi) \frac{\partial T}{\partial \phi}(\phi+\Delta \phi)\cos(\phi+\Delta\phi)
        -
        k(\phi) \frac{\partial T}
        {\partial \phi}(\phi)\cos(\phi) }{\Delta\phi}
    \\&=\frac{1}{r^2\cos(\phi)}
        \frac{\partial}{\partial\phi}\left( k\cos{\phi}\frac{\partial T}{\partial\phi} \right).
\end{align*}
$$

Noting that $x = \sin{\phi}$ the chain rule gives $\frac{\partial}{\partial\phi} = \frac{dx}{d\phi}\frac{\partial}{\partial x}$, where $\frac{dx}{d\phi} = \cos{\phi} = \sqrt{1-x^2}$ we can rewrite this expression as:

$$
    \nabla\cdot(k\nabla T) = \frac{1}{r^2}\frac{\partial}{\partial x} \left( k(1-x^2)\frac{\partial T}{\partial x} \right).
$$
