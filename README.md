# Active particles confined 2D
The aim of this project is to build simulations to describe the motion of active interacting particles under confinement in 2D. These simulations are based on Langevin equations and used the Euler-Mayurama algorithm. Two types of geometries will be used: either squared or circular.

Active interactive particles evolve in different geometries, such as circular or squared. The dynamics is given by two Langevins equations, one for the position $r(x,y)$ of the particles and for its orientation $\phi$

$$
\begin{align}
\frac{d}{dt}\bm{\tilde{r}} &= \tilde{v_s}\bm{\tilde{e}} + \sqrt{2\tilde{D}_t}\bm{\tilde{\xi_t}}\,,\\
\frac{d}{dt}\phi &= \sqrt{2\tilde{D}_e}\bm{\tilde{\xi_e}}\,,
\end{align}
$$
where $\bm{\tilde{e}} = (\cos(\phi),\sin(\phi))^{T}$, $\langle \xi_{t_i}(t')\xi_{t_j}(t) \rangle = \delta_{ij}\delta(t'-t)$ and $\langle \xi_{e_i}(t')\xi_{e_j}(t) \rangle = \delta_{ij}\delta(t'-t)$ two Gaussian white noises.