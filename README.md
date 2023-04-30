# Active particles confined 2D
The aim of this project is to build simulations describings the motion of active interacting particles under confinement in 2D. These simulations are based on Langevin equations and used the Euler-Mayurama algorithm. Two types of geometries will be used: either squared or circular.

Active interactive particles evolve in different geometries, such as circular or squared. The dynamics is given by two Langevins equations, one for the position $r(x,y)$ of the particles and for its orientation $\phi$

$$
\begin{align}
\frac{d}{dt}\mathbf{\tilde{r}} &= \tilde{v_s}\mathbf{\hat{e}} + \sqrt{2\tilde{D}_t}\tilde{\mathbf{\xi}_t}\,\\
\frac{d}{dt}\phi &= \sqrt{2\tilde{D}_e}\tilde{\mathbf{\xi}_e}\,
\end{align}
$$

where $\mathbf{\hat{e}} = (\cos(\phi),\sin(\phi))^{T}$, $\langle \xi_{t_i}(t')\xi_{t_j}(t) \rangle = \delta_{ij}\delta(t'-t)$ and $\langle \xi_{e_i}(t')\xi_{e_j}(t) \rangle = \delta_{ij}\delta(t'-t)$ two Gaussian white noises.
