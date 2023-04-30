# Active particles confined 2D
The aim of this project is to build simulations describings the motion of active interacting particles under confinement in 2D. These simulations are based on Langevin equations and used the Euler-Mayurama algorithm. Two types of geometries will be used: either squared or circular.

Active interactive particles evolve in different geometries, such as circular or squared. The dynamics is given by two Langevins equations, one for the position $\mathbf{\tilde{r}}(\tilde{x},\tilde{y})$ of the particles and for its orientation $\phi$

$$
\begin{align}
\frac{d}{d\tilde{t}}\mathbf{\tilde{r}} &= \tilde{v_s}\mathbf{\hat{e}} + \sqrt{2\tilde{D}_t}\tilde{\mathbf{\xi}_t}\,\\
\frac{d}{d\tilde{t}}\phi &= \sqrt{2\tilde{D}_e}\tilde{\mathbf{\xi}_e}\,
\end{align}
$$

where $\mathbf{\hat{e}} = (\cos(\phi),\sin(\phi))^{T}$, $\tilde{v_s}$ the self-propulsion, $\tilde{D_{t}}$ and $\tilde{D_{e}}$ the translational and rotational diffusivities, respectively. Moreover, $\langle \tilde{\xi_{t_i}}(\tilde{t}')\tilde{\xi_{t_j}}(\tilde{t}) \rangle = \delta_{ij}\delta(\tilde{t}'-\tilde{t})$ and $\langle \xi_{e_i}(\tilde{t}')\xi_{e_j}(\tilde{t}) \rangle = \delta_{ij}\delta(\tilde{t}'-\tilde{t})$ two Gaussian white noises.
