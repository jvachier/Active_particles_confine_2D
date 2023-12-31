# Active particles confined 2D
The aim of this project is to build simulations describings the motion of active interacting particles under confinement in 2D. These simulations are based on Langevin equations and used the Euler-Mayurama algorithm. Two types of geometries are used: either squared or circular.

Active interactive particles evolve in different geometries, such as circular or squared. The dynamics is given by two Langevins equations, one for the position $\mathbf{\tilde{r}}(\tilde{x},\tilde{y})$ of the particles and one for its orientation $\phi$

$$
\begin{align}
\frac{d}{d\tilde{t}}\mathbf{\tilde{r}} &= \tilde{v_s}\mathbf{\hat{e}} - \tilde{\nabla}_{\tilde{R}}(\tilde{LP}) + \sqrt{2\tilde{D}_t}\tilde{\mathbf{\xi}_t}\,\\
\frac{d}{d\tilde{t}}\phi &= \sqrt{2\tilde{D}_e}\tilde{\mathbf{\xi}_e}\,
\end{align}
$$

where $\mathbf{\hat{e}} = (\cos(\phi),\sin(\phi))^{T}$ is the orientational unit vector, $\tilde{v_s}$ is the self-propulsion, $\tilde{D_{t}}$ and $\tilde{D_{e}}$ are the translational and rotational diffusivities, respectively. Moreover, $\langle \tilde{\xi_{t_i}}(\tilde{t}')\tilde{\xi_{t_j}}(\tilde{t}) \rangle = \delta_{ij}\delta(\tilde{t}'-\tilde{t})$ and $\langle \tilde{\xi_{e_i}}(\tilde{t}')\tilde{\xi_{e_j}}(\tilde{t}) \rangle = \delta_{ij}\delta(\tilde{t}'-\tilde{t})$ are two Gaussian white noises. Moreover, the interactions between the particles is represented by using the Lennard-Jones  potential

$$
\tilde{LP} = 4\tilde{\epsilon}[(\frac{\tilde{\sigma}}{\tilde{R}})^{12} - (\frac{\tilde{\sigma}}{\tilde{R}})^{6}]\,
$$

where $\tilde{\epsilon}$ is the depth of the potential well, $\tilde{R}$ is the distance between two interacting particles. In this project, only the repulsive part of the potential is considered.


https://github.com/jvachier/Active_particles_confine_2D/assets/89128100/05b6ead5-0880-4eb5-8825-da8b687c240c

