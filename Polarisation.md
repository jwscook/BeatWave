$$
D(\omega_{choose}, k_{find})=0
$$

Find $k$

We know polarisation is in $E_x$ and $E_y=E_z=0$.
$$
\vec{\vec{D}}(\omega, k)\cdot \vec E=\vec0
$$
Find the nullspace of $\vec{\vec{D}}(\omega, k)$ gives you $\vec E$ and hence relative amplitudes of the components of the electric field.

For the extraordinary mode
$$
n^2=\frac{\epsilon_\bot^2-g^2}{\epsilon_\bot}\\
E_x\propto i g\\
E_y\propto -\epsilon_\bot
$$
where $n=kc/\omega$
$$
\epsilon_\bot = 1- \sum_s \frac{\Pi_s^2}{\omega^2 - \Omega_s^2},\\
g = -\sum_s \frac{\Pi_s^2\Omega_s}{\omega (\omega^2 - \Omega_s^2)}
$$
Hence find $k_\alpha$ and $k_\beta$ associated with $\omega_\alpha$ and $\omega_\beta = \omega_\alpha + \Omega_i$ respectively. 

Choose $\omega_\alpha=\omega_L+ \Delta\omega$ where $\Delta \omega$ is a tuning parameter to find the right mode and to fix the wavenumber at something sensible. $\omega_L\approx \sqrt{\frac{\Omega_e^2}{4}+\Pi_e^2} - \frac{\Omega_e}{2}$.

Now the other parameters:

Faraday's
$$
-i\omega\vec B_1 = - \vec k \times \vec E_1
$$
Gauss's law
$$
i\vec k \cdot \vec E_1 =  \sum_s\frac{q_s n_{s,1}}{\epsilon_0}\\
$$
Ampere's law:
$$
\mu_0 \vec J = i \vec k \times \vec {B}_1 + i \omega \mu_0 \epsilon_0 \vec E_1\\
\sum_s q_s n_{s,0} \vec {v}_{s,1} = \frac{i \vec k \times \vec {B}_1}{\mu_0} + i \omega \epsilon_0 \vec E_1\\
$$
Assume the ions are unperturbed in the initial conditions
$$
q_e n_{e,0} {v}_{e,1,y} = i  \omega \epsilon_0 E_{1,x}\\
q_e n_{e,0} {v}_{e,1,y} = -\frac{i k {B}_{1,z}}{\mu_0} + i \omega \epsilon_0 E_{1,y}\\
q_e n_{e,0} {v}_{e,1,z} = \frac{i k {B}_{1,y}}{\mu_0} + i \omega \epsilon_0 E_{1,z}\\
$$

$$
n_{e,1} = i k E_{x,1}\frac{\epsilon_0}{q_s}
$$





