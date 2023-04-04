# PHYS389 Project Outline - Cyclotron Simulation - Daniel Holmes

## Physics:

A cyclotron is a type of particle accelerator which works by using electromagnetic fields to accelerate charged particles along a spiral path, projecting them outwards from the centre of a cylindrical vacuum.

Looking to:
- Simulate charged particles in electric and magnetic fields.
- Verify energy conservation.

The acceleration, and subsequently velocity and position, can be calculated from the Lorentz force:

$$
\vec{F} = q(\vec{E}+\vec{v}\times\vec{B})
$$

Using that force is the rate of change of momentum:

$$
\implies \frac{d \vec{p}}{dt}=\frac{md\vec{v}}{dt} = q(\vec{E}+\vec{v}\times \vec{B})
$$

$$
\implies \frac{d\vec{v}}{dt}=\frac{q}{m}(\vec{E}+\vec{v} \times \vec{B})
$$

To have the electric field continually accelerate the charged particles, it will have to oscilate with the same frequency as the particles orbit. This can be achieved with a sinusoidal time dependence (t) at a fixed frequency ($\omega$), for example:

$$
\vec{E} = E_0 sin(\omega t+\phi)\hat{j}
$$

where $\phi$ is the phase shift and:

$$
\omega = \frac{qB}{m}
$$

$$
\implies T = \frac{2\pi m}{qB}
$$

Thus for the same type of particle, $\omega$ and T should be constant due to the fixed magnetic field (B).

Boundary conditions:
- Electric field between -L < y < L parallel to the y-axis.
- Magnetic field in the +z direction should be contained within the radii of the two D shaped cavities.



## Algorithim:

To solve the Lorentz equation numerically, we can use the Euler Forward Method:

$$
\frac{d\vec{v}}{dt}=\frac{q}{m}(\vec{E}+\vec{v}\times\vec{B})
$$

$$
\implies \frac{\Delta \vec{v}} {\Delta t}=\frac{q}{m} (\vec{E}+\vec{v}\times\vec{B})
$$

$$
\implies \vec{v_{n+1}}-\vec{v_{n}} = (\frac{q}{m}(\vec{E} + \vec{v}\times\vec{B})){\Delta t}
$$

$$
\implies \vec{v_{n+1}}=(\frac{q}{m}(\vec{E}+\vec{v} \times \vec{B})){\Delta t}+\vec{v_{n}}
$$

$$
\implies \vec{v_{n+1}} = (\vec{a}){\Delta t}+\vec{v_{n}}
$$

This will calculate the updated velocity from the acceleration caused by the Lorentz force. From the change in velocity the updated position can also be found.

Euler Forward is a simple integration method, which may result in inaccurate calculations over a large number of gyrations. For this reason other methods will also be used, such as Runge-Kutta and Verlet, to compare against each other.

## Testing:
- Simulate a charged particle in a constant magnetic field, does it achieve the expected orbit and period?
- Expand to include multiple particles and then an oscillating electric field.
- Analysing numerical data: first comparing different integration methods then comparing to analytical values - find % difference.
- Position of particles can be graphed over the time step to have a visual understanding of the path taken.
- Energy values can be calculated and graphed to test conservation laws.


## Initial Code & Class Ideas:
- Single particle class that updates particle properties, ie acceleration, velocity and position.
- Different methods or classes for different numerical methods: Euler Forward, Runge-Kutta, Verlet.
- Class or Python script that can initalise a randomised bunch of charged particles, Monte Carlo method could be used here.
- Class or Python script that saves numerical data to a file for each numerical method.
- Python script that can run the simulation from a couple simple lines, allowing easy changing of parameters.
- Use imports use as numpy and matplotlib to help with numerical calculations and plotting.
