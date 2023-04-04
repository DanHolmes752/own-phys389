import numpy as np
import matplotlib.pyplot as plt

"""
# ## Testing:
# - Simulate a charged particle in a constant magnetic field, does it achieve the expected orbit and period?
# - Expand to include multiple particles and then an oscillating electric field.
# - Analysing numerical data: first comparing different integration methods then comparing to analytical values - find % difference.
# - Position of particles can be graphed over the time step to have a visual understanding of the path taken.
# - Energy values can be calculated and graphed to test conservation laws.


# ## Initial Code & Class Ideas:
# - Single particle class that updates particle properties, ie acceleration, velocity and position.
# - Different methods or classes for different numerical methods: Euler Forward, Runge-Kutta, Verlet.
# - Class or Python script that can initalise a randomised bunch of charged particles, Monte Carlo method could be used here.
# - Class or Python script that saves numerical data to a file for each numerical method.
# - Python script that can run the simulation from a couple simple lines, allowing easy changing of parameters.
# - Use imports use as numpy and matplotlib to help with numerical calculations and plotting.



Let's prototype a simulation of a proton traveling in the +y-direction of a mag field in the +z-diection
"""



# Constants
charge = 1.6e-19                                        # Charge of proton
mass = 1.67e-27                                         # Mass of proton
MagFieldStrength = np.array([0, 0, 0.5])                # Magnetic field strength vector
initial_velocity = 1e4                                  # Initial velocity

dt = 1e-10                                              # Time step
t_max = 1e-6
t_list = np.arange(0, t_max, dt, dtype=float)           # Create time step list
print("time step list:", t_list)


# Initial position and velocity
position_list = np.array([[0, 0, 0]])                   # Starting position at origin as a vector
velocity_list = np.array([[0, initial_velocity, 0]])    # Starting velocity in y-direction
print("inital position:", position_list)
print("initial velcoity:", velocity_list)



# Loop over number of time steps to update position & velocity at each step
for x in range(0, len(t_list)):

    velocity = np.array([velocity_list[x]])
    position = np.array([position_list[x]])

    # Calculate the magnetic force and acceleration
    Force = charge * np.cross(velocity, MagFieldStrength)
    acceleration = Force / mass

    # Update velocity & position
    updated_velocity = velocity + acceleration*dt
    position = position + updated_velocity*dt

    # Append this to the velocity & position lists
    velocity_list = np.append(velocity_list, updated_velocity, axis=0)
    position_list = np.append(position_list, position, axis=0)
    print(x, "velocity list:", velocity_list)
    print(x, "position list:", position_list)

print(position_list[:,0])
#plt.scatter(position_list[:,0], position_list[:,1], color = "red")
plt.plot(position_list[:,0], position_list[:,1], color = "red")
plt.show()
