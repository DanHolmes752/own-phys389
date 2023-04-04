from Particle import Particle
from Proton import Proton
from Electron import Electron
import numpy as np


pos = np.array([1, 2, 3])
vel = np.array([0, 10, 0])
acc = np.array([0, 0, 0])


print("---------------------------------------------------------------------------------------------")

# Let's create an instance of a Proton object and let it print it's given attributes
pro_test = Proton(pos, vel, acc)
print(f"I am a {pro_test.name} at position {pro_test.position}, with velocity {pro_test.velocity} and acceleration {pro_test.acceleration}")
print(f"I also have mass {pro_test.mass} kg and charge {pro_test.charge} C")

print("---------------------------------------------------------------------------------------------")

# Let's create an instance of a Electron object and let it print it's given attributes
elec_test = Electron(pos*2, vel*3, acc)
print(f"I am an {elec_test.name} at position {elec_test.position}, with velocity {elec_test.velocity} and acceleration {elec_test.acceleration}")
print(f"I also have mass {elec_test.mass} kg and charge {elec_test.charge} C")

print("---------------------------------------------------------------------------------------------")