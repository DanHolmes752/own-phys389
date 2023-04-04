from Particle import Particle
from Simulation import Cyclotron
from cyclotronSimulation import EM_Cyclotron
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants
from testingSim import gyrationTests
from testingConservation import energyConservation

# Constants
dt = 1e-10
t_max = 5e-6
B = np.array([0, 0, 0.5])
numOfPro = 1
cyclotronRadius = 0.025 #0.0015
E = 0 #5e4
E_boundary = 0.05 * cyclotronRadius


### testingSim.py:

#gyrationTests.calculateData(dt = dt, t_max = t_max, B = B, numOfPro=numOfPro)
gyrationTests.calculate_EM_Data(dt = dt, t_max = t_max, B = B, numOfPro=numOfPro, cyclotronRadius = cyclotronRadius, E = E, E_boundary = E_boundary)
# gyrationTests.eulerPlot()
# gyrationTests.cromerPlot()
# gyrationTests.KuttaPlot()
gyrationTests.comparePlots()
#gyrationTests.sinTest(dt = dt, t_max = t_max, B = B)



### testingConservation.py:

# energyConservation.calcAngularFrequency(x_positions, t_list)
# energyConservation.eulerAngularFrequency(dt = dt, t_max = t_max)
# energyConservation.cromerAngularFrequency(dt = dt, t_max = t_max)
# energyConservation.kuttaAngularFrequency(dt = dt, t_max = t_max)
# energyConservation.compareAngular(dt = dt, t_max = t_max)

#     # energyConservation.eulerKineticTest(dt = dt, t_max = t_max)
#     # energyConservation.cromerKineticTest(dt = dt, t_max = t_max)
#     # energyConservation.kuttaKineticTest(dt = dt, t_max = t_max)
energyConservation.compareKineticTests(dt = dt,t_max = t_max)


