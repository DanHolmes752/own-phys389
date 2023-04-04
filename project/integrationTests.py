from Particle import Particle
from Proton import Proton
from Electron import Electron
import numpy as np
import matplotlib.pyplot as plt


def testAcc():
    proton1 = Proton( np.array([0, 1, 0]), np.array([0, 1000, 0]), np.array([0, 0, 0]) )

    print(proton1.position)

    for x in range(0,10):
        proton1.updateAcceleration(B=[0,0,1], E = [0, 50000, 0])
        proton1.updateEulerCromer
        print(proton1.acceleration)
        print(proton1.velocity)


def testEulerFor():
    proton1 = Proton( np.array([0, 1, 0]), np.array([0, 10, 0]), np.array([0, 0, 0]) )

    print(proton1.position)

    for x in range(0,10):
        proton1.updateEulerForward(dt = 1)
        print(proton1.position)

def testEulerCro():
    proton1 = Proton( np.array([0, 1, 0]), np.array([0, 10, 0]), np.array([0, 0, 0]) )

    print(proton1.position)

    for x in range(0,10):
        proton1.updateEulerCromer(dt = 1)
        print(proton1.position)

def testRungeKut():

    proton1 = Proton( np.array([0, 1, 0]), np.array([0, 10, 0]), np.array([0, 0, 0]) )

    print(proton1.position)

    for x in range(0,10):
        proton1.updateRungeKutta(dt = 1)
        print(proton1.position)

testAcc()
# testEulerFor()
# testEulerCro()