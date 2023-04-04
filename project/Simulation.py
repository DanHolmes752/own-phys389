from Particle import Particle
from Proton import Proton
from Electron import Electron
from particleBunch import particleBunch
import numpy as np
import matplotlib.pyplot as plt
import random
import copy


class Cyclotron:

    """
    A class to represent a constant magnetic field.

    Attributes:
    -----------
        B: 1D array, float entries.
           The magnetic field strength vector.
        dt: float
            The time step.
        t_max: float
            The total time that the simulation will be run over.

    Methods:
    --------
        runEulerSim(self, numOfPro):
            Runs the simulation using the Euler-Forward approximation, saving positional and KE data to a .npy file.
        runCromerSim(self, numOfPro):
            Runs the simulation using the Euler-Cromer approximation, saving positional and KE data to a .npy file.
        runKuttaSim(self, numOfPro):
            Runs the simulation using the Runge-Kutta approximation, saving positional and KE data to a .npy file.
    """

    def __init__(self, dt, t_max, B):

        """
        Constructs all the necessary attributes for the Cyclotron object.

        Parameters:
        -----------
            B : 1D array, float entries.
                The magnetic field strength vector.
            dt : float
                 The time step.
            t_max : float
                    The total time that the simulation will be run over.
        """

        self.B = B                                                      
        self.dt = dt                                                    
        self.t_max = t_max                                              
        self.t_list = np.arange(0, self.t_max, self.dt, dtype=float)    # List from 0 - maximum time in steps of dt.
    
    

    def runEulerSim(self, numOfPro):
        """
        Runs the simulation using the Euler-Forward approximation, saving positional and KE data to a .npy file.

        Parameters:
        -----------
            numOfPro (int) : Number of particles in the simulation.
        
        Returns:
        --------
            None

        """

        bunch = particleBunch.generateProtons(numOfPro, [0.0001, 0.0001, 0])                
        
        position_data = np.zeros((numOfPro, len(self.t_list), 3))                        # Size = layers,columns,rows = numProtons, num iterations, 3 coordinates xyz
        KE_data = np.zeros((numOfPro, len(self.t_list)))                                 # Size = columns,rows = numOfPro, no. of iterations

        # Starts simulation 
        for proton in bunch:

            position_list = np.array([proton.position])                                  # Array for proton positions
            KE_list = np.array([proton.kineticEnergy()])                                 # Array for proton KE's


            for j in range(0, len(self.t_list)-1):

                proton.updateAcceleration(self.B)
                proton.updateEulerForward(self.dt)

                KE = proton.kineticEnergy()

                position_list = np.append(position_list, [proton.position], axis=0)
                KE_list = np.append(KE_list, KE)


            position_data = np.append(position_data, [position_list], axis = 0)         # Appends position_list to a new layer in position_data
            KE_data = np.append(KE_data, [KE_list], axis = 0)                           # Appends KE_list to a new column in KE_data
            

        np.save("EF_particle_positions.npy", position_data)
        np.save("EF_particle_KE_data.npy", KE_data)


    def runCromerSim(self, numOfPro):
        """
        Runs the simulation using the Euler-Cromer approximation, saving positional and KE data to a .npy file.

        Parameters:
        -----------
            numOfPro (int) : Number of particles in the simulation.
        
        Returns:
        --------
            None
            
        """

        bunch = particleBunch.generateProtons(numOfPro, [0.0001, 0.0001, 0])                
        
        position_data = np.zeros((numOfPro, len(self.t_list), 3))                            # Size = layers,columns,rows = numProtons, num iterations, 3 coordinates xyz
        KE_data = np.zeros((numOfPro, len(self.t_list)))                                     # Size = columns,rows = numOfPro, no. of iterations

        # Starts simulation
        for proton in bunch:

            position_list = np.array([proton.position])                                  # Array for proton positions
            KE_list = np.array([proton.kineticEnergy()])                                 # Array for proton KE's


            for j in range(0, len(self.t_list)-1):

                proton.updateAcceleration(self.B)
                proton.updateEulerCromer(self.dt)

                KE = proton.kineticEnergy()

                position_list = np.append(position_list, [proton.position], axis=0)
                KE_list = np.append(KE_list, KE)


            position_data = np.append(position_data, [position_list], axis = 0)         # Appends position_list to a new layer in position_data
            KE_data = np.append(KE_data, [KE_list], axis = 0)                           # Appends KE_list to a new column in KE_data
            

        np.save("EC_particle_positions.npy", position_data)
        np.save("EC_particle_KE_data.npy", KE_data)


    
    def runKuttaSim(self, numOfPro):
        """
        Runs the simulation using the Runge-Kutta approximation, saving positional and KE data to a .npy file.

        Parameters:
        -----------
            numOfPro (int) : Number of particles in the simulation.
        
        Returns:
        --------
            None
            
        """

        bunch = particleBunch.generateProtons(numOfPro, [0.0001, 0.0001, 0])                
        
        position_data = np.zeros((numOfPro, len(self.t_list), 3))                       # Size = layers,columns,rows = numProtons, num iterations, 3 coordinates xyz
        KE_data = np.zeros((numOfPro, len(self.t_list)))                                # Size = columns,rows = numOfPro, no. of iterations

        # Starts simulation
        for proton in bunch:

            position_list = np.array([proton.position])                                  # Array for proton positions
            KE_list = np.array([proton.kineticEnergy()])                                 # Array for proton KE's


            for j in range(0, len(self.t_list)-1):

                proton.updateAcceleration(self.B)
                proton.updateRungeKutta(self.dt, proton, self.B)

                KE = proton.kineticEnergy()

                position_list = np.append(position_list, [proton.position], axis=0)
                KE_list = np.append(KE_list, KE)


            position_data = np.append(position_data, [position_list], axis = 0)         # Appends position_list to a new layer in position_data
            KE_data = np.append(KE_data, [KE_list], axis = 0)                           # Appends KE_list to a new column in KE_data
            

        np.save("RK_particle_positions.npy", position_data)
        np.save("RK_particle_KE_data.npy", KE_data)

                                                         


    

