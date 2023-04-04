from Particle import Particle
from Proton import Proton
from Electron import Electron
from particleBunch import particleBunch
import numpy as np
import matplotlib.pyplot as plt
import random
import copy
import scipy.constants


class EM_Cyclotron:

    """
    A class to represent a Cylcotron particle accelerator.

    Attributes:
    -----------
        dt: float
            The time step.
        t_max: float
            The total time that the simulation will be run over.
        B: 1D array, float entries.
           The magnetic field strength vector.
        E: 1D array, float entries.
           The electric field strength vector.

    Methods:
    --------
        runEulerSim(self, numOfPro, cyclotronRadius, E_boundary):
            Runs the simulation using the Euler-Forward approximation, saving positional and KE data to a .npy file.
        runCromerSim(self, numOfPro, cyclotronRadius, E_boundary):
            Runs the simulation using the Euler-Cromer approximation, saving positional and KE data to a .npy file.
        runKuttaSim(self, numOfPro, cyclotronRadius, E_boundary):
            Runs the simulation using the Runge-Kutta approximation, saving positional and KE data to a .npy file.
    """

    def __init__(self, dt, t_max, B, E):

        """
        Constructs all the necessary attributes for the Cyclotron object.

        Parameters:
        -----------
            dt : float
                 The time step.
            t_max : float
                The total time that the simulation will be run over.
            B : 1D array, float entries.
                The magnetic field strength vector.
            E : 1D array, float entries.
                The electric field strength vector.
        """

        self.B = B                                                      # Magnetic field strength.
        self.dt = dt                                                    # Time step.
        self.t_max = t_max                                              # Maximum time to reach.
        self.t_list = np.arange(0, self.t_max, self.dt, dtype=float)    # List from 0 - maximum time in steps of dt.
        self.E = E

    
    
    def runEulerSim(self, numOfPro, cyclotronRadius, E_boundary):
        """
        Runs the simulation using the Euler-Forward approximation, saving positional and KE data to a .npy file.

        Parameters:
        -----------
            numOfPro (int) : The number of particles in the simulation.
            cyclotronRadius (float) : The radial boundary of the cyclotron.
            E_boundary (float) : The boundary of the electric field in the y-axis.
        
        Returns:
        --------
            None

        """

        bunch = particleBunch.generateOneProton()
        #bunch = particleBunch.generateProtons(numOfPro, [0.001, 0.001, 0.001])                
        
        position_data = np.zeros((numOfPro, len(self.t_list), 3))        # Size = layers,columns,rows = numProtons, num iterations, 3 coordinates xyz
        KE_data = np.zeros((numOfPro, len(self.t_list)))                 # Size = columns,rows = numOfPro, no. of iterations

        w = (scipy.constants.elementary_charge * np.linalg.norm(self.B)) / scipy.constants.proton_mass

        # Starts simulation
        for proton in bunch:

            position_list = np.array([proton.position])                 # Array for proton positions
            KE_list = np.array([proton.kineticEnergy()])                # Array for proton KE's

            for j in range(0, len(self.t_list)-1):

                # Checks if inside cyclotron
                if proton.getRadius() < cyclotronRadius:
                    effective_B = self.B

                    # Checks if inside electric field
                    if proton.position[1] > -E_boundary and proton.position[1] < E_boundary:
                        effective_E = [0, self.E * np.sin(w * self.t_list[j]), 0]
                    
                    else:
                        effective_E = [0, 0, 0]

                else:
                    effective_B = [0, 0, 0]
                    effective_E = [0, 0, 0]

                proton.updateAcceleration(B = effective_B, E = effective_E)
                proton.updateEulerForward(dt = self.dt)

                KE = proton.kineticEnergy()

                position_list = np.append(position_list, [proton.position], axis=0)
                KE_list = np.append(KE_list, KE)


            position_data = np.append(position_data, [position_list], axis = 0)         # Appends particle specific positions to a new layer in position_data
            KE_data = np.append(KE_data, [KE_list], axis = 0)                           # Appends particle specific KE to a new column in KE_data
            
        np.save("EF_particle_positions.npy", position_data)
        np.save("EF_particle_KE_data.npy", KE_data)


    def runCromerSim(self, numOfPro, cyclotronRadius, E_boundary):
        """
        Runs the simulation using the Euler-Cromer approximation, saving positional and KE data to a .npy file.

        Parameters:
        -----------
            numOfPro (int) : The number of particles in the simulation.
            cyclotronRadius (float) : The radial boundary of the cyclotron.
            E_boundary (float) : The boundary of the electric field in the y-axis.
        
        Returns:
        --------
            None
        """

        bunch = particleBunch.generateOneProton()
        #bunch = particleBunch.generateProtons(numOfPro, [0.001, 0.001, 0.001])                
        
        position_data = np.zeros((numOfPro, len(self.t_list), 3))        # Size = layers,columns,rows = numProtons, num iterations, 3 coordinates xyz
        KE_data = np.zeros((numOfPro, len(self.t_list)))                 # Size = columns,rows = numOfPro, no. of iterations

        w = (scipy.constants.elementary_charge * np.linalg.norm(self.B)) / scipy.constants.proton_mass

        # Starts simulation
        for proton in bunch:

            position_list = np.array([proton.position])                 # Array for particle positions              
            KE_list = np.array([proton.kineticEnergy()])                # Array for particle KE           

            for j in range(0, len(self.t_list)-1):

                # Checks if inside cyclotron
                if proton.getRadius() < cyclotronRadius:
                    effective_B = self.B

                    # Checks if inside electric field
                    if proton.position[1] > -E_boundary and proton.position[1] < E_boundary:
                        effective_E = [0, self.E * np.sin(w * self.t_list[j]), 0]
                    
                    else:
                        effective_E = [0, 0, 0]

                else:
                    effective_B = [0, 0, 0]
                    effective_E = [0, 0, 0]

                proton.updateAcceleration(B = effective_B, E = effective_E)
                proton.updateEulerCromer(dt = self.dt)

                KE = proton.kineticEnergy()

                position_list = np.append(position_list, [proton.position], axis=0)
                KE_list = np.append(KE_list, KE)


            position_data = np.append(position_data, [position_list], axis = 0)         # Appends particle specific positions to a new layer in position_data
            KE_data = np.append(KE_data, [KE_list], axis = 0)                           # Appends particle specific KE to a new column in KE_data
            

        np.save("EC_particle_positions.npy", position_data)
        np.save("EC_particle_KE_data.npy", KE_data)


    def runKuttaSim(self, numOfPro, cyclotronRadius, E_boundary):
        """
        Runs the simulation using the Runge-Kutta approximation, saving positional and KE data to a .npy file.

        Parameters:
        -----------
            numOfPro (int) : The number of particles in the simulation.
            cyclotronRadius (float) : The radial boundary of the cyclotron.
            E_boundary (float) : The boundary of the electric field in the y-axis.
        
        Returns:
        --------
            None
        """

        bunch = particleBunch.generateOneProton()
        #bunch = particleBunch.generateProtons(numOfPro, [0.001, 0.001, 0.001])                
        
        position_data = np.zeros((numOfPro, len(self.t_list), 3))        # Size = layers,columns,rows = numProtons, num iterations, 3 coordinates xyz
        KE_data = np.zeros((numOfPro, len(self.t_list)))                 # Size = columns,rows = numOfPro, no. of iterations

        w = (scipy.constants.elementary_charge * np.linalg.norm(self.B)) / scipy.constants.proton_mass

        # Starts simulation
        for proton in bunch:

            position_list = np.array([proton.position])                     # Array for particle positions
            KE_list = np.array([proton.kineticEnergy()])                    # Array for particle KE


            for j in range(0, len(self.t_list)-1):

                # Checks if inside cyclotron
                if proton.getRadius() < cyclotronRadius:
                    effective_B = self.B

                    # Checks if inside electric field
                    if proton.position[1] > -E_boundary and proton.position[1] < E_boundary:
                        effective_E = [0, self.E * np.sin(w * self.t_list[j]), 0]
                    
                    else:
                        effective_E = [0, 0, 0]

                else:
                    effective_B = [0, 0, 0]
                    effective_E = [0, 0, 0]

                proton.updateAcceleration(B = effective_B, E = effective_E)
                proton.updateRungeKutta(dt = self.dt, particle = proton, B = effective_B, E = effective_E)

                KE = proton.kineticEnergy()

                position_list = np.append(position_list, [proton.position], axis=0)
                KE_list = np.append(KE_list, KE)


            position_data = np.append(position_data, [position_list], axis = 0)         # Appends particle specific positions to a new layer in position_data
            KE_data = np.append(KE_data, [KE_list], axis = 0)                           # Appends particle specific KE to a new column in KE_data


        np.save("RK_particle_positions.npy", position_data)
        np.save("RK_particle_KE_data.npy", KE_data)

                                                         


    

