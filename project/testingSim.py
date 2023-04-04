from Particle import Particle
from Simulation import Cyclotron
from cyclotronSimulation import EM_Cyclotron
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants



class gyrationTests:
    """
    Contains the methods related to testing the gyration of the particle in the simulaiton.

    Attributes:
    -----------
    None

    Methods:
    --------
    calculateData(dt, t_max, B, numOfPro):
        Runs the simulation for a constant magnetic field, using all numerical methods, to produce a fresh set of data.

    calculate_EM_Data(dt, t_max, B, numOfPro, cyclotronRadius, E, E_boundary):
        Runs the simulation for a complete cyclotron, using all numerical methods, to produce a fresh set of data.

    eulerPlot():
        Plots the gyration of the particle(s) in the x-y plane, using data stored from the Euler-Foward simulation.

    cromerPlot():
        Plots the gyration of the particle(s) in the x-y plane, using data stored from the Euler-Cromer simulation.

    kuttaPlot():
        Plots the gyration of the particle(s) in the x-y plane, using data stored from the Runge-Kutta simulation.

    comparePlots():
        Plots the gyrations from Euler-Forward, Euler-Cromer and Runge-Kutta methods onto the same axis.

    sinTest(dt, t_max, B):
        Plots the normalised x position of a particle against time.
        
    """

    def calculateData(dt, t_max, B, numOfPro):
        """
        Runs the simulation for a constant magnetic field, using all numerical methods, to produce a fresh set of data.

        Attributes:
        -----------
            dt (float) : The change in time / time step
            t_max (float) : Maximum time for the simulation to run
            numOfPro (int) : Number of particles in the simulation

        Returns:
        --------
            simulation (object) : The simulation object for a constant magnetic field
        
        """
        
        simulation = Cyclotron(dt=dt, t_max=t_max, B=B)
        simulation.runEulerSim(numOfPro = numOfPro)
        simulation.runCromerSim(numOfPro = numOfPro)
        simulation.runKuttaSim(numOfPro = numOfPro)

        print("All simulation versions were successfully executed, including data storage.")

        return simulation
    

    def calculate_EM_Data(dt, t_max, B, numOfPro, cyclotronRadius, E, E_boundary):
        """
        Runs the simulation for a complete cyclotron, using all numerical methods, to produce a fresh set of data.

        Attributes:
        -----------
            dt (float) : The change in time / time step
            t_max (float) : Maximum time for the simulation to run
            numOfPro (int) : Number of particles in the simulation
            cyclotronRadius (float) : Radius of the cyclotron bounds
            E (float) : Electric field strength vector
            E_boundary (float) : y-axis boundary which the E field is contained within

        Returns:
        --------
            simulation (object) : The simulation object for a complete cyclotron
        
        """
        
        simulation = EM_Cyclotron(dt=dt, t_max=t_max, B=B, E=E)
        
        simulation.runEulerSim(numOfPro = numOfPro, cyclotronRadius = cyclotronRadius, E_boundary = E_boundary)
        print("Euler Generated")

        simulation.runCromerSim(numOfPro = numOfPro, cyclotronRadius = cyclotronRadius, E_boundary = E_boundary)
        print("Cromer Generated")

        simulation.runKuttaSim(numOfPro = numOfPro, cyclotronRadius = cyclotronRadius, E_boundary = E_boundary)
        print("Kutta Generated")

        print("All simulation versions were successfully executed, including data storage.")

        return simulation


    def eulerPlot():
        """ Plots the gyration of the particle(s) in the x-y plane, using data stored from the Euler-Foward simulation. """

        positionData = np.load("EF_particle_positions.npy")

        for x in range(1, len(positionData)):
            plt.plot(positionData[x, :, 0], positionData[x, :, 1], color = "red")               

        plt.title("Euler-Forward Gyration Plot")
        plt.xlabel("X Position (m)")
        plt.ylabel("Y Position (m)")
        #plt.legend()
        plt.show()


    def cromerPlot():
        """ Plots the gyration of the particle(s) in the x-y plane, using data stored from the Euler-Cromer simulation. """

        positionData = np.load("EC_particle_positions.npy")

        for x in range(1, len(positionData)):
            plt.plot(positionData[x, :, 0], positionData[x, :, 1], color = "red")  

        plt.title("Euler-Cromer Gyration Plot")
        plt.xlabel("X Position (m)")
        plt.ylabel("Y Position (m)")
        #plt.legend()
        plt.show()


    def KuttaPlot():
        """ Plots the gyration of the particle(s) in the x-y plane, using data stored from the Runge-Kutta simulation. """

        positionData = np.load("RK_particle_positions.npy")

        for x in range(1, len(positionData)):
            plt.plot(positionData[x, :, 0], positionData[x, :, 1], color = "red")  

        plt.title("Runge-Kutta Gyration Plot")
        plt.xlabel("X Position (m)")
        plt.ylabel("Y Position (m)")
        #plt.legend()
        plt.show()

    def comparePlots():
        """ Plots the gyrations from Euler-Forward, Euler-Cromer and Runge-Kutta methods onto the same axis """

        eulerData = np.load("EF_particle_positions.npy")
        cromerData = np.load("EC_particle_positions.npy")
        kuttaData = np.load("RK_particle_positions.npy")

        labels = []
        for x in range(1, len(eulerData)):
            plt.plot(eulerData[x, :, 0], eulerData[x, :, 1], color = "Green") 
            plt.plot(cromerData[x, :, 0], cromerData[x, :, 1], color = "Blue")  
            plt.plot(kuttaData[x, :, 0], kuttaData[x, :, 1], color = "Purple")  

            if x == 1:  # append labels to list only once
                labels.append("Euler-Forward")
                labels.append("Euler-Cromer")
                labels.append("Runge-Kutta")


        plt.title("Simulation Proton Gyrations Using Three Numerical Approximations")
        plt.xlabel("X Position (m)")
        plt.ylabel("Y Position (m)")
        plt.legend(labels, numpoints=1)
        plt.show()

    # Testing conservation of energy via sin vs x position of one proton. Comparing accuracy of different methods to sin(x)
    def sinTest(dt, t_max, B):
        """
        Plots the normalised x position of a particle against time.
        Each numerical method is plotted along with y = sin(wt), where y is the analytical x position with a sinusoidal time dependency. 
        
        Parameters:
        -----------
            dt (float) : Change in time / time step
            t_max (float) : Max run time of the simulation
            B (1D array, float entries) : Magnetic field strength

        Returns:
        --------
            None

        """
        
        t_list = np.arange(0, t_max, dt, dtype=float)    # List from 0 - maximum time in steps of dt.
        
        eulerData = np.load("EF_particle_positions.npy")
        cromerData = np.load("EC_particle_positions.npy")
        KuttaData = np.load("RK_particle_positions.npy")

        # Normalising x positions between -1 and 1 to fit with sin function
        normalized_euler = 2*(eulerData[1, :, 0] - eulerData[1, :, 0].min())/(eulerData[1, :, 0].max() - eulerData[1, :, 0].min()) - 1
        normalized_cromer = 2*(cromerData[1, :, 0] - cromerData[1, :, 0].min())/(cromerData[1, :, 0].max() - cromerData[1, :, 0].min()) - 1
        normalized_kutta = 2*(KuttaData[1, :, 0] - KuttaData[1, :, 0].min())/(KuttaData[1, :, 0].max() - KuttaData[1, :, 0].min()) - 1

        B_mag = np.linalg.norm(B)
        omega = scipy.constants.elementary_charge * B_mag / scipy.constants.proton_mass

        # Producing the array of sin values for each time step
        y = -np.cos(omega*t_list)   

        plt.plot(t_list, normalized_euler, color = "Green", label = "Euler-Forward")            
        # plt.plot(t_list, normalized_cromer, color = "Blue", label = "Euler-Cromer")        
        # plt.plot(t_list, normalized_kutta, color = "Purple", label = "Runge-Kutta")           
        # plt.plot(t_list, y, color = "Black", label = "sin(wt)")           


        plt.title("x-positions against time for all numerical methods")
        plt.xlabel("Time (s)")
        plt.ylabel("X position (m)")
        plt.legend()
        plt.show()




