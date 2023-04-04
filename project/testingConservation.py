from Particle import Particle
from Simulation import Cyclotron
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants
from testingSim import gyrationTests

class energyConservation:
    """
    Contains the methods related to testing energy conservation.

    Attributes:
    -----------
        None

    Methods:
    --------
        calcAngularFrequency(x_positions, t_list): 
            Returns the angular frequency of each gyration of the particle.

        eulerAngularFrequency(dt, t_max):   
            Plots the angular frequency of each gyration against the number of gyrations of the particle using data from the Euler-Forward simulation.

        cromerAngularFrequency(dt, t_max):
            Plots the angular frequency of each gyration against the number of gyrations of the particle using data from the Euler-Cromer simulation.

        kuttaAngularFrequency(dt, t_max):
            Plots the angular frequency of each gyration against the number of gyrations of the particle using data from the Runge-Kutta simulation.

        compareAngular(dt, t_max):
            Plots the angular frequencies from each numerical method onto the same axis

        eulerKineticTest(dt, t_max):
            Plots the saved kinetic energy data, from the Euler-Forward simulation, against time.

        cromerKineticTest(dt, t_max):
            Plots the saved kinetic energy data, from the Euler-Cromer simulation, against time.

        kuttaKineticTest(dt, t_max):
            Plots the saved kinetic energy data, from the Runge-Kutta simulation, against time.

        compareKineticTests(dt, t_max):
            Plots the saved kinetic energy data, from all numerical methods, against time onto the same axis.

    """

    
    def calcAngularFrequency(x_positions, t_list):
        """
        Returns the angular frequency of each gyration of the particle.
        It also calculates the relative change in radius and prints it to the terminal.

        Parameters:
        -----------
            x_positions (1D array, float entries) : The array of x positions for each time step
            t_list (1D array, float entries) : The array of time steps

        Returns:
        --------   
            numOfGyrations (int) : The number of gyrations completed by the particle
            w (1D array, float entries) : The array of angular frequencies for each gyration
        
        """

        peaks = []
        w = []
        radii_points = []

        for i in range(1, len(t_list)-1):
            if x_positions[i] >= [x_positions[i-1]] and x_positions[i] >= [x_positions[i+1]]:
                peaks.append(t_list[i])
                radii_points.append(x_positions[i])                # find min and max to get difference -> error = diff/min

        for j in range(0, len(peaks)-1):
            T = peaks[j+1] - peaks[j]
            w.append(2*np.pi / T)                                  # calculating angular frequency and appending to array

        numOfGyrations =  [i for i in range(1, len(w)+1)]    # Finds number of gyrations of the particle via length of period array

        print("\nRelative change in angular frequency: ", (w[-1] - w[1]) / w[1])
        print("Relative change in radius: ", (radii_points[-1] - radii_points[1]) / radii_points[1])

        print(w)
        return numOfGyrations, w


    def eulerAngularFrequency(dt, t_max):
        """
        Plots the angular frequency of each gyration against the number of gyrations of the particle
        using data from the Euler-Forward simulation.

        Parameters:
        -----------
            dt (float) : The change in time / time step
            t_max (float) : The maximum time the simulation runs for

        Returns:
        --------   
            None
        
        """
        
        eulerData = np.load("EF_particle_positions.npy")
        t_list = np.arange(0, t_max, dt, dtype=float)       # List from 0 - maximum time in steps of dt.

        #print(eulerData[1,:,0])

        # for x in range(len(eulerData)):
        plotInfo = energyConservation.calcAngularFrequency(x_positions = eulerData[1, :, 0], t_list = t_list)    #
        plt.plot(plotInfo[0], plotInfo[1], color = "Blue")                                                        #


        plt.title("Euler-Forward Particle Angular Frequency")
        plt.xlabel("Number of Gyrations")
        plt.ylabel("Angular Frequency (Rad/s)")
        plt.show()      


    def cromerAngularFrequency(dt, t_max):
        """
        Plots the angular frequency, of each gyration, against the number of gyrations of the particle 
        using data from the Euler-Cromer simulation.

        Parameters:
        -----------
            dt (float) : The change in time / time step
            t_max (float) : The maximum time the simulation runs for

        Returns:
        --------   
            None
        
        """

        cromerData = np.load("EC_particle_positions.npy")
        t_list = np.arange(0, t_max, dt, dtype=float)        # List from 0 - maximum time in steps of dt.

        plotInfo = energyConservation.calcAngularFrequency(x_positions = cromerData[1, :, 0], t_list = t_list)   

        plt.plot(plotInfo[0], plotInfo[1], color = "Blue")  
        plt.title("Euler-Cromer Particle Angular Frequency")
        plt.xlabel("Number of Gyrations")
        plt.ylabel("Angular Frequency (Rad/s)")
        plt.show() 


    def kuttaAngularFrequency(dt, t_max):
        """
        Plots the angular frequency, of each gyration, against the number of gyrations of the particle 
        using data from the Runge-Kutta simulation.

        Parameters:
        -----------
            dt (float) : The change in time / time step
            t_max (float) : The maximum time the simulation runs for

        Returns:
        --------   
            None
        
        """

        kuttaData = np.load("RK_particle_positions.npy")
        t_list = np.arange(0, t_max, dt, dtype=float)    # List from 0 - maximum time in steps of dt.

        plotInfo = energyConservation.calcAngularFrequency(x_positions = kuttaData[1, :, 0], t_list = t_list)

        plt.plot(plotInfo[0], plotInfo[1], color = "Blue")  
        plt.title("Runge-Kutta Particle Angular Frequency")
        plt.xlabel("Number of Gyrations")
        plt.ylabel("Angular Frequency")
        plt.show() 


    def compareAngular(dt, t_max):
        """
        Plots the angular frequencies from each numerical method onto the same axis

        Parameters:
        -----------
            dt (float) : The change in time / time step
            t_max (float) : The maximum time the simulation runs for

        Returns:
        --------   
            None
        
        """

        eulerData = np.load("EF_particle_positions.npy")
        cromerData = np.load("EC_particle_positions.npy")
        kuttaData = np.load("RK_particle_positions.npy")

        t_list = np.arange(0, t_max, dt, dtype=float)    # List from 0 - maximum time in steps of dt.

        eulerFreq = energyConservation.calcAngularFrequency(x_positions = eulerData[1, :, 0], t_list = t_list)
        cromerFreq = energyConservation.calcAngularFrequency(x_positions = cromerData[1, :, 0], t_list = t_list)
        kuttaFreq = energyConservation.calcAngularFrequency(x_positions = kuttaData[1, :, 0], t_list = t_list)

        print("euler w: ", eulerFreq[0])
        print("cromer w: ", cromerFreq[0])
        print("kutta w: ", kuttaFreq[0])

        plt.plot(eulerFreq[0], eulerFreq[1], color = "Green", label = "Euler-Forward")  
        plt.plot(cromerFreq[0], cromerFreq[1], color = "Blue", label = "Euler-Cromer")  
        plt.plot(kuttaFreq[0], kuttaFreq[1], color = "Purple", label = "Runge-Kutta")  

        plt.title("Particle Angular Frequency For Each Gyration")
        plt.xlabel("Number of Gyrations")
        plt.ylabel("Angular Frequency (Rad/s)")
        plt.legend()
        plt.show() 

    
    def eulerKineticTest(dt, t_max):
        """
        Plots the saved kinetic energy data, from the Euler-Forward simulation, against time.

        Parameters:
        -----------
            dt (float) : The change in time / time step
            t_max (float) : The maximum time the simulation runs for

        Returns:
        --------   
            None
        
        """

        KE = np.load("EF_particle_KE_data.npy")
        t_list = np.arange(0, t_max, dt, dtype=float)    # List from 0 - maximum time in steps of dt.

        for x in range(1, len(KE)):
            plt.plot(t_list, KE[x], color = "Green")  

        plt.title("Euler-Forward Particle Kinetic Energy")
        plt.xlabel("Time (s)")
        plt.ylabel("Kinetic Energy (J)")
        plt.show() 

    
    def cromerKineticTest(dt, t_max):
        """
        Plots the saved kinetic energy data, from the Euler-Cromer simulation, against time.

        Parameters:
        -----------
            dt (float) : The change in time / time step
            t_max (float) : The maximum time the simulation runs for

        Returns:
        --------   
            None
        
        """

        KE = np.load("EC_particle_KE_data.npy")
        t_list = np.arange(0, t_max, dt, dtype=float)    # List from 0 - maximum time in steps of dt.

        for x in range(1, len(KE)):
            plt.plot(t_list, KE[x], color = "Green")  

        plt.title("Euler-Cromer Particle Kinetic Energy")
        plt.xlabel("Time (s)")
        plt.ylabel("Kinetic Energy (J)")
        plt.show() 


    def kuttaKineticTest(dt, t_max):
        """
        Plots the saved kinetic energy data, from the Runge-Kutta simulation, against time.

        Parameters:
        -----------
            dt (float) : The change in time / time step
            t_max (float) : The maximum time the simulation runs for

        Returns:
        --------   
            None
        
        """
         
        KE = np.load("RK_particle_KE_data.npy")
        t_list = np.arange(0, t_max, dt, dtype=float)    # List from 0 - maximum time in steps of dt.

        for x in range(1, len(KE)):
            plt.plot(t_list, KE[x], color = "Green")  

        plt.title("Runge-Kutta Particle Kinetic Energy")
        plt.xlabel("Time (s)")
        plt.ylabel("Kinetic Energy (J)")
        plt.show() 


    def compareKineticTests(dt, t_max):
        """
        Plots the saved kinetic energy data, from all numerical methods, against time onto the same axis.
        It also prints the relative change of the kinetic energy for each numerical method into the terminal.

        Parameters:
        -----------
            dt (float) : The change in time / time step
            t_max (float) : The maximum time the simulation runs for

        Returns:
        --------   
            None
        
        """
         
        EF_KE1 = np.load("EF_particle_KE_data.npy")
        EC_KE1 = np.load("EC_particle_KE_data.npy")
        RK_KE1 = np.load("RK_particle_KE_data.npy")

        # Calculating and printing the relative changes of kinetic energy for each numerical method
        # print((EF_KE[1, -1]-EF_KE[1, 1]) / EF_KE[1, 1])
        # print((EC_KE[1, -1]-EC_KE[1, 1]) / EC_KE[1, 1])
        # print((RK_KE[1, -1]-RK_KE[1, 1]) / RK_KE[1, 1])

        EF_KE = EF_KE1 / scipy.constants.elementary_charge
        EC_KE = EC_KE1 / scipy.constants.elementary_charge
        RK_KE = RK_KE1 / scipy.constants.elementary_charge

        print("EF last ke: ", EF_KE[1])
        print("EC last ke: ", EC_KE[1])
        print("RK last ke: ", RK_KE[1])


        t_list = np.arange(0, t_max, dt, dtype=float)    # List from 0 - maximum time in steps of dt.

        for x in range(1, len(EF_KE)):
            plt.plot(t_list, EF_KE[x], color = "Green", label = "Euler-Forward")  

            plt.plot(t_list, EC_KE[x], color = "Blue", label = "Euler-Cromer")  

            plt.plot(t_list, RK_KE[x], color = "Purple", label = "Runge-Kutta")  

        plt.title("Simulated Particle Kinetic Energy Against Time")
        plt.xlabel("Time (s)")
        plt.ylabel("Kinetic Energy (eV)")
        plt.legend()
        plt.show() 