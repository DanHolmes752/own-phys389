from Particle import Particle
from Proton import Proton
from Electron import Electron
import numpy as np
import matplotlib.pyplot as plt
import random
import copy
import scipy.constants


class particleBunch:
    """ Contains the methods related to initialising the particle(s) used in the simulation.
     
    Attributes
    ----------
    None
    
    Methods:
    --------
    generateOneProton():
        Initialises one proton object using the Proton subclass.

    generateProtons():
        Creates an array containing a specificed number of randomly initialised protons.

    generateElectrons():
        Creates an array containing a specificed number of randomly initialised electrons.

    """
    
    def generateOneProton():
        """
        Initialses one proton object using the Proton subclass, with specified position and velocity vectors.

        Parameters:
        -----------
            None
        
        Returns:
        --------
            protons (1D array) : A list containing the proton object
                                 Returned in an array for continuity with multiple particle simulations 
        
        """

        protons = []

        position = np.array( [0, 0, 0] )          
        velocity = np.array( [ 0, 1e4, 0 ] )        # Needs to be in y direction for B in z direction

        proton = Proton(position, velocity, np.array([[0, 0, 0]]))
        protons.append(proton)

        return protons



    def generateProtons(numOfPro, spawn_bounds):

        """
        Creates a protons array containing a specificed number of randomly initialised protons.

        Parameters:
        -----------
            numOfPro (Int) : The number of protons to be initialised
            spawn_bounds (1D array, float entries) : The boundaries of where the proton can be initialised
        
        Returns:
        --------
            protons (1D array) : Array of randomly initialised proton objects

        """

        protons = []
        for x in range(numOfPro):

            # Randomly assigns position and velocity within certain bounds
            position = np.array( [ random.uniform(-spawn_bounds[1], spawn_bounds[1]),
                                   random.uniform(-spawn_bounds[1], spawn_bounds[1]), 
                                   random.uniform(-spawn_bounds[1], spawn_bounds[1]) ] )          
            velocity = np.array( [ 0, random.uniform(1e2, 1e4), 0 ] )

            proton = Proton(position, velocity, np.array([[0, 0, 0]]))
            protons.append(proton)

        return protons
    


    def generateElectrons(numOfPro, spawn_bounds):
        """
        Creates an array containing a specificed number of randomly initialised electrons.

        Parameters:
        -----------
            numOfPro (Int) : The number of electrons to be initialised
            spawn_bounds (1D array, float entries) : The boundaries of where the electron can be initialised
        
        Returns:
        --------
            electrons (1D array) : Array of randomly initialised electron objects

        """
        

        electrons = []
        for x in range(numOfPro):

            # Randomly assigns position and velocity within certain bounds
            position = np.array( [ random.uniform(-spawn_bounds[1], spawn_bounds[1]), random.uniform(-spawn_bounds[1], spawn_bounds[1]), random.uniform(-spawn_bounds[1], spawn_bounds[1]) ] )          # can use spawn bounds here if you want
            velocity = np.array( [ 0, random.uniform(1e2, 1e4), 0 ] )

            electron = Electron(position, velocity, np.array([[0, 0, 0]]))
            electrons.append(electron)

        return electrons