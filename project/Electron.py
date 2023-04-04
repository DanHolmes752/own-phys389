from Particle import Particle
import scipy.constants


class Electron(Particle):
    """
    A subclass of the Particle class that represents an electron. 

    
    Attributes:
    -----------
        position : inhetited from Particle class
        velocity : inhetited from Particle class
        acceleration : inhetited from Particle class
        mass : float
               Mass of an electron
        charge : float
                 The charge of an electron

                 
    Methods:
    --------
        super().__init__ : inherits the attributes of Particle class while adding electron mass and charge

    """
     
    def __init__(self, position, velocity, acceleration):
        super().__init__(position, velocity, acceleration, 
                         mass = scipy.constants.electron_mass, 
                         charge = -scipy.constants.elementary_charge, 
                         name = "Electron"
                         )
