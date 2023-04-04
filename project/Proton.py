from Particle import Particle
import scipy.constants


class Proton(Particle):
    """
    A subclass of the Particle class that represents a proton. 

    
    Attributes:
    -----------
        position : inhetited from Particle class
        velocity : inhetited from Particle class
        acceleration : inhetited from Particle class
        mass : float
               Mass of a proton
        charge : float
                 The charge of a proton.


    Methods:
    --------
        super().__init__ : inherits the attributes of Particle class while adding proton mass and charge

    """

    def __init__(self, position, velocity, acceleration):
        super().__init__(position, velocity, acceleration, 
                         mass = scipy.constants.proton_mass, 
                         charge = scipy.constants.elementary_charge, 
                         name = "Proton"
                         )
