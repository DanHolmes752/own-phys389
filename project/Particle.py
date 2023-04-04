import numpy as np
import matplotlib as plt
import copy

class Particle:
    """
    A class to initialise and represent a charged particle.

    Attributes:
    -----------
    position :  1D array, float entries.
        Initial position vector of the particle.
    velocity : 1D array, float entires.
        Initial velocity vector of the particle.
    acceleration : 1D array, float entires.
        Initial acceleration vector of the particle.
    mass : float
        Mass of the particle.
    charge : float
        Charge of the particle.
    name : str
        The name of the particle.

    Methods:
    --------
    updateAcceleration(B): 
        Using the Lorentz Force equation to update the particle's acceleration vector.
    updateEulerForward(dt): 
        Updates the particles, position and velocity vectors using the Euler Forward integration method.
    updateEulerCromer(dt): 
        Updates the particles, velocity and position vectors using the Euler-Cromer integration method.
    kineticEnergy(): 
        Calculates and returns the particle's kinetic energy for it's current velocity.

    """

    def __init__(self, position, velocity, acceleration, mass, charge, name):
        """
        Constructs all the necessary attributes for the particle object.

        Parameters:
        -----------
        position :  1D array, float entries.
            Initial position vector of the particle.
        velocity : 1D array, float entires.
            Initial velocity vector of the particle.
        acceleration : 1D array, float entires.
            Initial acceleration vector of the particle.
        mass : float
            Mass of the particle.
        charge : float
            Charge of the particle.
        name : str
            The name of the particle.
        """


        self.position = np.array(position, dtype = float)
        self.velocity = np.array(velocity, dtype = float)
        self.acceleration = acceleration
        self.mass = mass
        self.charge = charge
        self.name = name
   


    def updateAcceleration(self, B, E):
        """
        Updates a particle's acceleration, calculated from the Lorentz force and F=ma.
        
        Parameters:
        -----------
            B : 1D Array, float entries.
                The magnetic field strength vector. 
            E : 1D array, float entries.
                The electric field strength vector.

        Returns:
        --------
            self.acceleration (1D array, float entries) : The updated acceleration vector.
        
        """

        # Calculate the magnetic force and acceleration
        Force = self.charge * (E + np.cross(self.velocity, B))
        self.acceleration = Force / self.mass

        return self.acceleration


    def updateEulerForward(self, dt):
        """
        Updates the particle's position and velocity vectors, via the Euler Forward integration method.
        
        Parameters:
        -----------
            dt : Float
                 A value for the change in time / time step.
        
        Returns:
        --------
            None

        """
        
        self.position = self.position + self.velocity * dt            
        self.velocity = self.velocity + self.acceleration * dt        
    

    def updateEulerCromer(self, dt):    
        """
        Updates the particle's position and velocity vectors, via the Euler-Cromer method. 
        
        Parameters:
        -----------
            dt : Float
                 A value for the change in time / time step.

        Returns:
        --------
            None

        """

        self.velocity = self.velocity + self.acceleration * dt
        self.position = self.position + self.velocity * dt


    def updateRungeKutta(self, dt, particle, B, E):
        """
        Updates the particle's position and velocity vectors, via the second order Runge-Kutta method. 
        
        Parameters:
        -----------
            dt : Float
                 A value for the change in time
            particle : object
                The particle object from the sub class that inherits from Particle
            B : 1D numpy array, float entries
                The magnetfic field strength vector
            E : 1D numpy array, float entries
                The electric field strength vector

        Returns:
        --------
            None

        """
        # Creates a copy of the particle to produce dt/2 data
        rkParticle = copy.copy(particle)     

        # Finds velocity and position for a dt/2 timestep, stores in copied particle object
        rkParticle.velocity = rkParticle.velocity + (dt/2) * rkParticle.updateAcceleration(B = B, E = E)
        rkParticle.position = rkParticle.position + (dt/2) * rkParticle.velocity

        # Updates actual particle velocity and position using acceleration calculated from dt/2 estimate
        self.velocity = self.velocity + dt*rkParticle.updateAcceleration(B = B, E = E) 
        self.position = self.position + dt*rkParticle.velocity
    
    
    def kineticEnergy(self):
        """ Returns the particles kinetic energy """

        v_mag_squared = (self.velocity[0]**2 + self.velocity[1]**2 + self.velocity[2]**2) 
        KE = 0.5 * self.mass * v_mag_squared
        return KE
    

    def getRadius(self):
        """ Returns the current radius magnitude of the particle """
        
        radius = np.sqrt(self.position[0]**2 + self.position[1]**2 + self.position[2]**2)
        return radius


    


