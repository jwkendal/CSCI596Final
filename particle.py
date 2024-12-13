import numpy as np

# Global Constants
global Tb, Lb, R, g, m
Tb = 294.261                # Temperature at Sea Level [K] (70 F)
Lb = -0.0065                # Temperature lapse rate [K/m]
R = 8.31432                 # Universal gas constant [N*m/mol*K]
g = 9.81                    # Gravity Constant [m/s^2]
n = 6.022*(10**23)          # Number of molecules in  mole
m = 4.799*(10**(-26))       # Mass of an air molecule [kg]
r = 1.8*(10**(-10))         # Radius of an air molecule [m]

class Particle:
    # Constructor
    def __init__(self, pid, x, y, z, vx, vy, vz, mass, temp, ke):
        # Set Properties
        # Particle ID Number
        self.pid = pid

        # Position (m)
        self.r = np.array([x, y, z])
        
        # Velocity (m/s)
        self.v = np.array([vx, vy, vz])
        
        # Acceleration (m/s^2) - set to 0
        #self.ax = 0
        #self.ay = 0
        #self.az = 0
        
        # Mass (kg)
        self.mass = mass

        # Temperature (K)
        self.temperature = temp

        # Kinetic Energy (J)
        self.kineticEnergy = ke
        
        # Collision flag to mark whether this particle has collided this iteration
        self.collided = False
        
        # List of particles it has collided with this iteration
        self.collisions = [-1]
        
        # List of particles this particle should check for a collision with this iteration
        # Any particle it has already been checked against will be removed
        self.toCheck = []
        
        # Temperature of molecule
    # end constructor
    
    # Calculate magnitude of velocity
    def getVelMag(self):
        vel = np.linalg.norm(self.v)
        return vel
    # end getVelMag
    
    # Calculate kinetic energy of particle
    def calcKinEnergy(self):
        vel = self.getVelMag()
        self.kineticEnergy = 0.5 * self.mass * (vel**2)
        
        # Update temperature
        self.calcTemperature()
    # end calcKinEnergy
    
    # Calculate temperature of a particle
    def calcTemperature(self):
        self.temperature = (2 * self.kineticEnergy * n) / (3 * R)
    # end calcTemperature
    
    # Add collision
    def addCollision(self, collisionPID):
        self.collided = True
        if -1 in self.collisions:
            self.collisions = []
        # end if
        self.collisions.append(collisionPID)
        self.toCheck.remove(collisionPID)
    # end addCollision
    
    # Clear collisions
    def clearCollisions(self, numParticles):
        self.collided = False
        self.collisions = [-1]
        self.toCheck = list(range(numParticles))
    # end clearCollisions
    
# end Particle Class