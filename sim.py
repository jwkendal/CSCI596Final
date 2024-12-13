import os
from particle import Particle
import numpy as np
from random import random as rand
from pandas import DataFrame
from pandas import concat

# Global Constants
global Tb, Lb, R, g, m
Tb = 294.261                # Temperature at Sea Level [K] (70 F)
Lb = -0.0065                # Temperature lapse rate [K/m]
R = 8.31432                 # Universal gas constant [N*m/mol*K]
g = 9.81                    # Gravity Constant [m/s^2]
n = 6.022*(10**23)          # Number of molecules in  mole
m = 4.799*(10**(-26))       # Mass of an air molecule [kg]
r = 1.8*(10**(-10))         # Radius of an air molecule [m]
d = 2*r

class Simulation:
    # Constructor
    def __init__(self, numParticles, size, timeStep, endTime, tempPct):
        # Set the number of particles to use in this simulation
        self.numParticles = numParticles
        
        # Set the size of the cube this simulation occurs in
        self.cubeSize = size
        
        # Set the time step between each update
        self.dt = timeStep
        self.time = 0
        self.endTime = endTime
        
        # Set the percentages and temperatures of molecules
        self.tempPct = tempPct
        
        # Create a list of Particles
        self.particles = [None] * self.numParticles
        
        # Calculate the starting positions/temperatures/velocities of each particle
        self.initParticles()
    # end constructor
    
    # Calculate the number of particles to use
    def calcNumParticles(self):
        # mols = (rho*V)/m
        self.numParticles = ((self.density*self.length*self.width*self.height)/m)*n
    # end calcNumParticles
    
    # Calculate distance between particles
    def calcParticleDist(self, particle1:Particle, particle2:Particle):
        # Get the distance vector between particles
        dist = np.subtract(particle1.r,particle2.r)
        
        # Calculate the magnitude of the distance vector
        distMag = np.linalg.norm(dist)
        
        return distMag
    # end calcParticleDist
    
    # Get normalize velocity vector
    def getNormVel(self):
        # Randomly generate 3 values
        x = rand()
        y = rand()
        z = rand()

        # Calculate the magnitude of the velocity vector
        randVel = np.array([x,y,z])
        mag = np.linalg.norm(randVel)

        # Get the normalized values of each component
        normVel = np.divide(randVel, mag)

        # Set some velocity components to negative
        for idx in range(3):
            val = rand()
            if val < 0.5:
                normVel[idx] = -1 * normVel[idx]
            # end if
        # end for

        return normVel
    # end getNormVel

    # Calculate starting positions/velocities of particles
    def initParticles(self):
        # Get set of possible starting coordinates in cube
        numPos = int(np.ceil(self.numParticles / 9))
        if numPos < 3:
            numPos = 3
        # end if
        startPos1D = np.linspace(0, self.cubeSize, numPos)
        
        # Make 3d array of possible starting points in matlab
        self.startPos = [None] * numPos**3
        for x in range(numPos):
            for y in range(numPos):
                for z in range(numPos):
                    idx = 9*x + 3*y + z
                    entry = [startPos1D[x], startPos1D[y], startPos1D[z]]
                    self.startPos[idx] = entry
                # end for
            # end for 
        # end for
        
        # Figure out which particles will be initialized to each temperature
        idx = 0
        usedParticles = 0
        self.startTemp = [None] * self.numParticles
        self.startKE = self.startTemp.copy()
        self.startVel = self.startTemp.copy()
        self.tempKE = self.tempPct.copy()
        for temp in self.tempPct:
            # Calculate the kinetic energy for this temperature
            kinEnergy = (3/2) * (R/n) * temp
            self.tempKE[temp] = kinEnergy

            # Calculate the velocity at this temperature
            vel = np.sqrt((3 * (R / n) * temp) / m)
            
            # Get the percentage of particles to set to this temp
            pct = self.tempPct[temp]*0.01
            
            # Get the number of particles based on percentage
            if idx == 0:
                numParticles = int(np.ceil(self.numParticles * pct))
            elif idx < len(self.tempPct.keys()) - 1:
                numParticles = int(np.floor(self.numParticles * pct))
            else:
                numParticles = int(self.numParticles - usedParticles)
            # end if/elseif/else
            
            for partIdx in range(usedParticles,usedParticles+numParticles):
                self.startTemp[partIdx] = temp
                self.startKE[partIdx] = kinEnergy
                self.startVel[partIdx] = vel
            # end for
            
            # Increment index and used particles count
            idx += 1
            usedParticles += numParticles
        # end for
        
        # Loop over every particle in the simulation
        for particleIdx in range(self.numParticles):
            # Calculate the normalized velocity vector
            velVec = [None] * 3
            idx = 0
            for velPart in self.getNormVel():
                velVec[idx] = self.startVel[particleIdx] * velPart
                idx += 1
            # end for

            # Make a new particle object
            newParticle = Particle(particleIdx, self.startPos[particleIdx][0], self.startPos[particleIdx][1], 
                                   self.startPos[particleIdx][2], velVec[0], velVec[1], velVec[2], m, self.startTemp[particleIdx], self.startKE[particleIdx])
            
            # Initialize list of PIDs it should check for collisions with
            newParticle.toCheck = list(range(self.numParticles))
            newParticle.toCheck.remove(newParticle.pid)
            
            # Add particle to list
            self.particles[particleIdx] = newParticle            
            
        # end for loop
    # end initParticles
    
    # Update particle position
    def updatePos(self, particle: Particle):
        # Get the current position and velocity
        pos = particle.r
        vel = particle.v
        
        outOfBounds = False
        
        # Calculate the new position
        posChange = np.multiply(vel, self.dt)
        newPos = np.add(pos, posChange)
        particle.r = newPos.copy()
       
        # If out of bounds, set flag
        OOBCheck = np.greater(newPos, self.cubeSize)
        OOBCheck2 = np.less(newPos, 0)
        check = np.logical_or(OOBCheck, OOBCheck2)
        indices = []
        if np.any(check):
            outOfBounds = True
            indices = np.where(check)
        # end if

        # Check if the new position is out of bounds
        if outOfBounds:
            self.handleBoundaryCondition(particle, indices[0])
        # end if
        
        if True:
            pass
    # end updatePos
    
    # Handle boundary condition
    def handleBoundaryCondition(self, particle: Particle, indices):
        # Reverse velocity on boundary dimension
        for idx in indices:
            if (particle.r[idx] > self.cubeSize):
                particle.r[idx] = self.cubeSize
                particle.v[idx] = -1 * particle.v[idx]
            # end if
            elif (particle.r[idx] < 0):
                particle.r[idx] = 0
                particle.v[idx] = -1 * particle.v[idx]
            # end elif
        # end for
    # end handleBoundaryCondition
    
    # Check for Particle Collisions
    def collisionCheck(self):
        # Loop over each particle
        for particle in self.particles:
            # Loop over all other particles
            for particle2PID in particle.toCheck:
                # Get particle
                particle2 = self.particles[particle2PID]
                
                # If same particle, skip
                if particle == particle2:
                    continue
                # Check if these two particles have already been collided with
                elif particle.collided and (particle2.pid in particle.collisions):
                    continue                
                # end if
                
                quickDist = np.abs(np.subtract(particle.r, particle2.r))
                check = np.greater(quickDist, d)
                # Do a quick check to see if these particles are nowhere near each other
                if np.any(check):
                    particle.toCheck.remove(particle2PID)
                    particle2.toCheck.remove(particle.pid)
                    continue
                # end if
                
                # Calculate the distance between these particles
                dist = self.calcParticleDist(particle, particle2)
                
                # If this distance is less than the 2*radius
                if dist <= d:                    
                    # Handle Particle collision
                    self.handleParticleCollision(particle, particle2)
                    
                    # Mark these particles as collided
                    particle.addCollision(particle2.pid)
                    particle2.addCollision(particle.pid)
                else:
                    particle.toCheck.remove(particle2PID)
                    particle2.toCheck.remove(particle.pid)
                # end if
        # end for
    # end
   
    # Calculate particle velocities post-collision
    def calcPostCollisionVel(self, particle1: Particle, particle2: Particle):
        # Calculate the normal vector between the particles
        distMag = 0
        normal = np.subtract(particle1.r, particle2.r)
        distMag = np.linalg.norm(normal)
        
        #if distMag < d:
        #    distMag = d
        # end if
        
        normal = np.divide(normal, distMag)
        
        # Calculate the relative velocity
        v_relative = np.subtract(particle1.v, particle2.v)
        
        # Calculate the relative velocity along the normal direction
        dotProd = np.dot(v_relative, normal)
        v_normal = np.multiply(normal, dotProd)
        
        # Calculate the new velocities of each particle
        particle1.v = np.subtract(particle1.v, v_normal)
        particle2.v = np.add(particle2.v, v_normal)
        
    # end calcPostCollisionVel
    
    # Particle Collision resolution
    def handleParticleCollision(self, particle1: Particle, particle2: Particle):
        # Save original total kinetic energy
        totalKE = particle1.kineticEnergy + particle2.kineticEnergy
        
        # Calculate new velocity for each particle
        self.calcPostCollisionVel(particle1, particle2)
        
        # Update kinetic energy/temperature for both particles
        particle1.calcKinEnergy()
        particle2.calcKinEnergy()
        
        # Assert that total kinetic energy was conserved
        assert(np.absolute(((particle1.kineticEnergy + particle2.kineticEnergy) - totalKE)) < 1*(10**-25))
    # end handleParticleCollision
   
    def simParticles(self):
        # Loop over molecule to update position
        for particle in self.particles:
            self.updatePos(particle)
        # end while
        
        # Check for particle collisions and handle them if necessary
        self.collisionCheck()
    # end simParticles
   
    # Run main simulation
    def runSimulation(self, table: DataFrame):
        self.time += self.dt

        while self.time <= self.endTime:
            # Simulate particles
            self.simParticles()
            
            # Update simulation table with new data
            table = self.updateSimTable(table)
            
            # Clear collision data from particles
            for particle in self.particles:
                particle.clearCollisions(self.numParticles)
            # end for
            
            # Update time
            self.time += self.dt
        # end while
        
        return table
    # end runSimulation
    
    # Add latest data to dataframe
    def updateSimTable(self, table: DataFrame):
        # Loop over molecules
        for particleIdx in range(self.numParticles):
            # Get particle
            curPar = self.particles[particleIdx]

            # Construct new table entry
            newEntry =     {"Time [s]": [self.time],
                            "PID": [int(curPar.pid)],
                            "Position-X [m]": [curPar.r[0]],
                            "Position-Y [m]": [curPar.r[1]],
                            "Position-Z [m]": [curPar.r[2]],
                            "Velocity-X [m]": [curPar.v[0]],
                            "Velocity-Y [m]": [curPar.v[1]],
                            "Velocity-Z [m]": [curPar.v[2]],
                            "Temperature [K]": [curPar.temperature],
                            "Kinetic Energy [J]": [curPar.kineticEnergy],
                            "Collision": [curPar.collided],
                            "Collided Particles": curPar.collisions}
            
            # Add entry to table
            table = concat([table, DataFrame(newEntry)])
        # end for

        return table
    # end updateSimTable

    # Write data to xyz file
    def writeXYZ(self, table: DataFrame, xyzFile):
        # Loop over time steps
        #for curTime in np.arange(0,self.endTime,self.dt):
        curTime = 0
        while curTime < self.endTime:
            # Find all entries from dataframe at this time
            rows = table.loc[table["Time [s]"] == curTime]

            startStr = str(self.numParticles) + "\nAir Particle Simulation\n"
            line = str.encode(startStr) 
            os.write(xyzFile, line)

            try:
                assert(rows.shape[0] == self.numParticles)
            except:
                curTime += self.dt
                continue
            # end

            # Get the postition, temperature and kinetic energy data
            # from these rows
            for rowIdx in range(self.numParticles):
                # Get Position Data
                xPos = rows.iloc[rowIdx]["Position-X [m]"]
                yPos = rows.iloc[rowIdx]["Position-Y [m]"]
                zPos = rows.iloc[rowIdx]["Position-Z [m]"]

                # Get Temperature
                temp = rows.iloc[rowIdx]["Temperature [K]"]

                # Get Kinetic Energy
                ke = rows.iloc[rowIdx]["Kinetic Energy [J]"]

                # Create string to write to xyz file
                dataStr = "H " + str(xPos) + " " + str(yPos) + " " + str(zPos) + "\n"
                line = str.encode(dataStr)
                os.write(xyzFile, line)
            # end for

            curTime += self.dt
        # end for
    # end writeXYZ
# end Simulation class