import os
from sim import Simulation
from pandas import DataFrame
from timeit import default_timer as timer
#####################################################################
# Read inputs in from file
start = timer()
inputVals = []
with open("./sim.in", "r") as input:    
    # Read each line in the file
    for line in input:
        # Get the line
        stripped = line.strip()
        
        # Get what comes after the equal sign and strip
        val = stripped.split("=")
        val = val[1].strip()
        
        # Add to list of input values
        inputVals.append(val)
    # end for
# end with

# Assign to correct value
numParticles = int(inputVals[0])
size = float(inputVals[1])          # [m]
timeStep = float(inputVals[2])      # [s]
endTime = float(inputVals[3])       # [s]
highTemp = float(inputVals[4])      # [K]
lowTemp =  float(inputVals[5])      # [K]
tempPct = {highTemp: float(inputVals[6]),
           lowTemp: float(inputVals[7])}
outputFile = inputVals[8]
xyzOutputFile = inputVals[9]

print("Simulation Parameters:")
print("     Num Particles:      %d" % numParticles)
print("     Simulation Size:    %f [m]" % size)
print("     Time Step:          %f [s]" % timeStep)
print("     End Time            %f [s]" % endTime)
print("     High Temperature:   %f [K]" % highTemp)
print("     Low Temperature:    %f [K]" % lowTemp)
print("     Percentages:        High - %f, Low - %f" % (tempPct[highTemp], tempPct[lowTemp]))
print("     Output File:        %s" % outputFile)
print("     XYZ File:           %s\n" % xyzOutputFile)

# Run Section
print("Initializing Simulation...")
sim = Simulation(numParticles, size, timeStep, endTime, tempPct)
print("Simulation Initialized\n")

print("Initialize Dataframe/Table...")
# Create dataframe to store
dataCols = {"Time [s]": [],
            "PID": [],
            "Position-X [m]": [],
            "Position-Y [m]": [],
            "Position-Z [m]": [],
            "Velocity-X [m]": [],
            "Velocity-Y [m]": [],
            "Velocity-Z [m]": [],
            "Temperature [K]": [],
            "Kinetic Energy [J]": [],
            "Collision": [],
            "Collided Particles": []}
simData = DataFrame(dataCols)

# Update dataframe with first frame data
simData = sim.updateSimTable(simData)
print("Dataframe and Table Initialized\n")

print("Run Simulation Loop")
# Begin simulation loop
simData = sim.runSimulation(simData)
print("Simulation Loop Complete")
print("Simulation Data Saved at: %s\n" % outputFile)

print("Write data to XYZ file format")
# Open output file
file = open(outputFile, "w")

# Open XYZ file
open(xyzOutputFile, "w")
xyzFile = os.open(xyzOutputFile, os.O_WRONLY)

# Write sim data to csv file
simData.to_csv(file, index=False)

# Write sim data to .xyz file format
sim.writeXYZ(simData, xyzFile)
print("XYZ File Created at: %s\n" % xyzOutputFile)

# Close files
file.close()
os.close(xyzFile)
end = timer()
print("Simulation Run Time: %f" % (end - start))