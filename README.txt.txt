This file will explain how to run the code, with brief outlines of each python script.

HOW TO RUN THE SIMULATION:

- testHub.py serves as the central script to run the code and any tests from. From here, you can provide the parameters for the simulation and 
  which case (cyclotorn or const magnetic field) to run. 
	
	- To just run the final simulation, call the "gyrationTests.calculate_EM_Data" method. This will run the simulation from "cyclotronSimulation.py" 
	  using all three numerical methods, saving the positional and KE data to .npy files.

	- The remainder of the listed functions are used to produce plots of the particle properties, i.e. positional, angular frequency and KE plots.

	- To run the simplified case of the constant magnetic field, you can either use the method above while inputting the constant "E" to equal 0 
	  OR you can use the "gyrationTests.calculateData" method, which runs the simplified simulation from "Simulation.py", again using all numerical
	  methods and saving data.



FILE DESCRIPTIONS:

- testHub.py : A central script where you run the simulation and any tests / plot functions from.
- cyclotronSimulation.py : The main file for the cyclotron simulation, contains the EM_ cyclotron class which contains the methods can calculate and save positional and KE data.
- Simulation.py : The file which was used to test the simplified case of the constant magnetic field. Nearly identical layout to above, just without E field code.
- Particle.py : Contains the Particle class used to initialise particle objects.
- Proton.py & Electron.py : Contains the subclasses Proton and Electron which inherit from parent class Particle. Used to initial proton and electron objects.
- particleBunch : Contains the particleBunch class which contains the methods that generate the single or bunch of particle used in the simulation.
- testingSim.py : Contains the gyrationTests class which contains the method used to test and plot the gyrations of the particle.
- energyConservation.py : Contains the testingConservation.py class which contains the methods which test and plot the angular frequencies and KE of the particle.
- prototype.py : The initial procedural code made to capture the scope of the project.

