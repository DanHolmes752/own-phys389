# Project Documentation

# 02.03.2023 - 09.03.23

## Constant Mag Field Prototype

Time steps that look to work are:
- dt = 1e-10
- t_max = 1e-6

Using B = 0.05, v0 = 1e4

Gives a decently small error which can be attributed to using a basic integration method (Euler-Cromer). Need to calculate time steps per period:

$$
T  = \frac {2 \pi m_p} {qB} = 1.3136*10^{-7}  s
$$

From this we can calculate that the number of time steps for one period is 1313.67 . Not very useful as the code T might be different to theory T, so lets find code T:

### Finding code T and testing numerical vs analytical
-  T is time for one period. Can find this via total time / number of orbits:
    - Have a counter every time it passes a certain coordinate to find no. of orbits.
- Using $ \omega = 2\pi / T$ we can find freuency of gyration and check against theory value from: $ \omega = qB / m$
- Graphical comparison is a valid and perfectly good way to test and compare.
    - Plot x-position against time, should be sinusoidal, and compare against sin(theory T).
    - Can fit sin plot to position time plot.
    - Plot x-pos against time, normalise to get 1>y>1, do inverse sine. This gives 2 pi t / T so we can find a vlue of T.
- Test script could get plots for orbit, KE and momentum conservation.

Next: 
- Implement this to object orientated code.
- Push documentation to doc folder instead of project folder.
- Try Runge-Kutta in the notes, should give more accurate values for less time steps as it involves more calculations. Potentially more computationally intensive. 



# 13.03.2023

## Simulation.py

- Now contains cyclotrons attributes and related methods of generating a bunch of protons and running the simulation.
- generateProtons function will be used to randomly initialise a bunch of protons within the bounds of the cyclotron.

## tesingSim.py

- Currently where you run the simulation from using OOP.
- Contains two functions:
    - orbitTest, which runs the sim and plots x against y positions.
    - sinTest, which runs the sim and plots x position against time step. This should show a sinusoidal plot, closer to sin(gyration T) the better the energy conservation.



# 20.03.2023

## Proton.py & Electron.py
- Created these scripts to include the sub classes: Proton and Electron which inherit from Particle.
    - Allows more effective assignment of objects, ie dont have to input specfic particle properties when making object now, just its position velocity and acceleration.

## Simulation.py
- Updated to incorperate Proton.py and Electron.py

## initTest.py
- Created a script to test for the correct construction of objects.
    - Creates proton and electron objects and prints out their respective (and known) attributes.

## integrationTests.py
- Created a script to test each integration method. Basic unit test.



# 24.03.2023

## Simulation.py
- Added serparate functions of "runSim()" for each numerical method.
- Added saving position data to .npy files for each numerical method.
- Cleaned up some comments and placement of variables.

## testingSim.py
- Updated functions to read saved data as opposed to running the sim each time to return the data.
- Separated orbitTest() into three functions for each numerical method.
- Improved sinTest to plot all sinusoidal positions against each other and sin($\omega$ t_max)



# 26.03.2023

## Particle.py
- Added Runge-Kutta method.

## Simulation.py
- Updated to save KE values for each particle.
- Changed use of np.empty to np.zeros
    - Made it clearer which layer was the one that was initially generated to provide array structure

## testingSim.py
- Updated script to capture all functions into a class named gyrationTests.
- Updated gyration plot functions to work with n-body system.
- Methods now calculate their own t_list rather can creating a Cyclotron object to do this.
    - same t_list values, calculated the same as in Cyclotron class.

## angularFreq_test.py
- Made another script with another testing class.
    - This class contains the methods that test for energy and $\omega$ conservation.

## testHub.py
- Made a central script to run all tests from.
    - Just made organising testing a lot quicker and clearer.

## Other To Do:
- Git commits of this stuff^ (including name change for AF test file) and clean up project folder
- Make $\omega$ tests n-body
- Add E field
- If time (probs not) add class that represents particle bunches



# 26.03.2023

## cyclotronSimulation.py
- This is a copy of Simulation so that the E field could be added without ruining the constant mag field script.
- Added boundaries to B field. Cyclotron now has a radius that the B field operates in.
- Added oscillating E field in the boundary 2L where L = 0.05*cyclotron radius

## Particle.py
- added getRadius method to find the particle radius from the centre.
- Updated updateAcceleration and updateRungeKutta methods to account for added E field.

## particleBunch.py
- Involves the particleBunch class that includes the generateProtons and generateElectrons methods.
- Can be updated to represent the bunch of particles to find average position, velocity and KE.

## ALL TESTING FILES
- Updated to account for new cyclotronSimulation script and added E field.

## Testing notes
- These constants work well to show the final orbits:
    - dt = 1e-10
    - t_max = 1e-6
    - B = np.array([0, 0, 0.5])
    - numOfPro = 1
    - cyclotronRadius = 0.0015
    - E = 5e4
    - E_boundary = 0.05 * cyclotronRadius