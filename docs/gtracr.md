# gtracr

## Goals for the project

- Create a simulation package that tracks cosmic rays at any location around the Earth at any zenith and azimuthal angle

  - This should be able to detect between allowed and forbidden trajectories
  - Also should be able to detect the location on Earth in which trajectories from the Earth came from (i.e. neutrino sources)
  - This should also be able to simulate the cosmic ray for other particles (electron, positron, antiproton)

- Accelerate the package using GPU cores

## Implementation

### Libraries we want

- Particle library

  - constructs a Particle object with mass, charge, ID (from PDG database)
  - The particles used are contained in a dictionary in which we can refer to within other parts of the code

- ParticleTrajectory library

  - Controls the trajectory of the particle based on its charge, mass, its location (longitude, latitude, altitude), and its path angle (zenith and azimuthal angle)
  - An intermediary between the actual integration and the Particle library
  - Various cutoffs / control parameters will be set here

- utility library

  - Has all utility functions and constants required
  - Things like conversion factors, coordinate transformations are contained here

- An integrator library

  - This should only do the integration
  - Currently a 4th-order Runge-Kutta algorithm is implemented

- A magnetic field library

  - Controls the function for Earth's magnetic field
  - uses the IGRF coefficients for the real magnetic field
  - Currently only dipole approximation is applied

### Evaluators

- A geomagnetic cutoff evaluator

  - evaluates the valid zenith / azimuth angles for a certain particle at a certain energy / rigidity
  - returns a heatmap of valid ranges as a function of zenith and azimuthal angle

### Tests

- some test function that can accurately test our code
- Test using a proton at 12 GeV from altitude of 500km at (0,0) in longitude and latitude
- Outputs a 3-D trajectory and a 2-D trajectory projected on (x,y) plane

### Relavent formulas

- Differential equation used for equation of motion of charged particles:
- $\dfrac{\mathrm{d} p_{r}}{\mathrm{d} t}=\dfrac{1}{m\gamma}\left(e\left(p_{\theta} B_{\phi}-p_{\phi} B_{\theta}\right)+\dfrac{p_{\theta}^{2}}{r}+\dfrac{p_{\phi}^{2}}{r}\right)$
- $\dfrac{\mathrm{d} p_{\theta}}{\mathrm{d} t}=\dfrac{1}{m\gamma}\left(e\left(p_{\phi} B_{r}-p_{r} B_{\phi}\right)-\dfrac{p_{r} p_{\theta}}{r}+\dfrac{p_{\phi}^{2}}{r \tan \theta}\right)$
- $\dfrac{\mathrm{d} p_{\phi}}{\mathrm{d} t}=\dfrac{1}{m\gamma}\left(e\left(p_{r} B_{\theta}-p_{\theta} B_{r}\right)-\dfrac{p_{r} p_{\phi}}{r}-\dfrac{p_{\theta} p_{\phi}}{r \tan \theta}\right)$
- $\dfrac{dr}{dt} = \dfrac{1}{m\gamma}\left(p_r\right)$
- $\dfrac{d\theta}{dt}=\dfrac{1}{m\gamma}\left(\dfrac{p_\theta}{r}\right)$
- $\dfrac{d\phi}{dt} = \dfrac{1}{m\gamma}\left(\dfrac{p_\phi}{r\sin\theta}\right)$
- 
- Vertical rigidity cutoff: $R_V = \dfrac{R_0\sin^4\theta}{r^2}$

  - This is a broad approximation and only applicable with the dipole approximation
  - We should instead evaluate cutoffs by testing the validity of the trajectory at every zenith and azimuth.

## Current Issues

- [ ] Code is not structured as nicely as I want it to be, more code organization can be performed

  - Using more classes to regroup certain aspects of the code
  - Understanding what we need / what we dont need as a function / class

- [x] Zenith and azimuthal angle is not implemented in the code

  - More about this below.

- [ ] We are using / borrowing some ideas from Luca Baldini's scripts; licensing will be an issue here

  - He uses the GPL (General Public License) from GNU
  - But it seems like he is really for free software usage (Source: [his "about me" page](http://osiris.df.unipi.it/~baldini/aboutme.html)), so asking him about the licensing issue may be a play here.

- [x] Not sure if particles are performed in a relativistic or non-relativistic sense.

  - RESOLVED: The mass component was supposed to be relativistic, so we added the Lorentz factor to account for this

- [x] There are some confusions between the starting and stopping altitude:

  - [x] Whats an appropriate altitude limit? (Currently starting altitude set to be 565km, stopping to be 2km)

    - RESOLVED: we set the starting altitude to be 1km and stopping altitude to be 500km. 1km seems to be an appropriate limit, but we should be able to do this with 0km (at sea level). the 500km is still arbitrary, but seems to work out fine (for now)

  - [x] We want to start from the detector location, not stop there. For this to work then, we should make our start altitude at the location of the detector, and the stop altitude to be something else.

    - RESOLVED: this is now properly implemented.

  - [x] ~~Would we even need a stop altitude?~~ Yes, this is required as a limit to where we should stop the trajectory integration process.

- [ ] Only dipole approximation currently applied

- [x] The DEs seem to work, however I am not sure why they work. This is important to check for my reference.

## Past Issues

- I accidentally implemented the momentum within my DE instead of velocity. To fix this, I let the conversions between velocity and momentum occur in the Trajectory class and Particle class, and apply relativistic assumptions as well. This now seems to work fine (i.e. looks like a particle trajectory :thumbsup:)

# Tasks

## Major Tasks

- [ ] Figure out things about the licensing with Baldini's scripts
- [x] Obtain a heatmap of allowed and forbidden trajectories at a single location and altitude and at all zenith and azimuthal angles (i.e. determine the geomagnetic cutoff for each location on Earth)

## Minor Tasks

- [x] Check out how the DE works
- [ ] Clean up the code structure
- [x] Get the zenith angle / azimuthal angle coordinate transformation working
- [x] Figure out a system to detect allowed vs forbidden trajectories
- [x] Get the geomagnetic cutoff code working

### Fixing some conventional issues

Naming conventions. These really need to be consistent. We should be following the PEP8 standard for Python, and the Google C++ standard for C++.

For naming conventions, the following will be used:

- Class names: CamelCase
- global variables: UPPER_CASE
- everything else (functions, variables, methods / member functions): snake_case

Additionally, I decided to abide to the following conventions for variables concerning classes in C++:

- members : variable_ (i.e. add underscore **after** the variable name)
- parameters for constructors: variable (the actual variable name)
- getters : variable() (i.e. the actual variable name)
- parameters in setters: _variable (i.e. add underscore **before** variable name)

### Restructuring the code

Anatoli also gave me a pep talk about standards when creating a package. The key thing is this: **_Each class should have one single function_** i.e., they should not be juggling many tasks all at once.

He gave me an analogy like a factory, where different sections work independently but have a secure communication line. The brain / control center communicates with certain sections when the brain wants certain things, and those sections in return communicate with lower sections that create the raw materials. But the lower sections **do not** communicate with the brain (generally).

In a similar fashion, code structure should obey the same sort of style.
