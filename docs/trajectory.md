# Implementation of Trajectory Algorithm

## Idea

- Creates a trajectory class that manages the particle trajectory of a given type, energy, location (longitude, latitude, altitude), and zenith / azimuthal angle
- has all of these as members
- there should be some checker where it checks if either...

  - particle has escaped (some escape condition)
  - particle is trapped (looping around Earth)
  - particle has reached Earth (r < Earth Radius)

- there should also be some way to compare between points, which is currently done by the TrajectoryPoint class

  - Should we make each ($\vec{r}, \vec{p}$) into a TrajectoryPoint for each iteration?

- There should also be some drawing functionality that draws the trajectory

  - This can also be done in a different class that uses the trajectory as an argument

### Current Issues

- [x] The zenith and azimuthal angle is not implemented properly. When using the same longitude, latitude, altitude, particle, and particle energy and only changing the zenith / azimuthal angle, the trajectory returned is not originating from the same initial location. 

  - can we fix this by making the trajectory class have all longitude / latitude / altitude / zenith / azimuth all as members and recording the initial location as some other object or something?
  - This now works, this is what we did:

    - We first created a simple test script that doesnt deal with all the member things so that we can focus on fixing this mistake.
    - We then translated the geodesic coordinates to geocentric coordinates:

      - $\begin{array} r = R_E + h \\ \theta = (90 - \lambda)(\frac{\pi}{180}) \\ \phi = \eta (\frac{\pi}{180})\end{array}$
      - here, $\theta\in[0, \pi]$ and $\phi\in[-\pi, \pi]$ with $\lambda\in[-90, 90]$ (with 90 being at the North Pole), and $\eta\in[-180, 180]$ (with 0 being at the Prime Meridian).

    - We then translated the coordinates (defined in the local tangent plane at some latitude and longitude) into the earth-centered, earth-fixed geocentric coordinates in ($r, \theta, \phi$) using a transformation matrix. This was obtained from [here](http://archive.psas.pdx.edu/CoordinateSystem/Latitude_to_LocalTangent.pdf). 

      - The equation used was here:

        - $\vec{r}_e = \vec{r_0} + \vec{R}_e^t \cdot \vec{r}_t$ where $\vec{R}_e^t = \left(\begin{array}{ccc} -\sin (\eta) & -\cos (\eta) \sin (\lambda) & \cos (\lambda) \cos (\eta) \\ \cos (\eta) & -\sin (\lambda) \sin (\eta) & \cos (\lambda) \sin (\eta) \\ 0 & \cos (\lambda) & \sin (\lambda) \end{array}\right)$ with $(\lambda, \eta)$ being the latitude (with $\lambda = 0$ defined at the geographic equator) and longitude (with $\eta=0$ defined at the Prime Meridian). Here, $\vec{r}_e$ is the coordinate vector defined in ECEF coordinates, $\vec{r}_0$ is the coordinate vector of the origin of the local tangent plane in terms of ECEF coordinates, and $\vec{r}_t$ is the coordinate vector defined in the local tangent plane.

      - The coordinate vectors in the local tangent coordinate plane are easily obtained from using the provided zenith and azimuthal angles, as using this, we can utilize spherical coordinate conversion equations to convert the unit vector of the velocity and position into Cartesian equivalents.

    - For the velocity transformation, we took the time derivative of the equation above, and since the latitude / longitude is constant throughout the particle propagation, we effectively reduced the equation as such: $\vec{v}_e = \vec{R}_e^t \cdot \vec{v}_t$. Here we want to simulate a particle that comes at some zenith / azimuth with the given particle energy.
    - Lastly we convert the transformed coordinate and velocity Cartesian components into spherical ones. This is easy for the coordinate, but a lot harder for velocity due to unit vector conversions as well. Nevertheless this was done relatively simply.

  - The results obtained is consistent with the 12 GV antiproton (fig 25) on Baldini's review paper
  - It is not consistent with the East-west effect on Baldini's paper, which is a bit weird. Hopefully this is due to my approximations (earth as dipole etc)

- Later, we will have to reorganize this into our original format with the trajectory classes, optionally changing some structuring up.



### Reorganizing the Trajectory class
- Previously we had 2 classes:
  - ParticleTrajectory class that creates some trajectory constructor at some given latitude, longitude, and altitude for some particle at some energy
    - This class creates the trajectory using the member getTrajectory with the zenith and azimuth as arguments
  - TrajectoryPoint class that records the latitude, longitude, and altitude at some point in the trajectory.

- What do we want to do now?
- Some ideas:
  - The main Trajectory class
    - required members:
      - particle: the label of the particle 
      - latitude, longitude: the geodesic coordinates of the detector
      - altitude: the altitude of the detector, i.e. the stopping altitude of the particle trajectory defined from sea level in the local tangent plane
        - these could alternatively be replaced with some location object, but im not so sure about the convenience / performance
      - zenith angle: the angle from the local zenith in which the particle comes from
      - azimuth angle: the angle with 0 being in the direction of the East hemisphere from the Prime Meridian in the local tangent plane coordinate system
      - rigidity / energy: the particle rigidity or energy when arriving at the detector
        - This should work in that if one is given, the other is evaluated vice versa
        - default energy = None, default rigidity = None
        - required since if we arent provided with the values, then we can use the fact that the parameter is None to then evaluate the other parameter
    - optional members:
      - escape altitude: the altitude in which the iteration will stop at, i.e. the particle has "escaped" earth's magnetic field 
        - default at 1000 km 
    - members without parameters:
      - energy / rigidity: as above
      - particleEscaped: boolean to decide if particle has effectively escaped or not
    - member function 1: get the transformation matrix
    - member function 2: get the initial location and the first location of the particle as TrajectoryPoints

    - member function 3: trajectory creator getTrajectory
      - optional members:
        - max steps: the maximum number of steps for the runge kutta integration (the overflow control parameter)
          - default at 10000 steps
        - step size: the step size for the RK integration
          - default 0.1 s
      - This should do the following:
        1. get the initial location
        2. transform the second location based on the zenith and azimuth angle
        3. append these as TrajectoryPoints
        4. get RK variables from RK integrator
        5. convert velocities into momenta
        6. either initialize TrajectoryPoint and append that or append the values themselves
        7. iterate until max_steps or until r < RE or r > (RE + escape altitude)
    - member function 4: get the x,y,z values of the trajectory from TrajectoryPoints
      - get x, y, z values from each TrajectoryPoint
      - append x, y, z into some class dictionary

  - Some TrajectoryPoint class
    - required members:
      - momentum ($p_r, p_\theta, p_\phi$) in spherical coordinates
        - This may not need to be implemented for now, we will see
      - latitude, longitude, altitude
      - all of these would have default values incase we want to initialize and give members values using other methods (except for particle)
    - member function 1: set from spherical coordinates
      - uses spherical coordinates ($r, \theta, \phi$) and transform them into latitude, longitude, and altitude values
    - member function 2: get r, theta, phi values from the members (a converter)
    - member function 3: get Cartesian coordinates from members (another converter)
    - member function 4: get Cartesian momentum from spherical momentum (yet another converter)
