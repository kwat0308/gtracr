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
- [ ] The zenith and azimuthal angle is not implemented properly. When using the same longitude, latitude, altitude, particle, and particle energy and only changing the zenith / azimuthal angle, the trajectory returned is not originating from the same initial     location. 
  - can we fix this by making the trajectory class have all longitude / latitude / altitude / zenith / azimuth all as members and recording the initial location as some other object or something?