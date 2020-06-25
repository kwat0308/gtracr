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
- [x] The zenith and azimuthal angle is not implemented properly. When using the same longitude, latitude, altitude, particle, and particle energy and only changing the zenith / azimuthal angle, the trajectory returned is not originating from the same initial     location. 
  - can we fix this by making the trajectory class have all longitude / latitude / altitude / zenith / azimuth all as members and recording the initial location as some other object or something?
  - This now works, this is what we did:
    - We first created a simple test script that doesnt deal with all the member things so that we can focus on fixing this mistake.
    - We then translated the geodesic coordinates to geocentric coordinates:
      - $\begin{array}r = R_E + h \\ \theta = (90 - \lambda)(\frac{\pi}{180}) \\ \phi = \eta (\frac{\pi}{180})\end{array}$
      - here, $\theta\in[0, \pi]$ and $\phi\in[-\pi, \pi]$ with $\lambda\in[-90, 90]$ (with 90 being at the North Pole), and $\eta\in[-180, 180]$ (with 0 being at the Prime Meridian). 
    - We then translated the coordinates (defined in the local tangent plane at some latitude and longitude) into the earth-centered, earth-fixed geocentric coordinates in ($r, \theta, \phi$) using a transformation matrix. This was obtained from [here](http://archive.psas.pdx.edu/CoordinateSystem/Latitude_to_LocalTangent.pdf). 
      - The equation used was here:
        - $\vec{r}_e = \vec{r_0} + \vec{R}_e^t \cdot \vec{r}_t$ where $\vec{R}_e^t = \left(\begin{array}{ccc}
-\sin (\eta) & -\cos (\eta) \sin (\lambda) & \cos (\lambda) \cos (\eta) \\
\cos (\eta) & -\sin (\lambda) \sin (\eta) & \cos (\lambda) \sin (\eta) \\
0 & \cos (\lambda) & \sin (\lambda)
\end{array}\right)$ with $(\lambda, \eta)$ being the latitude (with $\lambda = 0$ defined at the geographic equator) and longitude (with $\eta=0$ defined at the Prime Meridian). Here, $\vec{r}_e$ is the coordinate vector defined in ECEF coordinates, $\vec{r}_0$ is the coordinate vector of the origin of the local tangent plane in terms of ECEF coordinates, and $\vec{r}_t$ is the coordinate vector defined in the local tangent plane. 
      - The coordinate vectors in the local tangent coordinate plane are easily obtained from using the provided zenith and azimuthal angles, as using this, we can utilize spherical coordinate conversion equations to convert the unit vector of the velocity and position into Cartesian equivalents. 
    - For the velocity transformation, we took the time derivative of the equation above, and since the latitude / longitude is constant throughout the particle propagation, we effectively reduced the equation as such: $\vec{v}_e =  \vec{R}_e^t \cdot \vec{v}_t$. Here we want to simulate a particle that comes at some zenith / azimuth with the given particle energy. 
    - Lastly we convert the transformed coordinate and velocity Cartesian components into spherical ones. This is easy for the coordinate, but a lot harder for velocity due to unit vector conversions as well. Nevertheless this was done relatively simply. 
  - The results obtained is consistent with the 12 GV antiproton (fig 25) on Baldini's review paper
  - It is not consistent with the East-west effect on Baldini's paper, which is a bit weird. Hopefully this is due to my approximations (earth as dipole etc)
- Later, we will have to reorganize this into our original format with the trajectory classes, optionally changing some structuring up.
