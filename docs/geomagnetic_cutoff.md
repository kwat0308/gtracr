# Geomagnetic Cutoff

## What is it?

- a quantitative measure to determine which particle trajectories are allowed or forbidden

  - usually done by observing the rigidity of the particle
  - rigidity : the momentum per charge of the particle

## Algorithm

- [x] Chooses a single type of particle

  - for our sake, this would be a proton
  - for now, this can and will be changed in the near future
  - We should implement this so that the particle type can be changed easily (like some user input)

- [x] Choose some location

  - [x] Kamioka : where the Super-Kamiokande is

    - 36.434800, 137.276599 (in decimal degree notation)
    - 36°26'05.3"N 137°16'35.8"E (in decimal minute second (DMS) notation)

  - [x] South Pole : where IceCube is

    - 89.99 -63.453056 (decimal notation)
    - 89° 59' 24", -63° 27' 11" (in DMS)

  - we could implement this for more locations later on

  - To make this easier, we could store this data in some dictionary

- [x] Choose some energy for the particle

  - Test the waters
  - ~0.5GeV (low energy limit), ~1GeV, ~10GeV, ~100GeV should be good enough for now.
  - We should make a functionality where we can add more energies with ease when higher energies are supported

- [x] get the trajectories for each zenith and azimuthal angle

  - we want to start from high zenith angles to lower ones (higher zenith angle should be more allowed than lower zenith angle)

    - for azimuthal it doesnt really matter all that much

  - we need to have a way to convert the provided trajectory in longitude, latitude, altitude, zenith, and azimuth angle to spherical

  - currently we only have this for longitude, latitude, and altitude so we need this for zenith and azimuth as well

    - To do this, we can retain our current version or add on to it

  - Probably have to implement a double for-loop for this (for each zenith angle, for each azimuth angle)

- [x] In the end, we want to have a checking condition whether the particle can be a cosmic ray or not

  - this can be done by seeing if the particle touched the Earth sometime within its trajectory
  - Also seen when the particle has looped at least twice around the Earth

    - if there are 2 occurrances of the same longitude or latitude a total of 2 times

  - Other conditions may be needed, not sure what they would be

    - Either way we can integrate that on a later date, focus on those above for now.

  - we also want to return the rigidity of each **initial** trajectory, which will help us determnine the cutoff value for the rigidity as well

    - this will provide the quantitative measure for the geomagnetic cutoff

- [x] With the result (True of False), we can plot a 2-d heatmap

  - zenith and azimuth angles on each axis, and yes / no as the color axis
  - the title should specify the location and altitude and the particle type and energy
  - store it as a .png file for now, we can think of a better storage method later on

## Implementation

- create a new file that does all of the above

  - alternatively we can make this into a member function of the ParticleTrajectory class
  - But just to keep things separate for now, I think its better to make it its own file

    - we can always integrate it as a member function later on if we want to

- small things like conversions will be in the utility class like always

- we should also have a new file that can create a dictionary of locations based on longitude and latitude -

### Implementing the zenith and azimuthal angles

After drawing some geometrical figures, I think I have the equation to properly convert a 3-D vector ($s, \xi, \alpha$) defined on the local reference frame (sphere) to the geocentric coordinates ($r, \theta, \phi$) (given initial geocentric coordinates ($r_0, \theta_0, \phi_0$).)

Some definitions of the variables:

- $l$ is the altitude along the zenith from the tangent plane
- $\xi$ is the zenith angle defined locally
- $\alpha$ is the azimuthal angle defined locally
- ($r_0, \theta_0, \phi_0$) defines the location of the detector _after_ conversion from longitude and latitude coordinates to geocentric ones

  - $r_0$ here represents the Earth's radius

- ($r, \theta, \phi$) are our usual spherical coordinate system that will be used for the integration process

With these definitions, we have our conversion equation as follows:

- $(r \cos\phi)^2 = (r_0 \cos \phi_0)^2 + \left( \dfrac{l\cos\alpha}{\cos z}\right)^2 - 2r_0l \cos \alpha\cos\phi_0$

  - derived from Law of Cosines with the projection of the vectors onto the $(r, \theta)$-plane.

- $\theta = \theta_0 - \left(\dfrac{-l\cos \alpha \tan z}{r \cos \phi}\right)$

  - derived from Law of Sines and simpler triangle properties

- $\phi = \phi_o - \arctan \left(\dfrac{l\tan z \cos\alpha}{r_0\tan\theta_0}\right)$

  - derived from projection of 3-vector in local frame onto tangent plane with projection of tangent plane onto the $(r, \phi)$-plane
  - we observed that the tanget plane is at a right angle to the initial vector, and as such transforming the projected vector on the tanget plane to the Cartesian axes was not too hard.

So only $\phi$ can be evaluated completely using the variables provided. The other two variables $r, \theta$ must be evaluated from $\phi$. So we want to evaluate the coordinates like this:

1. Get ($r_0, \theta_0, \phi_0$) from our already existing conversion code for latitude and longitude to geocentric coordinates
2. Get $\phi$ from above
3. Evaluate for $r \cos\phi$ from the first equation
4. Use the results from 3\. to determine $\theta$
5. Divide the results from 3\. by $\cos \phi$ to get $r$.

Hopefully this allows a perfect conversion...

### Current issues

- [x] There should be a way to check each trajectory point to get a boolean yes or no for allowed / forbidden trajectories - Seems like this works now

  - My idea is to either:

    - Remove zenith and azimuthal angles as members of the TrajectoryPoint class and make some converter from horizontal coordinates to geodesic ones with zenith and azimuth as arguments

      - most promising one right now
      - This might work, the conversion is relatively simple:

        - The projection of the 3-D vector $\vec{s}$ at the initially given latitude and longitude $(L_o, L_a)$ is given as $s\sin \xi$ with $\xi$ being the zenith angle (angle from the zenith)
        - Using this, we can get the new latitude and longitude by utilizing the azimuthal angle $\alpha$, which would be $s\sin\xi\cos\alpha$ and $s\sin\xi\sin\alpha$ respectively.
        - Since we have the altitude $l$ instead of the 3-D vector, we replace $s$ with $\dfrac{l}{\cos\xi}$ instead (as the altitude is the projection of the 3-D vector onto the zenith).

  - remove those guys as members and create a new class that contains the trajectorypoint along with zenith and azimuthal

  - Leave them as members and create the converter within the class

    - least configuration necessary

- [x] The code works, however there are still some issues as the heatmaps all result in True (i.e. no forbidden trajectories at all). This can be caused by a few things:

  - Incorrect conversions between coordinates, although I feel as though this is already resolved pretty well.

    - This is still not resolved, there seems to be some overflow issue with the angular coordinates. Essentially we need to be able to get the latitude and longitude defined within [-90, 90] and [-180, 180] respectively. In order to do this, we want to divide by some integer so that we get the correct latitude / longitude.

      - perhaps we should divide them by 90 and 180 respectively?

  - Inappropriate termination conditions / checking conditions for cutoffs

    - right now the only factor is if the particle touches the earth again, but I am not even sure if this is implemented properly.

- Here the $\phi$ is now defined from [$-\pi, \pi$] instead of [$0, 2\pi$]. This allows easier conversions for longitude with $\phi$.

- There is still a problem with the cutoffs, seems like at any energy the particle doesnt reach the Earth at all, only for electrons at very low energies (~0.005GeV), which is definitely weird. I inputted the values from Baldini but they still do not give the correct trajectories. This is certainly weird and needs to be investigated later on.

### Now

- After fixing the Trajectory class, the problem with the weird cutoff behaviour was fixed. The geomagnetic cutoff is now completely integrated with the new Trajectory class.
- After doing a small test trial, it seems like this works as I intended it to, as the output for one plot has some cutoffs based on the zenith and azimuth angle.
- Now we need to perform some cleanup to make the outputs have more convenience for the user and add some utility functions.

### Initial Results

Here we present the initial results of the cutoffs. We performed this at 3 locations, namely:

- Kamioka (at latitude: 36.434800, longitude: 137.276599)
- IceCube / South Pole (at latitude -89.99, longitude: 0 (longitude doesnt matter))
- University of Alberta (latitude: 53.523230, longitude: -113.526319)

We did this with 4 different particles at 4 different rigidities ($R = \dfrac{pc}{Ze}$):

- Particles: positron, electron, proton, anti-proton
- Rigidities: 5, 10, 30, 50 (in GV)

### After showing results to Anatoli

- They werent visually appealing (black and white is a no go)
- We should reduce the number of plots required
- The colorbar was not used at all, so instead of presenting cutoff plots for each energy, we can do this all together in one plot with different shades + transparancy
- Optimization of the code is necessary, reduce the time of the actual process.

  - our goal is ~4 milliseconds, but 4 seconds is good for now.
  - This is now achieved, with ~30 seconds per iteration. Of course, this can be more optimized (look at the optimization markdown file for more details).

### Creating contour plots

The presentation of these rigidity cutoffs should be more elegant. An example would be from Honda's paper in 2002 (<https://arxiv.org/pdf/hep-ph/0203272.pdf>) where they use contour plots to indicate at which locations are the cutoffs for the rigidities (i.e. any rigidity **below** that is allowed).

So we should create some sort of plot that is similar to this, by utilizing color or whatever. Not too sure how to do this, but we should think about this for now.
