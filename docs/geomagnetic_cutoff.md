# Geomagnetic Cutoff

## What is it?

- a quantitative measure to determine which particle trajectories are allowed or forbidden

  - usually done by observing the rigidity of the particle
  - rigidity : the momentum per charge of the particle

## Algorithm

- [ ] Chooses a single type of particle

  - for our sake, this would be a proton
  - for now, this can and will be changed in the near future
  - We should implement this so that the particle type can be changed easily (like some user input)

- [ ] Choose some location

  - [ ] Kamioka : where the Super-Kamiokande is

    - 36.434800, 137.276599 (in decimal degree notation)
    - 36°26'05.3"N 137°16'35.8"E (in decimal minute second (DMS) notation)

  - [ ] South Pole : where IceCube is

    - 89.99 -63.453056 (decimal notation)
    - 89° 59' 24", -63° 27' 11" (in DMS)

  - we could implement this for more locations later on
  - To make this easier, we could store this data in some dictionary

- [ ] Choose some energy for the particle

  - Test the waters
  - ~0.5GeV (low energy limit), ~1GeV, ~10GeV, ~100GeV should be good enough for now.
  - We should make a functionality where we can add more energies with ease when higher energies are supported

- [ ] get the trajectories for each zenith and azimuthal angle

  - we want to start from high zenith angles to lower ones (higher zenith angle should be more allowed than lower zenith angle) 

    - for azimuthal it doesnt really matter all that much

  - we need to have a way to convert the provided trajectory in longitude, latitude, altitude, zenith, and azimuth angle to spherical
  - currently we only have this for longitude, latitude, and altitude so we need this for zenith and azimuth as well

    - To do this, we can retain our current version or add on to it

  - Probably have to implement a double for-loop for this (for each zenith angle, for each azimuth angle)

- [ ] In the end, we want to have a checking condition whether the particle can be a cosmic ray or not

  - this can be done by seeing if the particle touched the Earth sometime within its trajectory
  - Also seen when the particle has looped at least twice around the Earth

    - if there are 2 occurrances of the same longitude or latitude a total of 2 times

  - Other conditions may be needed, not sure what they would be

    - Either way we can integrate that on a later date, focus on those above for now.

  - we also want to return the rigidity of each **initial** trajectory, which will help us determnine the cutoff value for the rigidity as well

    - this will provide the quantitative measure for the geomagnetic cutoff

- [ ] With the result (True of False), we can plot a 2-d heatmap

  - zenith and azimuth angles on each axis, and yes / no as the color axis
  - the title should specify the location and altitude and the particle type and energy
  - store it as a .png file for now, we can think of a better storage method later on

## Implementation

- create a new file that does all of the above

  - alternatively we can make this into a member function of the ParticleTrajectory class
  - But just to keep things separate for now, I think its better to make it its own file

    - we can always integrate it as a member function later on if we want to

- small things like conversions will be in the utility class like always
- we should also have a new file that can create a dictionary of locations based on longitude and latitude
-
