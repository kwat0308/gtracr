# June 9th, 2020

## Goals for today
- [x] Solidify initial ideas for implementation
- [x] Get some abstract yet concrete steps for implementation
- [x] Start implementing the easier things (and figure out what to implement first)
    - This would be the integrator (most important and crux of the package)

## Implementation 

### More implementation ideas
- For our RK integrator, we need to use the Lorenz Force equation $F = \dfrac{d\vec{p}}{dt} = q(\vec{v} \times \vec{B} + \vec{E})$ as shown before
    - ~~we need to make this a differential equation somehow.~~ Look below.

- For the earth magnetic field class, maybe we want to use an ideal dipole field first so that it fits with the Lorenz equation first.
    - This would actually be better as the cutoff rigidity is much easier to implement in the dipole approximation (from the article by D.F.Smart). 
    - The geomagnetic cutoff rigidity is determined by the Störmer's equation: $R_{c}=\left[M \cos ^{4} \lambda\right]/\left\{r^{2}\left[1+\left(1-\sin \epsilon \sin \xi \cos ^{3} \lambda\right)^{1 / 2}\right]^{2}\right\}$where:
        - $R_c$ = geomagnetic cutoff rigidity [MV]
        - $M$ = magnetic dipole moment [Gcm^3]
        - $\lambda$ = latitude from magnetic equator
        - $\epsilon$ = angle from zenith direction (zenith direction = radial from position of dipole center)
        - $\xi$ = azimuthal angle measured clockwise from the direction to North magnetic pole
        - $r$ = distance from dipole center [cm]
    - they are all in CGS units so we need to convert them to SI units
    - we can obtain the magnetic dipole moment from IGRF database
    - If we only want the vertical cutoff rigidity, then we have an analytic solution for the ideal dipole field: $R_{v}=\frac{R_{0} \sin ^{4} \theta}{r^{2}}=\frac{R_{0} \cos ^{4} \lambda}{r^{2}}=\frac{R_{0}}{L^{2}}$ where $R_0 \approx 14.5GV$. 
        - the vertical cutoff rigidity is different from the geomagnetic cutoff rigidity in that it normalizes Störmer's equation to distance in Earth radii (units of $R_E$ from the dipole position) and obtains the cutoff value relative to the local zenith.
        - Not too sure if we only want this, this would be useful only if we are doing the simulation at a single location. If we want to map the whole Earth then this wouldnt be enough.

### Obtaining the differential equation
#### Assumptions:
- We are in a source-free region (outside the inside of the Earth) so that $\nabla \times \vec{B} = \vec{0}$
- no E-field present (E=0)
- non-relativistic particles (this would have to be changed later for sure)
- Earth's magnetic field is an ideal dipole field (first-order approximation of B-field)


#### Steps
- We first express the vector into its spherical coordinate equivalants (use physics definition of sph. coord.):
    - $\vec{p} = p_r\hat{r} + p_\theta \hat{\theta} + p_\phi \hat{\phi}$
    - $\vec{v} = v_r\hat{r} + v_\theta \hat{\theta} + v_\phi \hat{\phi}$
    - $\vec{B} = B_r\hat{r} + B_\theta \hat{\theta} + B_\phi \hat{\phi}$

- Now if we are working in an ideal dipole field, then we already have the analytical expression for the Earth's magnetic field:
    - $B_r = 2\left(\dfrac{R_E}{r}\right)^3 g_1^0 \cos\theta$
    - $B_\theta = \left(\dfrac{R_E}{r}\right)^3 g_1^0 \sin\theta$
    - $B_\phi = 0$
- Now just evaluate the cross product between velocity and magnetic field directly. Then we get:
    - $(\vec{v} \times \vec{B})_r = -B_\theta v_\phi$
    - $(\vec{v} \times \vec{B})_\theta = B_r v_\phi$
    - $(\vec{v} \times \vec{B})_\phi = B_\theta v_r - B_r v_\theta$
- Lastly equate this to each component in Lorenz Force equation, since in **non-relativistic regime** we have that $p = mv$. With 3 additional differential equations for the position, we then get 6 coupled differential equation as follows:
    - $\dfrac{dv_r}{dt} = \dfrac{-B_\theta v_\phi}{m}$
    - $\dfrac{dv_\theta}{dt} = \dfrac{B_rv_\phi}{m}$
    - $\dfrac{dv_\phi}{dt} = \dfrac{B_\theta v_r - B_r v_\theta}{m}$
    - $\dfrac{dr}{dt} = v_r$
    - $\dfrac{d\theta}{dt} = v_\theta$
    - $\dfrac{d\phi}{dt} = v_\phi$
- Finally from this, we can successfully evaluate the position $(r, \theta, \phi)$ and momentum $(p_r, p_\theta, p_\phi)$ of a **non-relativistic** particle in an **ideal dipole field**. 

#### Additional considerations:
- If we want to simulate particle trapping (drift), it would probably be necessary to use the guiding center approach to trace the particle trajectory. The method for this is on Walt's book at pages 11~22. 
- In these pages, there are also methods into adding on the E-field term, so if this would be necessary in the future this would be helpful.
- To do this for relativisitic particles, we use the relativistic equation $E^2 = p^2c^2 + m^2c^4$ where $m$ is the intrinsic mass of the particle, and $E$ is the total energy of the particle. 
    - Currently not too sure what / how to make this part into a differential equation after this, more looking into is necessary. 
    - If we are dealing with high energy particles (> 50GeV) then we could probably make the approximation $E \gg mc^2$ so that $pc \approx E$.
- Later on we might want to apply adaptive step-sizing for better performance
    - Check numerical recipes or ask Patrick

### Steps for implementation
- Create the RK integrator and see how that works first
    - Do this for some simple system first, then do this for our non-relativistic, ideal dipole field system
        - do this for an electron, then for some other particle for another test case
    - Check the output to the analytic solution for the simple case to see if our implementation is correct for the RK integration
- Then plot this on some 3-D thing (PlotLy? or simple mpl thing)
- We then create a proper plotter with (relatively) nice formatting
- Then we can create the trajectory tracer class and connect the plotter to the RK integrator
    - Here in this step we should also consider creating the magnetic field class and deem if its even necessary or not.
- Create the particle library that creates particle dictionaries.


## Progress for today
- Implemented a simple Particle class
    - has mass, pid, charge, name as members
    - we can evaluate things like energy or rigidity from this
    - A question for this: should I be using Planck system of units (i.e. c=1)? and why?
- Implemented initial version of RK-integrator 
    - I am not sure how to make it more compact since each RK variable depend on each other. If there was a way to do this I would do it. Keep looking i guess?
    - I should also think of a way to deal with a nicer way for inputs
- Tested this RK integrator out with equations above
    - Seems to output some trajectory-like thing, however without right units accuracy check is essentially useless.
    - I think this would eventually be put into the TrajectoryTracker class
- Changed package name
    - its now called ***gTraCR***, an acronym for ***A GPU-based Tracking simulation for Cosmic Rays***. Sounds pretty similar to before, however this has a nicer ring to it. 

## Next Steps
- Clean up the RK integrator code
- Create the TrajectoryTracker class 
    - This class would do the following:
        - know where to start (get initial conditions) based on zenith angle and stuff 
        - also know where to end (due to some cutoff / limiting value)
        - Run the integrator for *some* point and store that in an array (this means i need to fix RK so that it takes a single value instead of the max iteration , i.e. no iteration should occur in the RK integrator). But thats kinda weird as integration = loops, so I should really think of something else
        - It should also then just return the coordinate data for the plotter to plot
        - This would also use the Particle class so that we can use charge / mass info from different particles here instead of the actual integrator
- Create another script that creates a dictionary of particles