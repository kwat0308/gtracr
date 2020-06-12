---
title: 'June 6th, 2020'
created: '2020-06-05T20:45:16.982Z'
modified: '2020-06-07T23:45:42.526Z'
---

# June 6th, 2020

## Goals for today:
- [x] Look through the near-field environment section of Baldini's paper
- [x] Look through code / package from Baldini and try to run stuff
- [x] Think of a good name for my package and start creating project directory for this

## The geomagnetic field
- Outside Earth we have a source-free region so $\vec{J} = 0$. This implies that $\nabla \times \vec{B} = 0$, which allows us to describe the magnetic field outside Earth's surface as the negative gradient of a scalar potential: $\vec{B} = -\nabla \psi$. 
- The scalar potential can be described using spherical harmonics: 
  $\psi(r, \theta, \phi)= R_{E} \sum_{n=1}^{\infty}\left(\frac{R_{E}}{r}\right)^{n+1} \sum_{m-0}^{n}\left(g_{n}^{m} \cos m \phi+h_{n}^{m} \sin m \phi\right) P_{n}^{m}(\cos \theta)$
- We can always assume azimuthal symmetry.

### Ideal Dipole Field

- approximation where Earth is a perfect dipole: we truncate only to the lowest (dipole) term (n=1, m=0). 
- $g_1^0$: The field intensity at the equator on the Earth's surface
- We usually use the McIlwain parameter $L = \dfrac{r}{R_E \sin^2\theta}$ to describe the distance in which the field line crosses the equator
  - B-field is not constant on the same L-shells: $B(\theta) = \dfrac{B_0}{L^3} \dfrac{\sqrt{1+3\cos^2\theta}}{\sin^6\theta}$.
- Together with $B(\theta)$ and $L$, we call them the McIlwain coordinates, used to describe the motion of slow, charged particles in Earth's B-field

### Actual geomagnetic field

The actual geomagnetic field loses a lot of assumptions made in the ideal case:
- dipole is not perfectly aligned to the coordinate axis, its ~11 degrees off
- the center of the Earth is not the same as center of the dipole
- There is some asymmetry in the currents within Earth, so higher-order terms are necessary

Regardless, the McIlwain parameter is still useful to describe the Earth's B-field.

### Ray-tracing techniques
The easiest way to do this is to numerically solve the classical equation of motion for the trajectories. What we should do is that:
- Dont choose particles that arrive from some arbitrary source (inefficient)
- Instead generate particles that have the opposite charge of the detector and propagate them to infinity
- Utilize the Lorentz force $F = q\vec{v} \times \vec{B}$.

### Criteria / things to watch out for:
#### Geomagnetic cutoff:
- the cutoff in which particles can enter Earth's atmosphere
- The cutoff provides a separation between allowed and forbidden trajectories.
- Done by using the **vertical rigidity cutoff**, the minimum rigidity required for a charged particle to reach the point above Earth's surface from some local zenith at some altitude. 
  - This is evaluated by simulating particle trajectories of the opposite charge at the provided altitude and longitude / latitude. 
- Particles that have a rigidity below the vertical rigidity cutoff are shielded by the B-field from Earth. 
- For an ideal dipole field, this can be anatically determined to $R_v = \dfrac{14.5\text{GV}}{L^2}$. 
- Rigidity: ratio between momentum and charge of particle: $R = \dfrac{pc}{Ze}$.

#### East-West effect
- asymmetry in flux of particles due to its charge and their interaction with the B-field
- particles from one direction are (mostly) forbidden due to its behavior with the interaction with Earth's B-field
  - This allows the asymmetry in flux
- higher rigidity makes the trajectories more straight, thus diminishing the effect

#### Secondary radiation in low-Earth orbit
- particles can be trapped by interactions of the primary cosmic ray with the atmosphere
- these particles can bounce back and forth (i.e. drift) around the geomagnetic field lines between 2 mirror points.
- Trapped particles have a much higher positron ratio than from primary cosmic rays

#### Separating charges using geomagnetic field
- This essentially utilizes the vertical rigidity cutoff as well as the East-West effect. 
- We know that for positive charges, there are much more allowed trajectories from the West direction than from the East. The same can be applied for negative charges. Using this, we can simply distinguish the two particles
- We also need their energies in order to distinguish between them (?)
  - not so sure, look this up more.

## Next Steps
- Goals weren't met properly, I should focus a bit more on that. 
- I should definitely look into the code and try to run it and produce the output that Baldini got in his paper
- Then think of ways to integrate this into the project.




