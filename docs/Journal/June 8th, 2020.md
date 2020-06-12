# June 8th, 2020

## Goals for today

- [x] Create the package directory and put into GitHub
- [x] Implement / Think of how to implement the beginning steps of the project
- [ ] Get some ideas for the current steps of the project

## Package title name
After talking with Anatoli, I have a concrete idea on what my project is about. The short description is something like this:
> A GPU-based cosmic ray tracking simulation: Tracks cosmic rays at any location in Earth (in a 3-D sense) using real geomagnetic field data.

From this, I thought of the package name ***STraCRGPU***. This would be an acronym to ***Simulation for Tracking Cosmic Rays using / accelerated by GPUs***. Hopefully this is a good enough name...

Another name I came up with is ***GSimTraCR***, which is an acronym for ***A GPU-based Simulation Tracking package for Cosmic Rays***. This sounds better maybe?

So Ive decided to use the name **GSimTraCR**, sounds pretty nice and it does contain all the necessary contents for my project. Hopefully Anatoli approves of it :thumbsup: 

## Researched material for project
Anatoli gave me some things to research on as well:
> do some research if there are particle tracking libraries for GPU available, or what kind of literature exists on this topic. Also commercial /branded applications, like from nVidia or Intel. 

So it seems like I should look up some particle tracking libraries that utilizes GPUs / literature available for particle tracking. 

### Some literature that I found for particle tracking using GPUs:
- **Compass**: particle tracking algorithm that is optimized for GPUs
    - Used by people in LHCb
    - Uses CUDA
    - [ResearchGate Article](https://www.researchgate.net/publication/334435818_A_Parallel-Computing_Algorithm_for_High-Energy_Physics_Particle_Tracking_and_Decoding_Using_GPU_Architectures) (Full version available here)
- **SixTrackLib** : Purely Pythonic 6D single charged particle tracking library
    - Link to presentation of the program from [here](https://indico.cern.ch/event/833895/contributions/3577803/attachments/1927226/3190638/intro_sixtracklib_handouts.pdf).
    - Actually not really useful if we want to implment neutrinos (neutrinos have no charge) 
- **GEANT4** : The primary tool CERN uses to simulate particle collisions in accelerators. Uses OpenGL as its graphics library. Could we use / take away from this somehow as well?
- Symplectic multi-particle tracking model : some model that uses CUDA, not so sure about the package name
        - link: [INSPIRE-HEP](https://inspirehep.net/literature/1627309), [PDF](http://accelconf.web.cern.ch/ipac2017/papers/thpab027.pdf)

## Implementation (initial ideas)

### What should I include?
- a particle library that has identities for different particle types
    - mainly PID, charge, mass, (any others?)
- an RK-integrator that performs Runge-Kutta integration to track a single particle's propagation
    - This is the part where we may need to do this on some accelerator (C++, FORTRAN)
    - We use the Lorentz force equation for this $\vec{F} = \dfrac{d\vec{p}}{dt} = q(\vec{v} \times \vec{B} + \vec{E})$
        - from Martin Walt : Introduction to geomagnetically trapped radiation (on file explorer)
- A trajectory tracer that outputs some values that has the particle properties / locations at each discrete step
    - this class should control the actual RK integrator, much like **gTrajectoryTracer** in Baldini's script.
- a class that controls the earth magnetic field which utilizes / is utilized by **pyIGRF** (more on this below).
    - we might not really need this, just an idea incase we do need it
- a plotting script that creates some canvas for plotting 
    - This should create some figure / canvas that shows the 3-D plot of the trajectories in real time
    - Maybe some more scripts that augment the main plotting script?
    - either way the plotting scripts should be very well-versed and involved.
- Some CPU-GPU parallelization script
    - This would be the hardest part, I wouldnt have much clue as to how to do this
        - maybe **CuPy** instead of **NumPy**?
    - Probably dont need this as of yet, only required when we want to implement with more particles

### Some random ideas
- We could use [**pyIGRF**](https://pypi.org/project/pyIGRF/)  (International Geographic Reference Field) to model the magnetic field around Earth. This would help with particle tracking with Earth's magnetic field. Baldini's code also uses the FORTRAN version of this in his FORTRAN file for particle tracking (why FORTRAN :cry:).
    - Alternatively, I could assume an ideal dipole field and implement it myself. This would make it easier for me to add configurations but will result in less accurate results. 
- I think that the E-W effect actually comes into play just by the interconnection between the Earth's magnetic field with the charged particle trajectories. So maybe we dont actually have to physically implement this??
- We can use **VisPy**, a matplotlib-like module that creates 2-D / 3-D plots that are accelerated by GPUs. It seems like it has some hardware requirements, maybe we should stray off from this?
    - It requires OpenGL >= 2.1
    - [homepage](http://vispy.org/), [documentation](http://vispy.org/documentation.html), [GitHub](https://github.com/vispy/vispy)
    - But do we really need to **plot** using GPU?? Probably not, unless the output itself is very graphics heavy (there should be some way to augment this, so probably dont use this).
    - However this may be useful in the future when we have to implement ~1e6 particles at the same time. Then that will be graphically intense.
- I found the article in which Baldini refers to for numerical ray-tracing. I should have gotten this earlier...
    - link [here](https://www.sciencedirect.com/science/article/pii/S0273117705001997?via%3Dihub)
    - This article talks about how to implement the geomagnetic rigidity cutoff.



### Some things I learnt
- GPUs are primarily used due to its multitude of CPU cores (~1000) equivalent to 1 GPU core. So its good for parallalization programs.
    - However, GPUs dont have faster speed than CPUs, so its not well-versed in computation that requires high-speed performance (programs that cant be resolved with parallalization).


