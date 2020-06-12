---
title: 'June 5th, 2020'
---

# June 5th, 2020

Today I first looked at the paper that Anatoli had assigned me to look at (*Space-Based Cosmic-Ray and Gamma-Ray Detectors: a Review by Luca Baldini*). The paper is somewhat similar to a book, containing the following (in brief terms):
- Introduction to what cosmic rays are / history of them
- relavant variables associated with describing data obtained from cosmic ray (things like differential flux / intensity, energy, rigidity, grammage), as well as potential confusions associated with them
- how to interpret cosmic ray spectra
- the instrumentation used to detect cosmic rays / gamma-rays
- basic facts about atmosphere and the geomagnetic field of Earth that needs to be considered when detecting cosmic rays, and ways to accomodate for this in the model
- possible interactions of radiation (charged particles) with matter
  - energy loss
  - Coulomb scattering
  - pair-production
  - EM showers
  - hadrnoic showers
  - Cherenkov radiation
- Experimental techniques to detect cosmic rays

The facts I listed above are just the brief topics that are studied in this paper, I really need to read this carefully to understand whats going on.

## What do I need from this paper?
I think, just by skimming through this paper, that the first few pages concerning introduction to the relavent variables and interpretation of plots are definitely important. But for my project, this is not required just yet. I can always refer back to this paper if I need information for this.

What I need to read more about is **the things about the geomagnetic field and ray-tracing techniques that are implemented**. This directly relates to the tracing simulation project that I am currently doing. I should look into this more tomorrow. 

## The package

I also briefly skimmed through the package relavent to the paper. This contained a lot of programs that I wasnt sure about. A more thorough look-through is definitely required, as this would be my beginning steps to my project. 

## About my project

Here I just want to reiterate what my project is to allow me to collect my thoughts and summarize from the talk I had with Anatoli regarding this. 

### Short Description
Creating a GPU-based cosmic ray particle tracking simulation

### Long Description
Creating a simulation package that simulates the tracks that cosmic rays take upon reaching to Earth in spherical coordinates (a 3-D sense). This simulation should take into account the Earth's magnetic field, the charge of the particles, possible trajectories of cosmic rays etc. The tracking procedure should be somewhat similar to ray-tracing performed in optics. The package should be GPU-based to allow a very fast approach to running this simulation. 

### Steps
1. Set-up the code by using the trivial approach (B-field is constant and homogeneous, and we only simulate for one particle). This should be done in Python for now, and the package from Baldini should augment this.
2. Implement this for ~1e9 particles. This would probably result in some performance issues, and we may need to bring some parts of the code to an accelerator (C++ / FORTRAN).
3. Do this for the Earth's actual magnetic field. 
4. Add a different type of particle (muons first, then neutrinos)
5. Lastly, implement spin-polarization to the particles. 


### Next Steps
- Read through the paper about the near-earth environment part
- Look through and try to run the package from Baldini
- Think of a good name for the project. 
  - Possible options for now:
    - STraPCR : Simulation for Tracking Particles from Cosmic Rays
  - The project name should include:
    - tracking
    - cosmic rays
    - GPU