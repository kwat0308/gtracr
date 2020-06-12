---
title: 'June 7th, 2020'
created: '2020-06-07T20:32:00.149Z'
modified: '2020-06-08T15:44:14.533Z'
---

# June 7th, 2020

## Looking into Baldini's package

### Some small things I have noticed
- His package primarily uses **ROOT** in his modules, so I will download it and test his code from this.
  - However, I am not sure if I would want to use **ROOT** for my own package, this could be considered. But as of now I am not sure at all.
- To track particle trajectories, he uses Runge-Kutta integration (RK45). This would probably be doable. 
- In his article he mentions that "*Programs like IGRF numerically solves the classical equation of motion*", and he also uses this in one of his files from a FORTRAN program of IGRF. I looked up this and found a Python 3.7 equivalent from [here](https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html). We can definitely use this.
- Besides that, his project looks very well-implemented and solid, could we directly use his package instead of making our own from scratch? Or should we augment our package with the ideas that he has implemented?
- I still have no idea how to implement the vertical rigidity cutoff and the East-West effect... There's probably something about this in his code, I'll need to look into this more.
  - By looking into the package more, it seems like they have some sort of criterion for the trajectory of the particle that sets the cutoff value
  - As for the E-W effect, I am not too sure as it seems like they just drew 2 arbitrary trajectories that seem to recreate the E_W effect. 

### Main takeaways from his code
- The particle trajectories are controlled by the **TrajectoryTracer** class in **gTrajectoryTracer.py**. This seems to be the crux of the ray-tracing / particle trajectory tracing in his figures.
  - I would either want to directly use this class or take elements from it and make it my own. The latter may be preferred due to more control, but if we can just directly use his code (or if I am allowed to do so) then that would be good too. I should ask Anatoli about this.
- He also has a Python file named **gEarthSphere.py** that allows us to draw the Earth like in his paper. We should take some elements from this so that plotting would be easier.  

Key programs:
- **gTrajectoryTracer.py**: contains particle trajectory tracer
- **gParticle.py**: Creates a dictionary containing different particles (with its basic identities)
- **gEarthSphere.py**: contains creation of Earth map
- **fig_raytracing.py, fig_vertcutoff.py, fig_east_west.py**: Python files that create Figures 25, 26, and 28 in his paper. These are the relavent figures about the ray-tracing that we want.

The optional programs:
Of course, these optional programs are widely used within the programs above, but they could be replaced with something else (ex matplotlib instead of gCanvas). This would also help us avoid using ROOT. 
- **gCanvas.py**: contains creation of canvas in which plots are created in
- **tji2010x.for**: Seems like this contains implementation for RK4? This is a FORTRAN file, we should try to "rewrite" things from here to Python. 

### Things I should do
- Definitely try to recreate whatever he got and change things to understand his code better.
  - By doing so, I will be able to understand what to add / not add into my code and how to use this as a reference for my package
  - ~~For now, I am waiting for ROOT to download so that the python files can be run.~~ Now ROOT is properly installed.
  - I could also add some comments to help me understand the code better
  - I now ran the programs above and changed a few parameters to see whats going on.
- Honestly, I think we should think of an approach without ROOT, ROOT contains a lot of dependencies and we should make the package more universal than anything. Of course, if we have to, we should use ROOT, however, I think that we should try to find a way not to use it for now. 
- I also need to find a way to integrate GPU with this. Probably something like **CuPy** would work, however this needs more looking into. 

### After running the program...
I figured out a few things after I ran the three Python files above, as well with more questions...
- An issue in **gEarthSphere.py**: **ROOT.gPad** doesn't seem to have member **annotate**, which cause an error. However this was only used to plot the equitorial axes, and thus by commenting it out the problem was resolved.
- The particles change their behaviour based on the geomagnetic field when giving them different values of (latitude, longitude, altitude, zenith angle, energy). The particle type also heavily impacts the behaviour.
  - But just by changing these parameters, their behaviour seem to be impacted by the geomagnetic field somehow. I am still not sure how / where this is implemented. Should look into this. 
    - I think that it may be from the IGRF module, but I am honestly not sure yet.
  - I am still not sure how the vertical rigidity cutoff really works here (E-W effect comes directly from the particle behaviour with B-field). More looking into is necessary here.

## Questions and Next Steps
### Questions for myself
- How do the particles interact with the geomagnetic field in the program?
- How is the geomagnetic field implemented here?
- How does the vertical rigidity cutoff work here?
- How much from his project should I really take?
- Hopefully I dont have to read through the FORTRAN code...

### Questions for Anatoli
- Should we be using ROOT for my own project? : NO! He really doesnt like it...
- Could we use Baldini's code directly, or should I takeaway elements from his project instead? : No, I should do this from scratch and only use it as some reference
- Not really confident if I can make something this good of caliber...
- Clarify about the project so that I can make a good and appropriate package name : GPU-based tracking simulation for cosmic rays.

### Next Steps
- [x] Ask Anatoli about ROOT and if we can use his code directly
- [ ] If not, start to create the project
  - list off some ideas as to how to recreate his code just in Python
  - only take away things that I really need from his code
- [ ] Think of a good package name!!

