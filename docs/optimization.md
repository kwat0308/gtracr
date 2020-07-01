# Optimization

## Profiling
I ran a couple of profiling scripts onto the script **geomagnetic_cutoff.py**. There are some functions that require some optimization:
- the initializers, this takes like ~0.2 seconds for each initialization
- making the layout into a tight_layout (0.339 s)
- drawing the function onto a plot (~1.2 seconds)
- Saving the figure (~6.5 seconds)
- Magnetic field computation (takes ~20 seconds for the whole thing, ran ~3.6 million times)
- Runge_kutta computation (~246 seconds, ran 453 thousand times)
- the DEs (ran about 1.8 million times each, takes ~70 seconds)
- the weighted sum function (takes 6 seconds, ran ~2.7 million times)
- vmag_spherical (~42 seconds, ran ~5.4 million times)
- gamma (~41 seconds, ran ~5.4 million times)
- Trajectory initialization (takes about 12.8 seconds for 400 times)
- TrajectoryPoint initialization (~2.5 seconds for 453 thousand times)

Seems like we need to configure the runge kutta part of it, as it is ran the most times and takes the most time.

Other notes:
- **geomagnetic_cutoff.py** takes about 293 seconds
- **evalTrajectory** in Trajectory class takes total of 260 seconds, which is the majority of the time of the function

Some fixes that I have implemented:
- Effectively removed the gamma and vmag_spherical functions and wrote them explicitly in locations where they were used
  - This reduced the time by half!! Reducing Pythonic function calls really did help :thumbsup:
- Created classes for runge_kutta and magnetic field, and initialized them when the getTrajectory member function is called.
  - This didnt change anything... (maybe a few milliseconds??)

## Transitioning to C++
To further optimize the code, I brought the RungeKutta part (and subsequently the MagneticField part) to C++. This will greatly increase performance as the time taken from function calls in runtime will be effectively removed as this will all be done in compile time in C++.

### Current progress
- I was able to make this work within C++. I created a test function within c++ and ran a simple testor with this. The variables clearly change with each iteration. 
  - However, the values still do not update when brought into python. Further investigation is required.