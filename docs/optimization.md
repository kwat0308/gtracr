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

- After initial progress, the trajectory obtaining process for the geomagnetic cutoff took ~30 seconds (50 seconds less than before)!

  - This is still slow, however. The reason for this is caused by the fact that:

    - We initialize a new object (TrajectoryPoint) object at runtime, which takes ~6 seconds in total.
    - The actual trajectory obtaining process is performed within Python, which means that the loop for each iteration of the runge kutta is done within Python. We all know how slow Python loops are...

  - Possible remedies to do this:

    - Move the Trajectory class (and actually all classes) to C++. This will effectively make everything done in compile time, which will make everything a lot faster!
    - Append the variables themselves (i.e. do not create these TrajectoryPoint objects in each iteration). This will remove organization, which is greatly wanted, at the cost of performance...
    - Somehow restructure the code so that the whole iteration loop is contained within the RungeKutta class

      - An idea of this would be to use spherical coordinates instead of TrajectoryPoints so that it would not be necessary to access things from Python in C++

        - Of course, we could make the C++ Runge Kutta class access things from Python, but this would not be necessary, as we would only want geodesic information at the start and end of the trajectory.

- I was able to make this work within C++. I created a test function within c++ and ran a simple testor with this. The variables clearly change with each iteration.

  - ~~However, the values still do not update when brought into python. Further investigation is required.~~ This now works, it was an issue with the stepsize being a constant private member. This was changed and now the loop iterates as normal.

- Now the porting of the RungeKutta portion of the code works well with our current code structure.

# Second round of optimization

After the geomagnetic rigidity cutoffs work in a well-fashioned manner (with the dipole approximation), we decided to further optimize the C++ components of the code to increase performance. In order to do this, we profiled the C++ portion by creating a testing script that iterates over 250 times of constructing the same trajectory tracer instance, and evaluating it.

## The profiling tools used

- ~~callgrind (via valgrind) with kcachegrind as the visualization tool~~

  - this cannot be used without using a remote GUI that connects the Ubuntu desktop in the threadripper with my Windows laptop. Right now I am using VNC client.
  - the tools are very handy with lots of command options and nice methods for visualization (albeit not as nice as flamegraphs though).
  - links: [callgrind](https://www.valgrind.org/docs/manual/cl-manual.html), [kcachegrind](https://kcachegrind.github.io/html/Home.html)

callgrind with kcachegrind only shows the CPU profile when running the program, as as such is not useful for time profiling (i.e. optimization for performance).

So now we use a new tool called uftrace (for now) that is a function profiler that works without much overhead. The visual interface doesnt show much, however. This can be done better with either using WSL instead of the Linux VNC or implementing the flame graphs.

**The alternative**: We found a nice profiler called _orbit_ that is light weight and performs dynamic profiling. The usage of this in Linux is somewhat confusing, however, so I need to figure this out somehow. For Windows this seems like it may work, as we can launch it as a `.exe` file.

## Concerning optimization flags in setuptools

setuptools creates a nice framework to create a package that compiles C/C++ extensions in a well-ordered manner. However, when we compile this normally, there are a lot of compilation arguments set by default:

```
g++ -Wno-unused-result -Wsign-compare -DNDEBUG -fwrapv -O2 -Wall -Wstrict-prototypes -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -pipe -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -pipe -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /home/keito38/anaconda3/include -DNDEBUG -D_FORTIFY_SOURCE=2 -O2 -isystem /home/keito38/anaconda3/include -fPIC -I/home/keito38/anaconda3/include -Igtracr/lib/include -I/home/keito38/anaconda3/include/python3.7m -c gtracr/lib/src/Particle.cc -o build/temp.linux-x86_64-3.7/gtracr/lib/src/Particle.o
```

^ an example of one compilation of C++ class using g++ with setuptools

From [this stackexchange website](https://stackoverflow.com/questions/6928110/how-may-i-override-the-compiler-gcc-flags-that-setup-py-uses-by-default), they tell that the flags that come at the end completely override those preceding it. So a natural solution to our problem is to either:

- disable optimization flags given in `-O2` and override it with the lowest optimization flag `-O0`
- disable all default flags by adding `% OPT-''` in the beginning of the command `python setup.py ...`

  - this I am not sure how it works yet ([this stackexchange website](https://stackoverflow.com/questions/19779416/remove-all-default-compiler-arguments-in-setup-py) seems to have the same problem, and a possible remedy is suggested. this needs to be tried later on).

So it seems like if we apply any higher optimization (i.e. `-O3`) as a `extra-compiler-arg` to `setup.py`, then the higher optimization will override it.

A surefire way that will work (although very ugly) is through "undoing" all the flags presented within the compiler by providing "undoing" flags in `extra-comipiler-args`. This is tedious and ugly, but it would work.

### Another thing to note...

When working on the threadripper through remote SSH, the compilation arguments seem to be reduced. This is particularly strange as to why, maybe because the threadripper is purely Linux, whereas using WSL we dont have this purely on linux.

For reference, this is the compiler arguments presented in the threadripper:

```
g++ -Wno-unused-result -Wsign-compare -DNDEBUG -g -fwrapv -O2 -Wall -g -fstack-protector-strong -Wformat -Werror=format-security -g -fwrapv -O2 -Wdate-time -D_FORTIFY_SOURCE=2 -fPIC -I/home/keito/.local/lib/python3.8/site-packages/pybind11/include -Igtracr/lib/include -I/usr/include/python3.8 -c gtracr/lib/src/Particle.cc -o build/temp.linux-x86_64-3.8/gtracr/lib/src/Particle.o
```

### What do these compiler arguments mean???

Too many weird syntax to understand with these C++ compilers. Looking through each one on the web, the definitions are presented below:

- -Wno-unused-result: dont warn user if functions doesnt use return value
- -Wsign-compare: warn when incorrect results are obtained in conversion between sign and unsigned values
- -DNDEBUG : disable assertions
- -g : enable the debugger
- -fwrapv : do not optimize code segments that may result in overflow errors
- -Wall: enable all warnings
- -fstack-protector-strong: enables stack protection for vulnerable functions (i.e. inserts a guard variable to stack frame, that is checked if it is overwritten or not. Return warnings if it has been overwritten)
- -Wformat : warn when types with formatting are incorrect
- -Werror=format-security : change warnings to errors when formatting is incorrect
- -Wdate-time : warn when macros relating to time are encountered
- -D_FORTIFY_SOURCE=2 : level 2 checks to string and memory manipulation functions
- -fPIC : allow machine code to be used at any address
- -march=nocona -mtune=haswell : specifications based on the cpu being used (intel in this case)
- -ftree-vectorize : perform auto-vectorization
- -fno-plt : do not allow PLT to resolve external functions (in which address are not known) in runtime (i.e. make them all done in link time)
- -pipe: makes compilation process faster - use pipes instead of temporary files

Most of these warnings are obtained from [here](https://gcc.gnu.org/onlinedocs/gcc/Warning-Options.html)
