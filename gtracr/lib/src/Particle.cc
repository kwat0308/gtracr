/*
Utility class for cosmic ray particles
Members:
- name: the name of the particle (char*)
- pid: the particle id as in pdg (int)
- mass: the particle rest mass (double) [units of GeV / c^2]
- charge: particle's charge Z (int) [units of elementary charge]
- label: the shorthand name for the particle (char*)

The below members require the particle energy or one of the below members as additional information
- momentum: Particle momentum
- velocity: particle velocity
- rigidity: particle rigidity

*Default is set to Proton

Notes:
- PDGID obtained from here: http://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf
- The mass of the particles are also obtained from PDG

Example:
proton: proton = Particle("Proton", 2212, 0.938272, "p+")
*/

#include <math.h>
#include "constants.h"
#include "Particle.h"

// Constructors
// Default constructor, set to proton
// Default energy set to 1GeV
Particle::Particle()
    : nm{"Proton"}, pid{2212}, m{0.938272}, ch{1}, lbl{"p+"}
{
    p = sqrt(m * m + 1.);
    v = (p * constants::sc) /
        sqrt(
            (p * p) + (m * constants::sc) * (m * constants::sc));
    R = p / abs(ch);
}

// Construct with given initial configurations
// default energy set to 1GeV
Particle::Particle(const char *name, const int pdgid, const double &mass, const int charge, const char *label)
    : nm{name}, pid{pdgid}, m{mass}, ch{charge}, lbl{label}
{
    p = sqrt(m * m + 1.);
    v = (p * constants::sc) /
        sqrt(
            (p * p) + (m * constants::sc) * (m * constants::sc));
    R = p / abs(ch);
}

// Constructor with some provided energy
Particle::Particle(const char *name, const int pdgid, const double &mass, const int charge, const char *label, const double &energy)
    : nm{name}, pid{pdgid}, m{mass}, ch{charge}, lbl{label}
{
    p = sqrt(m * m + energy * energy);
    v = (p * constants::sc) /
        sqrt(
            (p * p) + (m * constants::sc) * (m * constants::sc));
    R = p / abs(ch);
}

// Destructor
Particle::~Particle()
{
    delete[] nm;
    delete[] lbl;
}

// copy constructor
Particle::Particle(const Particle &part)
    : nm{part.nm}, pid{part.pid}, m{part.m}, ch{part.ch}, lbl{part.lbl}
{
    p = part.p;
    v = part.v;
    R = part.R;
}

// copy assignment operator
Particle &Particle::operator=(const Particle &part)
{
    // const char* nam = part.nm;
    // const char* lb = part.lbl;
    // delete member pointers
    delete[] nm;
    delete[] lbl;
    nm = part.nm;
    pid = part.pid;
    m = part.m;
    ch = part.ch;
    lbl = part.lbl;

    p = part.p;
    v = part.v;
    R = part.R;
}

// setters
// void set_momentum_energy(const double& energy)
// {
//     p = sqrt(mass()*mass() + energy*energy);
// }
