// namespace for particle map in which more particles can be inserted
#ifndef __PARTICLEMAP_H_
#define __PARTICLEMAP_H_

#include <map>
#include <utility>
#include <string>
#include "Particle.h"

namespace particles {    
    // create particles
    Particle ep = Particle("positron", -11, 0.5109*(1e-3), 1, "e+");
    Particle em = Particle("electron", 11, 0.5109*(1e-3), -1, "e-");
    Particle pp = Particle("Proton", 2122, 0.937272, 1, "p+");
    Particle pm = Particle("anti-proton", -2122, 0.937272, -1, "p-");
};

// initialize map
extern std::map<std::string, Particle> partmap = {
    {particles::ep.label(), particles::ep},
    {particles::em.label(), particles::em},
    {particles::pp.label(), particles::pp},
    {particles::pm.label(), particles::pm},
};

#endif //__PARTICLEMAP_H_