// create location map
#ifndef __LOCATIONMAP_H_
#define __LOCATIONMAP_H_

#include <map>
#include <utility>
#include <string>
#include "Location.h"

namespace locations {    
    // create locations
    Location kamioka = Location("Kamioka", 36.434800, 137.276599, 0.);
    Location icecube = Location("IceCube",-89.99, 0., 0.);
    Location uofa = Location("UofA", 53.523230, -113.526319, 0.);
};

// initialize map
extern std::map<std::string, Location> locmap = {
    {locations::kamioka.name(), locations::kamioka},
    {locations::icecube.name(), locations::icecube},
    {locations::uofa.name(), locations::uofa}
};

#endif // __LOCATIONMAP_H_