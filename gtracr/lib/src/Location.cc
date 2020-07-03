// class that dictates locations around the globe
/*
A class of locations around the globe. Used to allow easy access to geodesic coordinates of locations of interest.
Members:
- name : the name of the location 
- latitude : the geographical latitude (0 = equator) in decimal degrees
- longitude : the geographical longitude (0 = prime meridian) in decimal degrees 
- altitude : the altitude above sea level of the location in km
*/
#include "Location.h"

// Constructor
// Default constructor at 0, 0
Location::Location()
    : lm{"Null Island"}, lat{0.}, lng{0.}, alt{0.}
{
}

// Constructor with default altitude
Location::Location(const char *name, const double &latitude, const double &longitude)
    : lm{name}, lat{latitude}, lng{longitude}, alt{0.}
{
}

// Constructor with variable altitude
Location::Location(const char *name, const double &latitude, const double &longitude, const double &altitude)
    : lm{name}, lat{latitude}, lng{longitude}, alt{altitude}
{
}

// copy constructor
Location::Location(const Location &loc)
    : lm{loc.lm}, lat{loc.lat}, lng{loc.lng}, alt{loc.alt}
{
}

// Copy assignment operator
Location &Location::operator=(const Location &loc)
{
    delete[] lm;
    lm = loc.lm;
    lat = loc.lat;
    lng = loc.lng;
    alt = loc.alt;
}
