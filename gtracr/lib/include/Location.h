// header file for Location class
#ifndef __LOCATION_H_
#define __LOCATION_H_

class Location
{
private:
    std::string ln;    // Location name
    double lat; // Latitude
    double lng; // Longitude
    double alt; // Altitude
public:
    // Constructors
    Location();  // Default
    Location(const std::string &, const double &, const double &); // altitude set to default
    Location(const std::string &, const double &, const double &, const double &);
    // destructors
    ~Location() {}
    // Copy constructor / operator
    Location(const Location&);
    Location &operator=(const Location&);
    // getters
    const std::string& name() {return ln;}
    const double latitude() {return lat;}
    const double longitude() {return lng;}
    const double altitude() {return alt;}

};

#endif //__LOCATION_H_