// header file for Location class
#ifndef __LOCATION_H_
#define __LOCATION_H_

class Location
{
private:
    std::string ln; // Location name
    double lat;     // Latitude
    double lng;     // Longitude
    double alt;     // Altitude
public:
    // Constructors
    Location();                                                                    // Default
    Location(const double &, const double &, const double &);                      // no location name
    Location(const std::string &, const double &, const double &, const double &); // with location name
    // destructors
    ~Location() {}
    // Copy constructor / operator
    Location(const Location &);
    Location &operator=(const Location &);
    // getters
    const std::string &name() { return ln; }
    const double &latitude() const { return lat; }
    const double &longitude() const { return lng; }
    const double &altitude() const { return alt; }
    // utility functions
    void print();
};

#endif //__LOCATION_H_