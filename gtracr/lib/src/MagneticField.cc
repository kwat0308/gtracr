// magnetic field class

#include <math.h>
#include "constants.h"
#include "MagneticField.h"

MagneticField::MagneticField()
    // : gcoeffs{new double[13 * 13]}, hcoeffs{new double[13 * 13]} 
{
    // gcoeffs[0] = -29404.8 * (1e-9);
    g10 = -29404.8 * (1e-9);
}
// components of magnetic fields in spherical coordinates
// r-component
const double MagneticField::Br(const double &r, const double &theta, const double &phi)
{
    return 2. * (constants::RE / r) * (constants::RE / r) * (constants::RE / r) * g10 * cos(theta);
}
// theta-component
const double MagneticField::Btheta(const double &r, const double &theta, const double &phi)
{
    return (constants::RE / r) * (constants::RE / r) * (constants::RE / r) * g10 * sin(theta);
}
// phi-component
const double MagneticField::Bphi(const double &r, const double &theta, const double &phi)
{
    return 0.;
}
