// magnetic field class

#include <math.h>
#include "constants.h"
#include "MagneticField.h"

MagneticField::MagneticField()
    : gcoeffs{new double[13 * 13]}, hcoeffs{new double[13 * 13]} 
{
    gcoeffs[0] = -29404.8 * (1e-9);
}

const double MagneticField::Br(const double &r, const double &theta, const double &phi)
{
    return 2. * (constants::RE / r) * (constants::RE / r) * (constants::RE / r) * gcoeffs[0] * cos(theta);
}

const double MagneticField::Btheta(const double &r, const double &theta, const double &phi)
{
    return (constants::RE / r) * (constants::RE / r) * (constants::RE / r) * gcoeffs[0] * sin(theta);
}

const double MagneticField::Bphi(const double &r, const double &theta, const double &phi)
{
    return 0.;
}
