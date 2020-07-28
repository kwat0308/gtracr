// header file for Magnetic Field
#ifndef __MAGNETICFIELD_H_
#define __MAGNETICFIELD_H_

#include <math.h>

#include "constants.h"

class MagneticField {
 private:
  // double* gcoeffs;
  // double* hcoeffs;
  double g10;  // mean value of the magnetic field at the magnetic
               //   equator

 public:
  // Constructor
  MagneticField() : g10{-2.94048 * (1e-5)} {}
  // Destructor
  // ~MagneticField() {delete[] gcoeffs; delete[] hcoeffs;}
  // MagneticField();
  // the radial component of the Earth's magnetic field
  const double Br(const double &r, const double &theta, const double &phi) {
    return 2. * (constants::RE / r) * (constants::RE / r) *
           (constants::RE / r) * g10 * cos(theta);
  }
  // the polar component of the Earth's magnetic field
  const double Btheta(const double &r, const double &theta, const double &phi) {
    return (constants::RE / r) * (constants::RE / r) * (constants::RE / r) *
           g10 * sin(theta);
  }
  // the azimuthal-component of the Earth's magnetic field
  const double Bphi(const double &r, const double &theta, const double &phi) {
    return 0.;
  }
};

#endif  // __MAGNETICFIELD_H_