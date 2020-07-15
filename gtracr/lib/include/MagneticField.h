// header file for Magnetic Field
#ifndef __MAGNETICFIELD_H_
#define __MAGNETICFIELD_H_

class MagneticField {
    private:
        // double* gcoeffs;
        // double* hcoeffs;
        double g10;  // mean value of the magnetic field at the magnetic equator

    public:
        // Constructor
        MagneticField();
        // Destructor
        // ~MagneticField() {delete[] gcoeffs; delete[] hcoeffs;}
        // MagneticField();
        // the radial component of the Earth's magnetic field
        const double Br(const double &r, const double &theta, const double &phi);
        // the polar component of the Earth's magnetic field
        const double Btheta(const double &r, const double &theta, const double &phi);
        // the azimuthal-component of the Earth's magnetic field
        const double Bphi(const double &r, const double &theta, const double &phi);
};

#endif // __MAGNETICFIELD_H_