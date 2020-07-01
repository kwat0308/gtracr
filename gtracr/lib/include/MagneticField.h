// header file for Magnetic Field
#ifndef __MAGNETICFIELD_H_
#define __MAGNETICFIELD_H_

class MagneticField {
    private:
        double* gcoeffs;
        double* hcoeffs;

    public:
        // Constructor
        MagneticField();
        // Destructor
        ~MagneticField() {delete[] gcoeffs; delete[] hcoeffs;}
        // MagneticField();
        const double Br(const double&, const double&, const double&);

        const double Btheta(const double&, const double&, const double&);

        const double Bphi(const double&, const double&, const double&);
};

#endif // __MAGNETICFIELD_H_