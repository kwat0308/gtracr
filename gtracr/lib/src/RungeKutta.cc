// Runge Kutta integrator class

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
// #include <pybind11/stl.h>
#include <vector>
#include <math.h>
// #include <iostream>
#include "constants.h"
#include "MagneticField.h"
#include "RungeKutta.h"

// Constructor
// Requires the charge and mass of the particle
RungeKutta::RungeKutta(const int charge, const double &mass)
    : bfield{new MagneticField()}, coeff{charge / mass}, h{0.01}
{
    // for (int i=0; i<8; ++i){
    //     rk1[i] = 0.0;
    //     rk2[i] = 0.0;
    //     rk3[i] = 0.0;
    //     rk4[i] = 0.0;
    // }
}

// Constructor
// Requires the charge and mass of the particle and stepsize
RungeKutta::RungeKutta(const int charge, const double &mass, const double &stepSize)
    : bfield{new MagneticField()}, coeff{charge / mass}, h{stepSize}
{
    // for (int i=0; i<8; ++i){
    //     rk1[i] = 0.0;
    //     rk2[i] = 0.0;
    //     rk3[i] = 0.0;
    //     rk4[i] = 0.0;
    // }
}

// Destructor
// RungeKutta::~RungeKutta()
// {
//     delete[] rk1;
//     delete[] rk2;
//     delete[] rk3;
//     delete[] rk4;
// }

// ODEs, inputs are given as numpy arrays of the form (time, r, theta, phi, vr, vtheta, vphi)

// dvrdt
double RungeKutta::dvrdt(double r, double theta, double phi, double vr, double vtheta, double vphi)
{
    double dvrdt = ((coeff / gamma(vr, vtheta, vphi, r, theta)) * (vtheta * bfield->Bphi(r, theta, phi) -
                                                                   bfield->Btheta(r, theta, phi) * vphi)) +
                   (r * vtheta * vtheta) + (r * vphi * vphi * sin(theta) * sin(theta));
    return dvrdt;
}

//dvthetadt
double RungeKutta::dvthetadt(double r, double theta, double phi, double vr, double vtheta, double vphi)
{
    double dvthetadt = (((coeff / gamma(vr, vtheta, vphi, r, theta)) * (vphi * bfield->Br(r, theta, phi) -
                                                                        bfield->Bphi(r, theta, phi) * vr)) /
                        r) +
                       (vphi * vphi * sin(theta) * cos(theta)) - (2. * vr * vtheta / r);
    return dvthetadt;
}

//dvphidt
double RungeKutta::dvphidt(double r, double theta, double phi, double vr, double vtheta, double vphi)
{
    double dvphidt = (((coeff / gamma(vr, vtheta, vphi, r, theta)) * (vr * bfield->Btheta(r, theta, phi) -
                                                                      bfield->Br(r, theta, phi) * vtheta)) /
                      (r * sin(theta))) -
                     ((2. * vphi * vr) / r) + ((2. * vphi * vtheta) / tan(theta));
    return dvphidt;
}

double RungeKutta::velocity(double vr, double vtheta, double vphi, double r, double theta)
{
    double velocity = sqrt((vr * vr) + (r * r * vtheta * vtheta) + (r * r * sin(theta) * sin(theta) * vphi * vphi));
    return velocity;
}

double RungeKutta::gamma(double vr, double vtheta, double vphi, double r, double theta)
{
    double beta = velocity(vr, vtheta, vphi, r, theta) / constants::sc;
    double gamma = 1. / sqrt(1. - beta * beta);
    return gamma;
}

// evaluate one step of the RK integration
// The for loop of the evaluation should be brought into C++ in the near future
std::vector<double> RungeKutta::evaluate(std::vector<double> values)
{
    // pybind11::buffer_info buf = values.request();

    // double *ptr = (double *)buf.ptr;

    // double* rk1 = new double[7];    // stores 1st RungeKutta variables
    // double* rk2 = new double[7];    // stores 2nd RungeKutta variables
    // double* rk3 = new double[7];    // stores 3rd RungeKutta variables
    // double* rk4 = new double[7];    // stores 4th RungeKutta variables

    // double t = ptr[0];
    // double r = ptr[1];
    // double th = ptr[2];
    // double ph = ptr[3];
    // double vr = ptr[4];
    // double vth = ptr[5];
    // double vph = ptr[6];
    // for (double val: values) {
    //     std::cout << val << std::endl;
    // }

    double t = values[0];
    double r = values[1];
    double th = values[2];
    double ph = values[3];
    double vr = values[4];
    double vth = values[5];
    double vph = values[6];

    // std::cout << t << ' ' << r << ' ' << th << ' ' << ph << ' ' << vr << ' ' << vth << ' ' << vph << ' ' << std::endl;

    // delete[] ptr;

    double k1 = h * vr;
    double l1 = h * (vth / r);
    double m1 = h * (vph / (r * sin(th)));
    double a1 = h * dvrdt(r, th, ph, vr, vth, vph);
    double b1 = h * dvthetadt(r, th, ph, vr, vth, vph);
    double c1 = h * dvphidt(r, th, ph, vr, vth, vph);

    // std::cout << k1 << ' ' << l1 << ' ' << m1 << ' ' << a1 << ' ' << b1 << ' ' << c1 << std::endl;

    double k2 = h * (vr + 0.5 * a1);
    double l2 = h * ((vth + 0.5 * b1) / (r + 0.5 * k1));
    double m2 = h * ((vph + 0.5 * c1) / ((r + 0.5 * k1) * (sin(th + 0.5 * l1))));
    double a2 = h * dvrdt(r + 0.5 * k1, th + 0.5 * l1, ph + 0.5 * m1,
                          vr + 0.5 * a1, vth + 0.5 * b1, vph + 0.5 * c1);
    double b2 = h * dvthetadt(r + 0.5 * k1, th + 0.5 * l1, ph + 0.5 * m1,
                              vr + 0.5 * a1, vth + 0.5 * b1, vph + 0.5 * c1);
    double c2 = h * dvphidt(r + 0.5 * k1, th + 0.5 * l1, ph + 0.5 * m1,
                            vr + 0.5 * a1, vth + 0.5 * b1, vph + 0.5 * c1);

    // std::cout << k2 << ' ' << l2 << ' ' << m2 << ' ' << a2 << ' ' << b2 << ' ' << c2 << std::endl;

    double k3 = h * (vr + 0.5 * a2);
    double l3 = h * ((vth + 0.5 * b2) / (r + 0.5 * k2));
    double m3 = h * ((vph + 0.5 * c2) / ((r + 0.5 * k2) * (sin(th + 0.5 * l2))));
    double a3 = h * dvrdt(r + 0.5 * k2, th + 0.5 * l2, ph + 0.5 * m2,
                          vr + 0.5 * a2, vth + 0.5 * b2, vph + 0.5 * c2);
    double b3 = h * dvthetadt(r + 0.5 * k2, th + 0.5 * l2, ph + 0.5 * m2,
                              vr + 0.5 * a2, vth + 0.5 * b2, vph + 0.5 * c2);
    double c3 = h * dvphidt(r + 0.5 * k2, th + 0.5 * l2, ph + 0.5 * m2,
                            vr + 0.5 * a2, vth + 0.5 * b2, vph + 0.5 * c2);

    // std::cout << k3 << ' ' << l3 << ' ' << m3 << ' ' << a3 << ' ' << b3 << ' ' << c3 << std::endl;

    double k4 = h * (vr + a3);
    double l4 = h * ((vth + b3) / (r + k3));
    double m4 = h * ((vph + c3) / ((r + k3) * (sin(th + l3))));
    double a4 = h * dvrdt(r + k3, th + l3, ph + m3, vr + a3, vth + b3,
                          vph + c3);
    double b4 = h * dvthetadt(r + k3, th + l3, ph + m3, vr + a3, vth + b3,
                              vph + c3);
    double c4 = h * dvphidt(r + k3, th + l3, ph + m3, vr + a3, vth + b3,
                            vph + c3);

    // std::cout << k4 << ' ' << l4 << ' ' << m4 << ' ' << a4 << ' ' << b4 << ' ' << c4 << std::endl;

    // std::vector<double> vec(7, 0.);

    // auto nextvals = pybind11::array_t<double>(7);

    // pybind11::buffer_info newbuf = nextvals.request();

    // double* nextptr = (double*) newbuf.ptr;

    values[1] = r + (1. / 6.) * k1 + (1. / 3.) * k2 + (1. / 3.) * k3 + (1. / 6.) * k4;
    values[2] = th + (1. / 6.) * l1 + (1. / 3.) * l2 + (1. / 3.) * l3 + (1. / 6.) * l4;
    values[3] = ph + (1. / 6.) * m1 + (1. / 3.) * m2 + (1. / 3.) * m3 + (1. / 6.) * m4;
    values[4] = vr + (1. / 6.) * a1 + (1. / 3.) * a2 + (1. / 3.) * a3 + (1. / 6.) * a4;
    values[5] = vth + (1. / 6.) * b1 + (1. / 3.) * b2 + (1. / 3.) * b3 + (1. / 6.) * b4;
    values[6] = vph + (1. / 6.) * c1 + (1. / 3.) * c2 + (1. / 3.) * c3 + (1. / 6.) * c4;
    values[0] = t + h;

    //     for (int i=0; i<7; ++i) {
    //     newbuf.ptr[i] = nextptr[i];
    // }

    // delete[] nextptr;

    // std::cout << h << std::endl;

    // // std::cout << vec[0] << std::endl;

    // for (double val: vec) {
    //     std::cout << val << '\t';
    // }

    return values;
}
// pybind11::array_t<double> newarr(7) {t, r, th, ph, vr, vth, vph};

// auto nextvals = pybind11::array_t<double>(7);

// pybind11::buffer_info newbuf = nextvals.request();

// double* nextptr = (double*) newbuf.ptr;

// for (int i=0; i<7; ++i) {
//     newptr[i] = nextptr[i];
// }

// double *newptr = (double *) newbuf.ptr,

// newptr[0] = t;
// newptr[1] = r;
// newptr[2] = th;
// newptr[3] = ph;
// newptr[4] = vr;
// newptr[5] = vth;
// newptr[6] = vph;

// ptr[0] = t;
// ptr[1] = r;
// ptr[2] = th;
// ptr[3] = ph;
// ptr[4] = vr;
// ptr[5] = vth;
// ptr[6] = vph;

// return nextvals;
// return newvals;
// return values;
// rk1[0] = h * vr;
// rk1[1] = h * vth / r;
// rk1[2] = h * vph / (r * sin(th));
// rk1[3] = sqrt((vr*vr) + (r * r* vth*vth) + (r * sin(th) * vph*r * sin(th) * vph));
// rk1[4] = 1. / (sqrt(1 - (rk1[3] / C)*(rk1[3] / C)));
// rk1[5] = h * dvrdt(t, r, th, ph, vr, vth, vph,rk1[4]);
// rk1[6] = h * dvthetadt(t, r, th, ph, vr, vth, vph,rk1[4]);
// rk1[7] = h * dvphidt(t, r, th, ph, vr, vth, vph,rk1[4]);

// rk2[0] = h * vr + 0.5 * rk1[5];
// rk2[1] = h * vth + 0.5 * rk1[6] / (r + 0.5 * rk1[0]);
// rk2[2] = h * vph + 0.5 * rk1[7] / ((r + 0.5 * rk1[0]) * (sin(th + 0.5 * rk1[1])));
// rk2[3] = sqrt((vr + 0.5 * rk1[5])*(vr + 0.5 * rk1[5]) + (r + 0.5 * rk1[0]) *
//                                     (vth + 0.5 * rk1[6])*(r + 0.5 * rk1[0]) *
//                                     (vth + 0.5 * rk1[6]) +
//             (r + 0.5 * rk1[0]) * sin(th + 0.5 * rk1[1]) *
//                 (vph + 0.5 * rk1[7])*(r + 0.5 * rk1[0]) * sin(th + 0.5 * rk1[1]) *
//                 (vph + 0.5 * rk1[7]));
// rk2[4] = 1. / (sqrt(1 - (rk2[3] / C)*(rk2[3] / C)));
// rk2[5] = h * dvrdt(t + 0.5 * h, r + 0.5 * rk1[0], th + 0.5 * rk1[1], ph + 0.5 * rk1[2],
//             vr + 0.5 * rk1[5], vth + 0.5 * rk1[6], vph + 0.5 * rk1[7], rk2[4]);
// rk2[6] = h * dvthetadt(t + 0.5 * h, r + 0.5 * rk1[0], th + 0.5 * rk1[1], ph + 0.5 * rk1[2],
//                 vr + 0.5 * rk1[5], vth + 0.5 * rk1[6], vph + 0.5 * rk1[7], rk2[4]);
// rk2[7] = h * dvphidt(t + 0.5 * h, r + 0.5 * rk1[0], th + 0.5 * rk1[1], ph + 0.5 * rk1[2],
//                 vr + 0.5 * rk1[5], vth + 0.5 * rk1[6], vph + 0.5 * rk1[7], rk2[4]);

// rk3[0] = h * vr + 0.5 * rk2[5];
// rk3[1] = h * vth + 0.5 * rk2[6] / (r + 0.5 * rk2[0]);
// rk3[2] = h * vph + 0.5 * rk2[7] / ((r + 0.5 * rk2[0]) * (sin(th + 0.5 *rk2[1])));
// rk3[3] = sqrt((vr + 0.5 * rk2[5])*(vr + 0.5 * rk2[5]) + ((r + 0.5 * rk2[0]) *
//                                     (vth + 0.5 * rk2[6])*(r + 0.5 * rk2[0]) *
//                                     (vth + 0.5 * rk2[6])) +
//             ((r + 0.5 * rk2[0]) * sin(th + 0.5 *rk2[1]) *
//                 (vph + 0.5 * rk2[7])*(r + 0.5 * rk2[0]) * sin(th + 0.5 *rk2[1]) *
//                 (vph + 0.5 * rk2[7])));
// rk3[4] = 1. / (sqrt(1 - (rk3[3] / C)*(rk3[3] / C)));
// rk3[5] = h * dvrdt(t + 0.5 * h, r + 0.5 * rk2[0], th + 0.5 *rk2[1], ph + 0.5 * rk2[2],
//             vr + 0.5 * rk2[5], vth + 0.5 * rk2[6], vph + 0.5 * rk2[7], rk3[4]);
// rk3[6] = h * dvthetadt(t + 0.5 * h, r + 0.5 * rk2[0], th + 0.5 *rk2[1], ph + 0.5 * rk2[2],
//                 vr + 0.5 * rk2[5], vth + 0.5 * rk2[6], vph + 0.5 * rk2[7], rk3[4]);
// rk3[7] = h * dvphidt(t + 0.5 * h, r + 0.5 * rk2[0], th + 0.5 *rk2[1], ph + 0.5 * rk2[2],
//                 vr + 0.5 * rk2[5], vth + 0.5 * rk2[6], vph + 0.5 * rk2[7], rk3[4]);

// rk4[0] = h * vr + rk3[5];
// rk4[1] = h * vth + rk3[6] / (r + rk3[0]);
// rk4[2] = h * vph + rk3[7] / ((r + rk3[0]) * (sin(th + rk3[1])));
// rk4[3] = sqrt((vr + rk3[5])*(vr + rk3[5]) + ((r + rk3[0]) * (vth + rk3[6])*(r + rk3[0]) * (vth + rk3[6])) +
//             ((r + rk3[0]) * sin(th + rk3[1]) * (vph + rk3[7])*(r + rk3[0]) * sin(th + rk3[1]) * (vph + rk3[7])));
// rk4[4] = 1. / (sqrt(1 - (rk4[3] / C)*(rk4[3] / C)));
// rk4[5] = h * dvrdt(t + h, r + rk3[0], th + rk3[1], ph + rk3[2], vr + rk3[5], vth + rk3[6],
//             vph + rk3[7], rk4[4]);
// rk4[6] = h * dvthetadt(t + h, r + rk3[0], th + rk3[1], ph + rk3[2], vr + rk3[5], vth + rk3[6],
//                 vph + rk3[7], rk4[4]);
// rk4[7] = h * dvphidt(t + h, r + rk3[0], th + rk3[1], ph + rk3[2], vr + rk3[5], vth + rk3[6],
//                 vph + rk3[7], rk4[4]);

// k = (1. / 6.) * rk1[0] + (1. / 3.) * rk2[0] + (1. / 3.) * rk3[0] + (1. / 6.) * rk4[0]
// l = (1. / 6.) * rk1[1] + (1. / 3.) * rk2[1] + (1. / 3.) * rk3[1] + (1. / 6.) * rk4[1]
// m = (1. / 6.) * rk1[2] + (1. / 3.) * rk2[2] + (1. / 3.) * rk3[2] + (1. / 6.) * rk4[2]
// a = (1. / 6.) * rk1[5] + (1. / 3.) * rk2[5] + (1. / 3.) * rk3[5] + (1. / 6.) * rk4[5]
// b = (1. / 6.) * rk1[6] + (1. / 3.) * rk2[6] + (1. / 3.) * rk3[6] + (1. / 6.) * rk4[6]
// c = (1. / 6.) * rk1[7] + (1. / 3.) * rk2[7] + (1. / 3.) * rk3[7] + (1. / 6.) * rk4[7]

// double newarr[7] = {0.,0.,0.,0.,0.,0.,0.};
// double *newptr = newarr;

// auto newvals = pybind11::array_t<double>(7);

// pybind11::buffer_info newbuf = newvals.request();

// double* newptr = (double*) newbuf.ptr;

// newptr[1] = r + (1. / 6.) * rk1[0] + (1. / 3.) * rk2[0] + (1. / 3.) * rk3[0] + (1. / 6.) * rk4[0];
// newptr[2] = th + (1. / 6.) * rk1[1] + (1. / 3.) * rk2[1] + (1. / 3.) * rk3[1] + (1. / 6.) * rk4[1];
// newptr[3] = ph + (1. / 6.) * rk1[2] + (1. / 3.) * rk2[2] + (1. / 3.) * rk3[2] + (1. / 6.) * rk4[2];
// newptr[4] = vr + (1. / 6.) * rk1[5] + (1. / 3.) * rk2[5] + (1. / 3.) * rk3[5] + (1. / 6.) * rk4[5];
// newptr[5] = vth + (1. / 6.) * rk1[6] + (1. / 3.) * rk2[6] + (1. / 3.) * rk3[6] + (1. / 6.) * rk4[6];
// newptr[6] = vph + (1. / 6.) * rk1[7] + (1. / 3.) * rk2[7] + (1. / 3.) * rk3[7] + (1. / 6.) * rk4[7];
// newptr[0] = t + h;

///
// // dvrdt
// const double &RungeKutta::dvrdt(const pybind11::array_t<double> &values, const double &gamma)
// {
//     return (coeff / gamma) * (values[5] * bfield.Bphi(values[1], values[2], values[3]) -
//                               bfield.Btheta(values[1], values[2], values[3]) * values[6]) +
//            (values[5] * values[5]) / values[1] + (values[6] * values[6]) / values[1];
// }

// //dvthetadt
// const double &RungeKutta::dvthetadt(const pybind11::array_t<double> &values, const double &gamma)
// {
//     return (coeff / gamma) * (values[6] * bfield.Br(values[1], values[2], values[3]) -
//                               bfield.Bphi(values[1], values[2], values[3]) * values[4]) -
//            (values[4] * values[5]) / values[1] + (values[6] * values[6]) / (values[1] * tan(values[2]));
// }

// //dvphidt
// const double &RungeKutta::dvphidt(const pybind11::array_t<double> &values, const double &gamma)
// {
//     return (coeff / gamma) * (values[4] * bfield.Btheta(values[1], values[2], values[3]) -
//                               bfield.Br(values[1], values[2], values[3]) * values[5]) -
//            (values[4] * values[6]) / values[1] - (values[5] * values[6]) / (values[1] * tan(values[2]));
// }
//