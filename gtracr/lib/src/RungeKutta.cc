// Runge Kutta integrator class

// #include <pybind11/pybind11.h>
// #include <pybind11/numpy.h>
// #include <pybind11/stl.h>
// #include <vector>
#include <array>
#include <math.h>
#include <iostream>
#include "constants.h"
#include "MagneticField.h"
#include "RungeKutta.h"

// Constructor
// default
RungeKutta::RungeKutta()
    : bfield_{MagneticField()}, charge_{1}, mass_{0.938}, h{0.01}
{
}

// Requires the charge and mass of the particle
RungeKutta::RungeKutta(const int charge, const double &mass)
    : bfield_{MagneticField()}, charge_{charge}, mass_{mass}, h{0.01}
{
}

// Constructor
// Requires the charge and mass of the particle and stepsize
RungeKutta::RungeKutta(const int charge, const double &mass, const double &step_size)
    : bfield_{MagneticField()}, charge_{charge}, mass_{mass}, h{step_size}
{
}

// copy constructor
RungeKutta::RungeKutta(const RungeKutta &rk)
    : charge_{rk.charge_}, mass_{rk.mass_}, h{rk.h}
{
    bfield_ = rk.bfield_;
}

// copy assignment operator
RungeKutta &RungeKutta::operator=(const RungeKutta &rk)
{
    bfield_ = rk.bfield_;
    charge_ = rk.charge_;
    mass_ = rk.mass_;
    h = rk.h;
    return *this;
}

// ODEs, inputs are given as numpy arrays of the form (time, r, theta, phi, vr, vtheta, vphi)
// dvrdt
double RungeKutta::dvr_dt(double r, double theta, double phi, double vr, double vtheta, double vphi)
{   
    double charge_per_mass = ((-charge_ * constants::ELEMENTARY_CHARGE) / 
                                    (mass_* constants::KG_PER_GEVC2 * gamma(vr, vtheta, vphi, r, theta)));
    double lorentz_term = charge_per_mass * 
                            (vtheta * bfield_.Bphi(r, theta, phi) - bfield_.Btheta(r, theta, phi) * vphi);
    // double accel_term = (r * vtheta * vtheta) + (r * vphi * vphi * sin(theta) * sin(theta));
    double accel_term = (vtheta*vtheta / r) + (vphi*vphi / r);
    double dvr_dt =  lorentz_term + accel_term;
    return dvr_dt;
}

//dvthetadt
double RungeKutta::dvtheta_dt(double r, double theta, double phi, double vr, double vtheta, double vphi)
{
    double charge_per_mass = ((-charge_ * constants::ELEMENTARY_CHARGE) / 
                                    (mass_* constants::KG_PER_GEVC2 * gamma(vr, vtheta, vphi, r, theta)));
    double lorentz_term = charge_per_mass * 
                            (vphi * bfield_.Br(r, theta, phi) - bfield_.Bphi(r, theta, phi) * vr);
    // double accel_term = (vphi * vphi * sin(theta) * cos(theta)) - (2. * vr * vtheta / r);
    double accel_term = ((vphi*vphi) / (r*tan(theta))) - (vr*vtheta / r);
    double dvtheta_dt =  lorentz_term + accel_term;
    return dvtheta_dt;
}

//dvphidt
double RungeKutta::dvphi_dt(double r, double theta, double phi, double vr, double vtheta, double vphi)
{
    double charge_per_mass = ((-charge_ * constants::ELEMENTARY_CHARGE) / 
                                    (mass_* constants::KG_PER_GEVC2 * gamma(vr, vtheta, vphi, r, theta)));
    double lorentz_term = charge_per_mass * 
                            (vr * bfield_.Btheta(r, theta, phi) - bfield_.Br(r, theta, phi) * vtheta);
    // double accel_term = ((2. * vphi * vr) / r) + ((2. * vphi * vtheta) / tan(theta));
    double accel_term = (vr*vphi / r) + ((vtheta*vphi) / (r*tan(theta)));
    double dvphi_dt = lorentz_term - accel_term;
    return dvphi_dt;
}

// magnitude of velocity from spherical coordinates
double RungeKutta::velocity(double vr, double vtheta, double vphi, double r, double theta)
{
    double velocity = sqrt((vr * vr) + (r * r * vtheta * vtheta) + (r * r * sin(theta) * sin(theta) * vphi * vphi));
    return velocity;
}

// lorentz factor
double RungeKutta::gamma(double vr, double vtheta, double vphi, double r, double theta)
{
    double beta = velocity(vr, vtheta, vphi, r, theta) / constants::SPEED_OF_LIGHT;
    double gamma = 1. / sqrt(1. - beta * beta);
    return gamma;
}

// evaluate one step of the RK integration
// The for loop of the evaluation should be brought into C++ in the near future
// Runge Kutta evaluation is done inversely (reverse charge) since we want to do backtracking
std::array<double, 7>& RungeKutta::evaluate(std::array<double, 7>& vec)
{

    double t = vec[0];
    double r = vec[1];
    double theta = vec[2];
    double phi = vec[3];
    double vr = vec[4];
    double vtheta = vec[5];
    double vphi = vec[6];

    // std::cout << t << ' ' << r << ' ' << theta << ' ' << phi << ' ' << vr << ' ' << vtheta << ' ' << vphi << ' ' << std::endl;

    // delete[] ptr;

    // double k1 = vr;
    // double l1 = (vtheta / r);
    // double m1 = (vphi / (r * sin(theta)));
    // double a1 = dvr_dt(r, theta, phi, vr, vtheta, vphi);
    // double b1 = dvtheta_dt(r, theta, phi, vr, vtheta, vphi);
    // double c1 = dvphi_dt(r, theta, phi, vr, vtheta, vphi);

    // // std::cout << k1 << ' ' << l1 << ' ' << m1 << ' ' << a1 << ' ' << b1 << ' ' << c1 << std::endl;

    // double k2 = (vr + h * 0.5 * a1);
    // double l2 = ((vtheta + h * 0.5 * b1) / (r + h * 0.5 * k1));
    // double m2 = ((vphi + h * 0.5 * c1) / ((r + h * 0.5 * k1) * (sin(theta + h * 0.5 * l1))));
    // double a2 = dvr_dt(r + h * 0.5 * k1, theta + h * 0.5 * l1, phi + h * 0.5 * m1,
    //                       vr + h * 0.5 * a1, vtheta + h * 0.5 * b1, vphi + h * 0.5 * c1);
    // double b2 = dvtheta_dt(r + h * 0.5 * k1, theta + h * 0.5 * l1, phi + h * 0.5 * m1,
    //                           vr + h * 0.5 * a1, vtheta + h * 0.5 * b1, vphi + h * 0.5 * c1);
    // double c2 = dvphi_dt(r + h * 0.5 * k1, theta + h * 0.5 * l1, phi + h * 0.5 * m1,
    //                         vr + h * 0.5 * a1, vtheta + h * 0.5 * b1, vphi + h * 0.5 * c1);

    // // std::cout << k2 << ' ' << l2 << ' ' << m2 << ' ' << a2 << ' ' << b2 << ' ' << c2 << std::endl;

    // double k3 = (vr + h * 0.5 * a2);
    // double l3 = ((vtheta + h * 0.5 * b2) / (r + h * 0.5 * k2));
    // double m3 = ((vphi + h * 0.5 * c2) / ((r + h * 0.5 * k2) * (sin(theta + h * 0.5 * l2))));
    // double a3 = dvr_dt(r + h * 0.5 * k2, theta + h * 0.5 * l2, phi + h * 0.5 * m2,
    //                       vr + h * 0.5 * a2, vtheta + h * 0.5 * b2, vphi + h * 0.5 * c2);
    // double b3 = dvtheta_dt(r + h * 0.5 * k2, theta + h * 0.5 * l2, phi + h * 0.5 * m2,
    //                           vr + h * 0.5 * a2, vtheta + h * 0.5 * b2, vphi + h * 0.5 * c2);
    // double c3 = dvphi_dt(r + h * 0.5 * k2, theta + h * 0.5 * l2, phi + h * 0.5 * m2,
    //                         vr + h * 0.5 * a2, vtheta + h * 0.5 * b2, vphi + h * 0.5 * c2);

    // // std::cout << k3 << ' ' << l3 << ' ' << m3 << ' ' << a3 << ' ' << b3 << ' ' << c3 << std::endl;

    // double k4 = (vr + h * a3);
    // double l4 = ((vtheta + h * b3) / (r + h * k3));
    // double m4 = ((vphi + h * c3) / ((r + h * k3) * (sin(theta + h * l3))));
    // double a4 = dvr_dt(r + h * k3, theta + h * l3, phi + h * m3, vr + h * a3, vtheta + h * b3,
    //                       vphi + h * c3);
    // double b4 = dvtheta_dt(r + h * k3, theta + h * l3, phi + h * m3, vr + h * a3, vtheta + h * b3,
    //                           vphi + h * c3);
    // double c4 = dvphi_dt(r + h * k3, theta + h * l3, phi + h * m3, vr + h * a3, vtheta + h * b3,
    //                         vphi + h * c3);

    // // vec[1] = r + (h / 6.) * (k1 + 2. * k2 + 2. * k3 + k4);
    // // vec[2] = theta + (h / 6.) * (l1 + 2. * l2 + 2. * l3 + l4);
    // // vec[3] = phi + (h / 6.) * (m1 + 2. * m2 + 2. * m3 + m4);
    // // vec[4] = vr + (h / 6.) * (a1 + 2. * a2 + 2. * a3 + a4);
    // // vec[5] = vtheta + (h / 6.) * (b1 + 2. * b2 + 2. * b3 + b4);
    // // vec[6] = vphi + (h / 6.) * (c1 + 2. * c2 + 2. * c3 + c4);
    // vec[1] += (h / 6.) * (k1 + 2. * k2 + 2. * k3 + k4);
    // vec[2] += (h / 6.) * (l1 + 2. * l2 + 2. * l3 + l4);
    // vec[3] += (h / 6.) * (m1 + 2. * m2 + 2. * m3 + m4);
    // vec[4] += (h / 6.) * (a1 + 2. * a2 + 2. * a3 + a4);
    // vec[5] += (h / 6.) * (b1 + 2. * b2 + 2. * b3 + b4);
    // vec[6] += (h / 6.) * (c1 + 2. * c2 + 2. * c3 + c4);
    // vec[0] += h;

    vec[1] = r + h * (vr);
    vec[2] = theta + h * (vtheta / r);
    vec[3] = phi + h * (vphi / (r*sin(theta)));
    vec[4] = vr + h * dvr_dt(r, theta, phi, vr, vtheta, vphi);
    vec[5] = vtheta + h * dvtheta_dt(r, theta, phi, vr, vtheta, vphi);
    vec[6] = vphi + h * dvphi_dt(r, theta, phi, vr, vtheta, vphi);
    vec[0] = t + h;


    // for (double val:vec) {
    //     std::cout << val << ' ';
    // }
    // std::cout << std::endl;

    return vec;
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
//     return (charge_per_mass / gamma) * (values[5] * bfield.Bphi(values[1], values[2], values[3]) -
//                               bfield.Btheta(values[1], values[2], values[3]) * values[6]) +
//            (values[5] * values[5]) / values[1] + (values[6] * values[6]) / values[1];
// }

// //dvthetadt
// const double &RungeKutta::dvthetadt(const pybind11::array_t<double> &values, const double &gamma)
// {
//     return (charge_per_mass / gamma) * (values[6] * bfield.Br(values[1], values[2], values[3]) -
//                               bfield.Bphi(values[1], values[2], values[3]) * values[4]) -
//            (values[4] * values[5]) / values[1] + (values[6] * values[6]) / (values[1] * tan(values[2]));
// }

// //dvphidt
// const double &RungeKutta::dvphidt(const pybind11::array_t<double> &values, const double &gamma)
// {
//     return (charge_per_mass / gamma) * (values[4] * bfield.Btheta(values[1], values[2], values[3]) -
//                               bfield.Br(values[1], values[2], values[3]) * values[5]) -
//            (values[4] * values[6]) / values[1] - (values[5] * values[6]) / (values[1] * tan(values[2]));
// }
//