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

// ODEs based on relativistic Lorentz force equation with auxiliary terms from acceleration in spherical
// coordinates
// drdt
double RungeKutta::dr_dt(double pr)
{
    return pr;
}

// dtheta dt
double RungeKutta::dtheta_dt(double r, double ptheta)
{
    return ptheta / r;
}
// dphidt
double RungeKutta::dphi_dt(double r, double theta, double pphi)
{
    return pphi / (r * sin(theta));
}
// dvrdt
double RungeKutta::dpr_dt(double r, double theta, double phi, double pr, double ptheta, double pphi)
{   
    double lorentz_term = (charge_ * constants::ELEMENTARY_CHARGE) *  
                            (ptheta * bfield_.Bphi(r, theta, phi) - bfield_.Btheta(r, theta, phi) * pphi);
    // double auxiliary_terms = (r * ptheta * ptheta) + (r * pphi * pphi * sin(theta) * sin(theta));
    double auxiliary_terms = ((ptheta*ptheta) / r) + ((pphi*pphi) / r);
    double dpr_dt =  lorentz_term + auxiliary_terms;
    return dpr_dt;
}

//dpthetadt
double RungeKutta::dptheta_dt(double r, double theta, double phi, double pr, double ptheta, double pphi)
{
    double lorentz_term = (charge_ * constants::ELEMENTARY_CHARGE) *  
                            (pphi * bfield_.Br(r, theta, phi) - bfield_.Bphi(r, theta, phi) * pr);
    // double auxiliary_terms = (pphi * pphi * sin(theta) * cos(theta)) - (2. * pr * ptheta / r);
    double auxiliary_terms = ((pphi*pphi) / (r*tan(theta))) - ((pr*ptheta) / r);
    double dptheta_dt =  lorentz_term + auxiliary_terms;
    return dptheta_dt;
}

//dpphidt
double RungeKutta::dpphi_dt(double r, double theta, double phi, double pr, double ptheta, double pphi)
{
    double lorentz_term = (charge_ * constants::ELEMENTARY_CHARGE) *  
                            (pr * bfield_.Btheta(r, theta, phi) - bfield_.Br(r, theta, phi) * ptheta);
    // double auxiliary_terms = ((2. * pphi * pr) / r) + ((2. * pphi * ptheta) / tan(theta));
    double auxiliary_terms = ((pr*pphi) / r) + ((ptheta*pphi) / (r*tan(theta)));
    double dpphi_dt = lorentz_term - auxiliary_terms;
    return dpphi_dt;
}

// lorentz factor
double RungeKutta::gamma(double pr, double ptheta, double pphi)
{
    double momentum = sqrt((pr * pr) + (ptheta * ptheta) + (pphi* pphi));  // momentum magnitude
    double momentum_ratio = momentum / (mass_ * constants::SPEED_OF_LIGHT);  // p / mc
    double gamma = sqrt(1. + (momentum_ratio * momentum_ratio));
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
    double pr = vec[4];
    double ptheta = vec[5];
    double pphi = vec[6];

    // evaluate relativistic mass here
    double rel_mass = mass_ * gamma(pr, ptheta, pphi) * constants::KG_PER_GEVC2;

    // std::cout << gamma(pr, ptheta, pphi) << std::endl;

    // std::cout << t << ' ' << r << ' ' << theta << ' ' << phi << ' ' << vr << ' ' << vtheta << ' ' << vphi << ' ' << std::endl;

    // delete[] ptr;

    double r_k1 = h * dr_dt(pr);
    double theta_k1 = h * dtheta_dt(r, ptheta);
    double phi_k1 = h * dphi_dt(r, theta, pphi);
    double pr_k1 = h * dpr_dt(r, theta, phi, pr, ptheta, pphi);
    double ptheta_k1 = h * dptheta_dt(r, theta, phi, pr, ptheta, pphi);
    double pphi_k1 = h * dpphi_dt(r, theta, phi, pr, ptheta, pphi);

    // std::cout << r_k1 << ' ' << theta_k1 << ' ' << phi_k1 << ' ' << pr_k1 << ' ' << ptheta_k1 << ' ' << pphi_k1 << std::endl;

    double r_k2 = h * dr_dt(pr + 0.5 * pr_k1);
    double theta_k2 = h * dtheta_dt(r + 0.5 * r_k1, ptheta + 0.5 * ptheta_k1);
    double phi_k2 = h * dphi_dt(r + 0.5 * r_k1, theta + 0.5 * theta_k1, pphi + 0.5 * pphi_k1);
    double pr_k2 = h * dpr_dt(r + 0.5 * r_k1, theta + 0.5 * theta_k1, phi + 0.5 * phi_k1,
                          pr + 0.5 * pr_k1, ptheta + 0.5 * ptheta_k1, pphi + 0.5 * pphi_k1);
    double ptheta_k2 = h * dptheta_dt(r + 0.5 * r_k1, theta + 0.5 * theta_k1, phi + 0.5 * phi_k1,
                             pr + 0.5 * pr_k1, ptheta + 0.5 * ptheta_k1, pphi + 0.5 * pphi_k1);
    double pphi_k2 = h * dpphi_dt(r + 0.5 * r_k1, theta + 0.5 * theta_k1, phi + 0.5 * phi_k1,
                            pr + 0.5 * pr_k1, ptheta + 0.5 * ptheta_k1, pphi + 0.5 * pphi_k1);

    // std::cout << r_k2 << ' ' << theta_k2 << ' ' << phi_k2 << ' ' << pr_k2 << ' ' << ptheta_k2 << ' ' << pphi_k2 << std::endl;

    double r_k3 = h * dr_dt(pr + 0.5 * pr_k2);
    double theta_k3 = h * dtheta_dt(r + 0.5 * r_k2, ptheta + 0.5 * ptheta_k2);
    double phi_k3 = h * dphi_dt(r + 0.5 * r_k2, theta + 0.5 * theta_k2, pphi + 0.5 * pphi_k2);
    double pr_k3 = h * dpr_dt(r + 0.5 * r_k2, theta + 0.5 * theta_k2, phi + 0.5 * phi_k2,
                          pr + 0.5 * pr_k2, ptheta + 0.5 * ptheta_k2, pphi + 0.5 * pphi_k2);
    double ptheta_k3 = h * dptheta_dt(r + 0.5 * r_k2, theta + 0.5 * theta_k2, phi + 0.5 * phi_k2,
                              pr + 0.5 * pr_k2, ptheta + 0.5 * ptheta_k2, pphi + 0.5 * pphi_k2);
    double pphi_k3 = h * dpphi_dt(r + 0.5 * r_k2, theta + 0.5 * theta_k2, phi + 0.5 * phi_k2,
                            pr + 0.5 * pr_k2, ptheta + 0.5 * ptheta_k2, pphi + 0.5 * pphi_k2);

    // std::cout << r_k3 << ' ' << theta_k3 << ' ' << phi_k3 << ' ' << pr_k3 << ' ' << ptheta_k3 << ' ' << pphi_k3 << std::endl;

    double r_k4 = h * dr_dt(pr + pr_k3);
    double theta_k4 = h * dtheta_dt(r + r_k3, ptheta + ptheta_k3);
    double phi_k4 = h * dphi_dt(r + r_k3, theta + theta_k3, pphi + pphi_k3);
    double pr_k4 = h * dpr_dt(r + r_k3, theta + theta_k3, phi + phi_k3, pr + pr_k3, ptheta + ptheta_k3,
                          pphi + pphi_k3);
    double ptheta_k4 = h * dptheta_dt(r + r_k3, theta + theta_k3, phi + phi_k3, pr + pr_k3, ptheta + ptheta_k3,
                              pphi + pphi_k3);
    double pphi_k4 = h * dpphi_dt(r + r_k3, theta + theta_k3, phi + phi_k3, pr + pr_k3, ptheta + ptheta_k3,
                            pphi + pphi_k3);

    vec[1] += (1. / (6. * rel_mass)) * (r_k1 + 2. * r_k2 + 2. * r_k3 + r_k4);
    vec[2] += (1. / (6. * rel_mass)) * (theta_k1 + 2. * theta_k2 + 2. * theta_k3 + theta_k4);
    vec[3] += (1. / (6. * rel_mass)) * (phi_k1 + 2. * phi_k2 + 2. * phi_k3 + phi_k4);
    vec[4] += (1. / (6. * rel_mass)) * (pr_k1 + 2. * pr_k2 + 2. * pr_k3 + pr_k4);
    vec[5] += (1. / (6. * rel_mass)) * (ptheta_k1 + 2. * ptheta_k2 + 2. * ptheta_k3 + ptheta_k4);
    vec[6] += (1. / (6. * rel_mass)) * (pphi_k1 + 2. * pphi_k2 + 2. * pphi_k3 + pphi_k4);
    vec[0] += h;

    // vec[1] = r + h * (vr);
    // vec[2] = theta + h * (ptheta / r);
    // vec[3] = phi + h * (vphi / (r*sin(theta)));
    // vec[4] = vr + h * dvr_dt(r, theta, phi, vr, vtheta, vphi);
    // vec[5] = vtheta + h * dvtheta_dt(r, theta, phi, vr, vtheta, vphi);
    // vec[6] = vphi + h * dvphi_dt(r, theta, phi, vr, vtheta, vphi);
    // vec[0] = t + h;


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
// rr_k1[0] = h * vr;
// rr_k1[1] = h * vth / r;
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