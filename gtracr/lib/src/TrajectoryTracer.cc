// Runge Kutta integrator class

#include <vector>
#include <map>
#include <string>
#include <array>
#include <math.h>
#include <iostream>
#include "constants.h"
#include "MagneticField.h"
#include "TrajectoryTracer.h"

// Constructor
// default
TrajectoryTracer::TrajectoryTracer()
    : bfield_{MagneticField()}, charge_{1}, mass_{0.938}, escape_radius_{10.*constants::RE}, stepsize_{1e-5}, max_iter_{100000}
{
    particle_escaped_ = false;
}

// Requires the charge and mass of the particle
TrajectoryTracer::TrajectoryTracer(const int charge, const double &mass, const double& escape_radius=10.*constants::RE, const double &stepsize=1e-5, const int max_iter=100000)
    : bfield_{MagneticField()}, charge_{charge}, mass_{mass}, escape_radius_{escape_radius}, stepsize_{stepsize}, max_iter_{max_iter}
{
    particle_escaped_ = false;
}

// copy constructor
TrajectoryTracer::TrajectoryTracer(const TrajectoryTracer &traj_tracer)
    : charge_{traj_tracer.charge_}, mass_{traj_tracer.mass_}, escape_radius_{traj_tracer.escape_radius_}, stepsize_{traj_tracer.stepsize_}, max_iter_{traj_tracer.max_iter_}
{
    bfield_ = traj_tracer.bfield_;
    particle_escaped_ = false;
}

// copy assignment operator
TrajectoryTracer &TrajectoryTracer::operator=(const TrajectoryTracer &traj_tracer)
{
    bfield_ = traj_tracer.bfield_;
    charge_ = traj_tracer.charge_;
    mass_ = traj_tracer.mass_;
    escape_radius_ = traj_tracer.escape_radius_;
    stepsize_ = traj_tracer.stepsize_;
    max_iter_ = traj_tracer.max_iter_;
    particle_escaped_ = false;
    return *this;
}

// ODEs based on relativistic Lorentz force equation with auxiliary terms from acceleration in spherical
// coordinates
// drdt
double TrajectoryTracer::dr_dt(double pr)
{
    return pr;
}

// dtheta dt
double TrajectoryTracer::dtheta_dt(double r, double ptheta)
{
    return ptheta / r;
}
// dphidt
double TrajectoryTracer::dphi_dt(double r, double theta, double pphi)
{
    return pphi / (r * sin(theta));
}
// dvrdt
double TrajectoryTracer::dpr_dt(double r, double theta, double phi, double pr, double ptheta, double pphi)
{   
    double lorentz_term = (charge_ * constants::ELEMENTARY_CHARGE) *  
                            (ptheta * bfield_.Bphi(r, theta, phi) - bfield_.Btheta(r, theta, phi) * pphi);
    // double auxiliary_terms = (r * ptheta * ptheta) + (r * pphi * pphi * sin(theta) * sin(theta));
    double auxiliary_terms = ((ptheta*ptheta) / r) + ((pphi*pphi) / r);
    double dpr_dt =  lorentz_term + auxiliary_terms;
    return dpr_dt;
}

//dpthetadt
double TrajectoryTracer::dptheta_dt(double r, double theta, double phi, double pr, double ptheta, double pphi)
{
    double lorentz_term = (charge_ * constants::ELEMENTARY_CHARGE) *  
                            (pphi * bfield_.Br(r, theta, phi) - bfield_.Bphi(r, theta, phi) * pr);
    // double auxiliary_terms = (pphi * pphi * sin(theta) * cos(theta)) - (2. * pr * ptheta / r);
    double auxiliary_terms = ((pphi*pphi) / (r*tan(theta))) - ((pr*ptheta) / r);
    double dptheta_dt =  lorentz_term + auxiliary_terms;
    return dptheta_dt;
}

//dpphidt
double TrajectoryTracer::dpphi_dt(double r, double theta, double phi, double pr, double ptheta, double pphi)
{
    double lorentz_term = (charge_ * constants::ELEMENTARY_CHARGE) *  
                            (pr * bfield_.Btheta(r, theta, phi) - bfield_.Br(r, theta, phi) * ptheta);
    // double auxiliary_terms = ((2. * pphi * pr) / r) + ((2. * pphi * ptheta) / tan(theta));
    double auxiliary_terms = ((pr*pphi) / r) + ((ptheta*pphi) / (r*tan(theta)));
    double dpphi_dt = lorentz_term - auxiliary_terms;
    return dpphi_dt;
}

// lorentz factor
double TrajectoryTracer::gamma(double pr, double ptheta, double pphi)
{
    double momentum = sqrt((pr * pr) + (ptheta * ptheta) + (pphi* pphi));  // momentum magnitude
    double momentum_ratio = momentum / (mass_ * constants::SPEED_OF_LIGHT);  // p / mc
    double gamma = sqrt(1. + (momentum_ratio * momentum_ratio));
    return gamma;
}

// evaluate one step of the RK integration
// The for loop of the evaluation should be brought into C++ in the near future
// Runge Kutta evaluation is done inversely (inverse charge) since we want to do backtracking
std::array<double, 7>& TrajectoryTracer::rk_step(std::array<double, 7>& vec)
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

    double r_k1 = stepsize_ * dr_dt(pr);
    double theta_k1 = stepsize_ * dtheta_dt(r, ptheta);
    double phi_k1 = stepsize_ * dphi_dt(r, theta, pphi);
    double pr_k1 = stepsize_ * dpr_dt(r, theta, phi, pr, ptheta, pphi);
    double ptheta_k1 = stepsize_ * dptheta_dt(r, theta, phi, pr, ptheta, pphi);
    double pphi_k1 = stepsize_ * dpphi_dt(r, theta, phi, pr, ptheta, pphi);

    // std::cout << r_k1 << ' ' << theta_k1 << ' ' << phi_k1 << ' ' << pr_k1 << ' ' << ptheta_k1 << ' ' << pphi_k1 << std::endl;

    double r_k2 = stepsize_ * dr_dt(pr + 0.5 * pr_k1);
    double theta_k2 = stepsize_ * dtheta_dt(r + 0.5 * r_k1, ptheta + 0.5 * ptheta_k1);
    double phi_k2 = stepsize_ * dphi_dt(r + 0.5 * r_k1, theta + 0.5 * theta_k1, pphi + 0.5 * pphi_k1);
    double pr_k2 = stepsize_ * dpr_dt(r + 0.5 * r_k1, theta + 0.5 * theta_k1, phi + 0.5 * phi_k1,
                          pr + 0.5 * pr_k1, ptheta + 0.5 * ptheta_k1, pphi + 0.5 * pphi_k1);
    double ptheta_k2 = stepsize_ * dptheta_dt(r + 0.5 * r_k1, theta + 0.5 * theta_k1, phi + 0.5 * phi_k1,
                             pr + 0.5 * pr_k1, ptheta + 0.5 * ptheta_k1, pphi + 0.5 * pphi_k1);
    double pphi_k2 = stepsize_ * dpphi_dt(r + 0.5 * r_k1, theta + 0.5 * theta_k1, phi + 0.5 * phi_k1,
                            pr + 0.5 * pr_k1, ptheta + 0.5 * ptheta_k1, pphi + 0.5 * pphi_k1);

    // std::cout << r_k2 << ' ' << theta_k2 << ' ' << phi_k2 << ' ' << pr_k2 << ' ' << ptheta_k2 << ' ' << pphi_k2 << std::endl;

    double r_k3 = stepsize_ * dr_dt(pr + 0.5 * pr_k2);
    double theta_k3 = stepsize_ * dtheta_dt(r + 0.5 * r_k2, ptheta + 0.5 * ptheta_k2);
    double phi_k3 = stepsize_ * dphi_dt(r + 0.5 * r_k2, theta + 0.5 * theta_k2, pphi + 0.5 * pphi_k2);
    double pr_k3 = stepsize_ * dpr_dt(r + 0.5 * r_k2, theta + 0.5 * theta_k2, phi + 0.5 * phi_k2,
                          pr + 0.5 * pr_k2, ptheta + 0.5 * ptheta_k2, pphi + 0.5 * pphi_k2);
    double ptheta_k3 = stepsize_ * dptheta_dt(r + 0.5 * r_k2, theta + 0.5 * theta_k2, phi + 0.5 * phi_k2,
                              pr + 0.5 * pr_k2, ptheta + 0.5 * ptheta_k2, pphi + 0.5 * pphi_k2);
    double pphi_k3 = stepsize_ * dpphi_dt(r + 0.5 * r_k2, theta + 0.5 * theta_k2, phi + 0.5 * phi_k2,
                            pr + 0.5 * pr_k2, ptheta + 0.5 * ptheta_k2, pphi + 0.5 * pphi_k2);

    // std::cout << r_k3 << ' ' << theta_k3 << ' ' << phi_k3 << ' ' << pr_k3 << ' ' << ptheta_k3 << ' ' << pphi_k3 << std::endl;

    double r_k4 = stepsize_ * dr_dt(pr + pr_k3);
    double theta_k4 = stepsize_ * dtheta_dt(r + r_k3, ptheta + ptheta_k3);
    double phi_k4 = stepsize_ * dphi_dt(r + r_k3, theta + theta_k3, pphi + pphi_k3);
    double pr_k4 = stepsize_ * dpr_dt(r + r_k3, theta + theta_k3, phi + phi_k3, pr + pr_k3, ptheta + ptheta_k3,
                          pphi + pphi_k3);
    double ptheta_k4 = stepsize_ * dptheta_dt(r + r_k3, theta + theta_k3, phi + phi_k3, pr + pr_k3, ptheta + ptheta_k3,
                              pphi + pphi_k3);
    double pphi_k4 = stepsize_ * dpphi_dt(r + r_k3, theta + theta_k3, phi + phi_k3, pr + pr_k3, ptheta + ptheta_k3,
                            pphi + pphi_k3);

    vec[1] += (1. / (6. * rel_mass)) * (r_k1 + 2. * r_k2 + 2. * r_k3 + r_k4);
    vec[2] += (1. / (6. * rel_mass)) * (theta_k1 + 2. * theta_k2 + 2. * theta_k3 + theta_k4);
    vec[3] += (1. / (6. * rel_mass)) * (phi_k1 + 2. * phi_k2 + 2. * phi_k3 + phi_k4);
    vec[4] += (1. / (6. * rel_mass)) * (pr_k1 + 2. * pr_k2 + 2. * pr_k3 + pr_k4);
    vec[5] += (1. / (6. * rel_mass)) * (ptheta_k1 + 2. * ptheta_k2 + 2. * ptheta_k3 + ptheta_k4);
    vec[6] += (1. / (6. * rel_mass)) * (pphi_k1 + 2. * pphi_k2 + 2. * pphi_k3 + pphi_k4);
    vec[0] += stepsize_;

    // for (double val:vec) {
    //     std::cout << val << ' ';
    // }
    // std::cout << std::endl;

    return vec;
}


// evaluate the full runge kutta algorithm (i.e. this contains the loop)
// input is a 7-vector that contains the initial values for the integration processs
// in the form [t0, r0, theta0, phi0, pr0, ptheta0, pphi0]
std::map<std::string, std::vector<double> > TrajectoryTracer::evaluate(std::array<double, 7>& initial_values)
{


    // a container that holds the variables throughout each step
    // variables are given in the order [t, r, theta, phi, pr, ptheta, pphi]
    // this can also be replaced with 7 doubles by moving the function rk_step directly
    // into the for loop
    std::array<double, 7> traj_vector = initial_values;

    // start the integration process
    for (int i=0; i < max_iter_; ++i) {

        // append to arrays first
        // to do this we need to convert spherical to cartesian

        // first rename the variables for readability
        // these can probably be const but lets leave that for now
        double t = traj_vector[0];
        double r = traj_vector[1];
        double theta = traj_vector[2];
        double phi = traj_vector[3];
        double pr = traj_vector[4];
        double ptheta = traj_vector[5];
        double pphi = traj_vector[6];


        // convert the coordinates
        double x = r * sin(theta) * cos(phi);
        double y = r * sin(theta) * sin(phi);
        double z = r * cos(theta);

        // convert the momentum
        double px = pr * sin(theta) * cos(phi) 
                        + r * ptheta * cos(theta) * cos(phi) 
                            - r * pphi * sin(theta) * sin(phi);
        double py = pr * sin(theta) * sin(phi) 
                        + r * ptheta * cos(theta) * sin(phi) 
                            + r * pphi * sin(theta) * cos(phi);
        double pz = pr * cos(theta) 
                        - r * ptheta * sin(theta);

        // finally append the values
        
        time_arr.push_back(t);
        x_arr.push_back(x);
        y_arr.push_back(y);
        z_arr.push_back(z);
        px_arr.push_back(px);
        py_arr.push_back(py);
        pz_arr.push_back(pz);

        // evaluate a runge kutta step
        // return the next iteration of values 
        traj_vector = rk_step(traj_vector);

        // break condition depending on value of r
        // this is set based on if particle has "escaped" 
        // or if the particle has reached back to earth
        // i.e. an allowed or forbidden trajectory

        const double &radius = traj_vector[1];

        // an allowed trajectory
        if (radius > constants::RE + escape_radius_) {
            particle_escaped_ = true;
            // std::cout << "Allowed Trajectory!" << std::endl;
            break;
        }
        
        // a forbidden trajectory
        if (radius < constants::RE) {
            // std::cout << "Forbidden Trajectory!" << std::endl;
            break;
        }
        
    }

    // convert the final values of the trajectory into a std::vector
    // to put this into our map
    // dont want the time component, so start from 2nd component of 
    // trajectory vector
    std::vector<double> final_values (traj_vector.begin()+1, traj_vector.end());

    // create a map and return the map
    std::map<std::string, std::vector<double> > value_map = {
        {"t", time_arr},
        {"x", x_arr},
        {"y", y_arr},
        {"z", z_arr},
        {"px", px_arr},
        {"py", py_arr},
        {"pz", pz_arr},
        {"final_values", final_values}
    };

    return value_map;

}