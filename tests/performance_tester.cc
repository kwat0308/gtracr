// runs ~300 iterations of the TrajectoryTracer object 
// to evaluate the performance of the C++ portion
// of the package
#include <array>
#include <vector>
#include <string>
#include <map>
#include "constants.h"
#include "TrajectoryTracer.h"

int main()
{
    int max_iter = 250;  // max number of iterations
    // contain everything within a for loop
    // since we imitate whats done in Python but just doing it all in C++
    for (int i=0; i<max_iter; ++i) {
        // initialize trajectory tracer
        // use default (proton, stepsize=1e-5, maxsteps=10000)
        TrajectoryTracer traj_tracer = TrajectoryTracer();

        // define initial values
        std::array<double, 7> initial_values {
            0.,    // t0
            constants::RE, constants::pi / 2., 0.,   // r0, theta0, phi0
            10.*constants::KG_M_S_PER_GEVC, 0., 0.   // pr0, ptheta0, pphi0
        };
        // obtain the result as a map (as it is currently)
        std::map<std::string, std::vector<double> > result_map = traj_tracer.evaluate(initial_values);
    }
}