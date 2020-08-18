
#include <vector>
#include <iostream>
#include "RungeKutta.h"

int main()
{   
    int max_iter = 10;

    RungeKutta *RKI = new RungeKutta(1, 1., 0.01);

    // std::vector vec = {0., 1., 0., 0., 12., 1., 1.};
    std::vector<double> vec (7,1.);
    vec[0] = 0.;
    vec[1] = 6372.;
    vec[4] = 12.;

    for (double val : vec) {
            std::cout << val << std::endl;
        }


    for (int i=0; i<max_iter; ++i) {
        vec = RKI->evaluate(vec);
        // vec = newvec;
        for (double val : vec) {
            std::cout << val << std::endl;
        }
        
    }
}