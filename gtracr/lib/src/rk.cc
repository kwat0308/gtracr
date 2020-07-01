#include <math.h>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // for STL container type conversions

namespace py = pybind11;

// RE = 6371.2;
// g10 = -29404.8 * (1e-9);

const double& B_r(const float &r, const float &theta, const float &phi)
{
    return 2. * (6371.2 / r) * (6371.2 / r) * (6371.2 / r) * -29404.8 * (1e-9) * cos(theta);
}

const double& B_theta(const float &r, const float &theta, const float &phi)
{
    return (6371.2 / r) * (6371.2 / r) * (6371.2 / r) * -29404.8 * (1e-9) * sin(theta);
}

const double& B_phi(const float &r, const float &theta, const float &phi)
{
    return 0.;
}

const double& gamma(const double& vr, const double& vtheta, const double& vphi, const double& r, const double& theta)
{
    const double& v = sqrt((vr*vr) + (r*vtheta*r*vtheta) + (r*sin(theta)*vphi*r*sin(theta)*vphi));
    return 1. / sqrt(1 - (v / 299792458)*(v / 299792458));

}

// # obtained from D.F.Smart, M.A.Shea, Sept. 1, 2004
// # radial velocity DE
const double& dvrdt(const double& t, const double& r, const double& theta, const double& phi, const double& vr, const double& vtheta, const double& vphi, const double& coeff){
    const double& term1 = (coeff / gamma(vr, vtheta, vphi, r, theta)) * (vtheta * B_phi(r, theta, phi) -
                     B_theta(r, theta, phi) * vphi);
    const double& term2 = (vtheta*vtheta) / r;
    const double& term3 = (vphi*vphi) / r;

    return term1 + term2 + term3;
}
    


// # theta component velocity DE
const double& dvthetadt(const double& t, const double& r, const double& theta, const double& phi, const double& vr, const double& vtheta, const double& vphi, const double& coeff){
    const double& term1 = (coeff / gamma(vr, vtheta, vphi, r, theta)) * (vphi * B_r(r, theta, phi) - B_phi(r, theta, phi) * vr);
    const double& term2 = (vr * vtheta) / r;
    const double& term3 = (vphi*vphi) / (r * tan(theta));
    return term1 - term2 + term3;
}

// # phi comp vel/. DE
const double& dvphidt(const double& t, const double& r, const double& theta, const double& phi, const double& vr, const double& vtheta, const double& vphi, const double& coeff){
    const double& term1 = (coeff / gamma(vr, vtheta, vphi, r, theta)) * (B_theta(r, theta, phi) * vr - B_r(r, theta, phi) * vtheta);
    const double& term2 = (vr * vphi) / r;
    const double& term3 = (vtheta * vphi) / (r * tan(theta));
    return term1 - term2 - term3;
}

// # evaluate 4th-order Runge Kutta with 6 coupled differential equations
// # this code can be reduced greatly by removing certain arguments, however by
// # making the code more general those arguments have to stay there.
// # Edit 2: this only performs one iteration of RK
// # the position DEs are replaced with the spherical definition for velocity
std::vector<double> evaluate(const double& mass, const int charge, const double& h, const double& t, const double& r, const double& th, const double& ph, const double& vr, const double& vth, const double& vph){

    const double& coeff = charge / mass;
    // # get the RK terms
    // # k,l,m are for drdt, dthdt, dphdt
    // # p,q,s are for momenta
    const double& k1 = h * vr;
    const double& l1 = h * vth / r;
    const double& m1 = h * vph / (r * sin(th));
    const double& a1 = h * dvrdt(t, r, th, ph, vr, vth, vph, coeff);
    const double& b1 = h * dvthetadt(t, r, th, ph, vr, vth, vph, coeff);
    const double& c1 = h * dvphidt(t, r, th, ph, vr, vth, vph, coeff);

    const double& k2 = h * vr + 0.5 * a1;
    const double& l2 = h * vth + 0.5 * b1 / (r + 0.5 * k1);
    const double& m2 = h * vph + 0.5 * c1 / ((r + 0.5 * k1) * (sin(th + 0.5 * l1)));
    const double& a2 = h * dvrdt(t + 0.5 * h, r + 0.5 * k1, th + 0.5 * l1, ph + 0.5 * m1,
                   vr + 0.5 * a1, vth + 0.5 * b1, vph + 0.5 * c1, coeff);
    const double& b2 = h * dvthetadt(t + 0.5 * h, r + 0.5 * k1, th + 0.5 * l1, ph + 0.5 * m1,
                       vr + 0.5 * a1, vth + 0.5 * b1, vph + 0.5 * c1, coeff);
    const double& c2 = h * dvphidt(t + 0.5 * h, r + 0.5 * k1, th + 0.5 * l1, ph + 0.5 * m1,
                     vr + 0.5 * a1, vth + 0.5 * b1, vph + 0.5 * c1, coeff);

    const double& k3 = h * vr + 0.5 * a2;
    const double& l3 = h * vth + 0.5 * b2 / (r + 0.5 * k2);
    const double& m3 = h * vph + 0.5 * c2 / ((r + 0.5 * k2) * (sin(th + 0.5 * l2)));
    const double& a3 = h * dvrdt(t + 0.5 * h, r + 0.5 * k2, th + 0.5 * l2, ph + 0.5 * m2,
                   vr + 0.5 * a2, vth + 0.5 * b2, vph + 0.5 * c2, coeff);
    const double& b3 = h * dvthetadt(t + 0.5 * h, r + 0.5 * k2, th + 0.5 * l2, ph + 0.5 * m2,
                       vr + 0.5 * a2, vth + 0.5 * b2, vph + 0.5 * c2, coeff);
    const double& c3 = h * dvphidt(t + 0.5 * h, r + 0.5 * k2, th + 0.5 * l2, ph + 0.5 * m2,
                     vr + 0.5 * a2, vth + 0.5 * b2, vph + 0.5 * c2, coeff);

    const double& k4 = h * vr + a3;
    const double& l4 = h * vth + b3 / (r + k3);
    const double& m4 = h * vph + c3 / ((r + k3) * (sin(th + l3)));
    const double& a4 = h * dvrdt(t + h, r + k3, th + l3, ph + m3, vr + a3, vth + b3,
                   vph + c3, coeff);
    const double& b4 = h * dvthetadt(t + h, r + k3, th + l3, ph + m3, vr + a3, vth + b3,
                       vph + c3, coeff);
    const double& c4 = h * dvphidt(t + h, r + k3, th + l3, ph + m3, vr + a3, vth + b3,
                     vph + c3, coeff);

    // # get the weighted sum of each component
    const double& k = (1. / 6.) * k1 + (1. / 3.) * k2 + (1. / 3.) * k3 + (1. / 6.) * k4;
    const double& l = (1. / 6.) * l1 + (1. / 3.) * l2 + (1. / 3.) * l3 + (1. / 6.) * l4;
    const double& m = (1. / 6.) * m1 + (1. / 3.) * m2 + (1. / 3.) * m3 + (1. / 6.) * m4;
    const double& a = (1. / 6.) * a1 + (1. / 3.) * a2 + (1. / 3.) * a3 + (1. / 6.) * a4;
    const double& b = (1. / 6.) * b1 + (1. / 3.) * b2 + (1. / 3.) * b3 + (1. / 6.) * b4;
    const double& c = (1. / 6.) * c1 + (1. / 3.) * c2 + (1. / 3.) * c3 + (1. / 6.) * c4;
    // # iterate by weighted sum of stepsize
    // r = r + k;
    // th = th + l;
    // ph = ph + m;
    // vr = vr + a;
    // vth = vth + b;
    // vph = vph + c;
    // t = t + h;

    std::vector<double> vec {t+h, r+k, th+l, ph+m, vr+a, vth+b, vph+c};

    return vec;

}




// declarations
// const std::vector<double>& evaluate(const double&, const int, const double&, const double&, const double&, const double& , const double&, const double&, const double&, const double&);

// Python binding
PYBIND11_MODULE(runge_kutta, m)
{
    m.doc() = "Runge Kutta Evaluator";
    m.def("evaluate", &evaluate, "evaluates runge_kutta");
    // return m.ptr();
}