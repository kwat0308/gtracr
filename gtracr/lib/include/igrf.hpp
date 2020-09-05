#ifndef __IGRF_HPP_
#define __IGRF_HPP_

#include <array>
#include <fstream>
#include <iostream>
#include <sstream>
// #include <cstdio>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>

#include "MagneticField.hpp"

using json = nlohmann::json;

namespace igrf_const {
// flag for external variables, not used
constexpr int IEXT = 0;
// ?
constexpr int RECL = 81;
// ?
constexpr int MAXINBUFF = RECL + 14;
// ?
constexpr int MAXREAD = MAXINBUFF - 2;
// max number of models
constexpr int MAXMOD = 25;
// max path and filename length
constexpr int PATHLEN = MAXREAD;

constexpr int MAXDEG = 13;
constexpr int MAXCOEFF = (MAXDEG * (MAXDEG + 2) + 1);

constexpr double MAXEPOCH = 2020.00;

};  // namespace igrf_const

/*
The IGRF model of Earth's magnetic field.

This code is derived from the MagneticField class that contains the Earth's
magnetic field as a dipole field.

This code supports dates (year, month, and day) from 1900 to 2025 by utilizing
interpolation and extrapolation methods with the spherical harmonic coefficients
obtained from IGRF13.COF in the data/ directory.

This code is derived from the geomag 7.0 software from NGDC (see geomag70.c
for the original code). The code is obtained from here :
https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html

Class Members
--------------
- sdate_ (double) :
    The starting date for the magnetic field model
- nmax_ (int) :
    The number of truncation for the Taylor series expansion of the evaluation
- model_index (int) :
    The index correlating to the model closest to the provided sdate.
of the magnetic field.
- bfield_ (struct of doubles) :
    The main field values in Cartesian coordinates.
*/

class IGRF : public MagneticField {
 private:
  // members
  double sdate_;    // start date
  int nmax_;        // number of truncation for specific model
  int model_index;  // index of the model closest to sdate_

  // key members for interpolation
  double epoch1_;
  double epoch2_;
  int nmain1_;
  int nmain2_;
  int nsv1_;
  int nsv2_;

  // containers for coefficients
  std::array<double, igrf_const::MAXCOEFF> gh1_arr;   // coeff of first model
  std::array<double, igrf_const::MAXCOEFF> gh2_arr;  // coeff of second model
  std::array<double, igrf_const::MAXCOEFF> ghsv1_arr;   // coeff of first sv model
  std::array<double, igrf_const::MAXCOEFF> ghsv2_arr;  // coeff of second sv model
  std::array<double, igrf_const::MAXCOEFF> gh_arr;     // coefficients
  std::array<double, igrf_const::MAXCOEFF>
      ghsv_arr;  // secular variation coeffs

  // bfield values in cartesian coordiantes
  struct {
    double x;  // northern
    double y;  // eastern
    double z;  // vertically down
  } bfield_;

  // temporary bfield value for ROC
  struct {
    double x;
    double y;
    double z;
  } bfield_temp_;

  // rate of change of bfield
  struct {
    double xdot;
    double ydot;
    double zdot;
  } bfield_sv_;

  // control variables
  int nmodel;             // number of models
  int igdgc;              // flag to choose between geodesic vs geocentric
  double minyr, maxyr;    // min / max year (what is this used for?)
  double minalt, maxalt;  // min / max alt (what is this used for?)

  /****************************************************************************/
  /*                                                                          */
  /*                           Subroutine getshc                              */
  /*                                                                          */
  /****************************************************************************/
  /*                                                                          */
  /*     Reads spherical harmonic coefficients from the specified             */
  /*     model into an array.                                                 */
  /*                                                                          */
  /*     Input:                                                               */
  /*           stream     - Logical unit number                               */
  /*           iflag      - Flag for SV equal to ) or not equal to 0          */
  /*                        for designated read statements                    */
  /*           strec      - Starting record number to read from model         */
  /*           nmax_of_gh - Maximum degree and order of model                 */
  /*                                                                          */
  /*     Output:                                                              */
  /*           gh1 or 2   - Schmidt quasi-normal internal spherical           */
  /*                        harmonic coefficients                             */
  /*                                                                          */
  /*     FORTRAN                                                              */
  /*           Bill Flanagan                                                  */
  /*           NOAA CORPS, DESDIS, NGDC, 325 Broadway, Boulder CO.  80301     */
  /*                                                                          */
  /*     C                                                                    */
  /*           C. H. Shaffer                                                  */
  /*           Lockheed Missiles and Space Company, Sunnyvale CA              */
  /*           August 15, 1988                                                */
  /*                                                                          */
  /****************************************************************************/
  void getshc(const std::string &fname);
  /****************************************************************************/
  /*                                                                          */
  /*                           Subroutine interpsh                            */
  /*                                                                          */
  /****************************************************************************/
  /*                                                                          */
  /*     Interpolates linearly, in time, between two spherical harmonic       */
  /*     models.                                                              */
  /*                                                                          */
  /*     Input:                                                               */
  /*           date     - date of resulting model (in decimal year)           */
  /*           dte1     - date of earlier model                               */
  /*           nmax1    - maximum degree and order of earlier model           */
  /*           gh1      - Schmidt quasi-normal internal spherical             */
  /*                      harmonic coefficients of earlier model              */
  /*           dte2     - date of later model                                 */
  /*           nmax2    - maximum degree and order of later model             */
  /*           gh2      - Schmidt quasi-normal internal spherical             */
  /*                      harmonic coefficients of internal model             */
  /*                                                                          */
  /*     Output:                                                              */
  /*           gha or b - coefficients of resulting model                     */
  /*           nmax     - maximum degree and order of resulting model         */
  /*                                                                          */
  /*     FORTRAN                                                              */
  /*           A. Zunde                                                       */
  /*           USGS, MS 964, box 25046 Federal Center, Denver, CO.  80225     */
  /*                                                                          */
  /*     C                                                                    */
  /*           C. H. Shaffer                                                  */
  /*           Lockheed Missiles and Space Company, Sunnyvale CA              */
  /*           August 17, 1988                                                */
  /*                                                                          */
  /****************************************************************************/
  void interpsh(double date, int gh);
  /****************************************************************************/
  /*                                                                          */
  /*                           Subroutine extrapsh                            */
  /*                                                                          */
  /****************************************************************************/
  /*                                                                          */
  /*     Extrapolates linearly a spherical harmonic model with a              */
  /*     rate-of-change model.                                                */
  /*                                                                          */
  /*     Input:                                                               */
  /*           date     - date of resulting model (in decimal year)           */
  /*           dte1     - date of base model                                  */
  /*           nmax1    - maximum degree and order of base model              */
  /*           gh1      - Schmidt quasi-normal internal spherical             */
  /*                      harmonic coefficients of base model                 */
  /*           nmax2    - maximum degree and order of rate-of-change model    */
  /*           gh2      - Schmidt quasi-normal internal spherical             */
  /*                      harmonic coefficients of rate-of-change model       */
  /*                                                                          */
  /*     Output:                                                              */
  /*           gha or b - Schmidt quasi-normal internal spherical             */
  /*                    harmonic coefficients                                 */
  /*           nmax   - maximum degree and order of resulting model           */
  /*                                                                          */
  /*     FORTRAN                                                              */
  /*           A. Zunde                                                       */
  /*           USGS, MS 964, box 25046 Federal Center, Denver, CO.  80225     */
  /*                                                                          */
  /*     C                                                                    */
  /*           C. H. Shaffer                                                  */
  /*           Lockheed Missiles and Space Company, Sunnyvale CA              */
  /*           August 16, 1988                                                */
  /*                                                                          */
  /****************************************************************************/
  void extrapsh(double date, int gh);
  // not too sure if we need functions below
  /****************************************************************************/
  /*                                                                          */
  /*                           Subroutine dihf                                */
  /*                                                                          */
  /****************************************************************************/
  /*                                                                          */
  /*     Computes the geomagnetic d, i, h, and f from x, y, and z.            */
  /*                                                                          */
  /*     Input:                                                               */
  /*           x  - northward component                                       */
  /*           y  - eastward component                                        */
  /*           z  - vertically-downward component                             */
  /*                                                                          */
  /*     Output:                                                              */
  /*           d  - declination                                               */
  /*           i  - inclination                                               */
  /*           h  - horizontal intensity                                      */
  /*           f  - total intensity                                           */
  /*                                                                          */
  /*     FORTRAN                                                              */
  /*           A. Zunde                                                       */
  /*           USGS, MS 964, box 25046 Federal Center, Denver, CO.  80225     */
  /*                                                                          */
  /*     C                                                                    */
  /*           C. H. Shaffer                                                  */
  /*           Lockheed Missiles and Space Company, Sunnyvale CA              */
  /*           August 22, 1988                                                */
  /*                                                                          */
  /****************************************************************************/
  // std::array<double, 4> dihf(int gh);
  // julday (convert day/month/year to decimal day)
  // double julday();

 public:
  /*
   Initialize the IGRF model. The coefficients are imported and are stored
   in the corresponding arrays.

   Parameters
   -----------
   - fname (std::string) :
         The path to the .COF file.
     - sdate (double) :
         The current year, month, and day in decimal days
   */
  IGRF(const std::string &fname, const double sdate);
  /****************************************************************************/
  /*                                                                          */
  /*                           Subroutine shval3                              */
  /*                                                                          */
  /****************************************************************************/
  /*                                                                          */
  /*     Calculates field components from spherical harmonic (sh)             */
  /*     models.                                                              */
  /*                                                                          */
  /*     Input:                                                               */
  /*           igdgc     - indicates coordinate system used; set equal        */
  /*                       to 1 if geodetic, 2 if geocentric                  */
  /*           latitude  - north latitude, in degrees                         */
  /*           longitude - east longitude, in degrees                         */
  /*           elev      - WGS84 altitude above ellipsoid (igdgc=1), or       */
  /*                       radial distance from earth's center (igdgc=2)      */
  /*           a2,b2     - squares of semi-major and semi-minor axes of       */
  /*                       the reference spheroid used for transforming       */
  /*                       between geodetic and geocentric coordinates        */
  /*                       or components                                      */
  /*           nmax      - maximum degree and order of coefficients           */
  /*           iext      - external coefficients flag (=0 if none)            */
  /*           ext1,2,3  - the three 1st-degree external coefficients         */
  /*                       (not used if iext = 0)                             */
  /*                                                                          */
  /*     Output:                                                              */
  /*           x         - northward component                                */
  /*           y         - eastward component                                 */
  /*           z         - vertically-downward component                      */
  /*                                                                          */
  /*     based on subroutine 'igrf' by D. R. Barraclough and S. R. C. Malin,  */
  /*     report no. 71/1, institute of geological sciences, U.K.              */
  /*                                                                          */
  /*     FORTRAN                                                              */
  /*           Norman W. Peddie                                               */
  /*           USGS, MS 964, box 25046 Federal Center, Denver, CO.  80225     */
  /*                                                                          */
  /*     C                                                                    */
  /*           C. H. Shaffer                                                  */
  /*           Lockheed Missiles and Space Company, Sunnyvale CA              */
  /*           August 17, 1988                                                */
  /*                                                                          */
  /****************************************************************************/
  void shval3(int igdgc, double flat, double flon, double elev, int gh,
              int iext = 0, int ext1 = 0, int ext2 = 0, int ext3 = 0);
  /*
  Obtain the values of the IGRF magnetic field components from the 3-vector in
  spherical coordinates.

  Parameters
  -----------
  - r (const double&) : the radial component of the 3-vector in km
  - theta (const double&) : the polar component of the 3-vector
  - phi (const double&) : the azimuthal component of the 3-vector

  Returns
  -------
  - values (std::array<double, 3>) :
      Array containing the radial, polar, and azimuthal components of the
      magnetic field.
  */
  std::array<double, 3> values(const double &r, const double &theta,
                               const double &phi);
  /*
  The starting date of the model in decimal days
  */
  inline const double &sdate() { return sdate_; }
  /*
  The number of truncation for the chosen model
  */
  inline const int nmax() { return nmax_; }
  /*
  The Northern, Eastern, and Vertical components (i.e. the Cartesian components)
   of the resulting magnetic field at the provided (r, theta, phi) values in
   IGRF.values(). Mainly used for debugging purposes.
  */
  inline const std::array<double, 3> &cartesian_values() {
    return std::array<double, 3>{bfield_.x, bfield_.y, bfield_.z};
  }
};  // IGRF

#endif  // __IGRF_HPP_