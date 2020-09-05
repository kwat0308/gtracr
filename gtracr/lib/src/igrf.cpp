#include "igrf.hpp"

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
IGRF::IGRF(const std::string& fname, const double sdate)
    : model_index{0}, nmodel{24}, igdgc{3}, sdate_{sdate} {
  readmdfile(fname);

  /* Pick model */
  for (int i = 0; i < nmodel; i++) {
    if (sdate_ < yrmax_arr[i]) {
      model_index = i;
      break;
    }
    if (i == nmodel) {
      model_index = nmodel - 1; /* if beyond end of last model use last model */
    }
  }

  /* Get altitude min and max for selected model. */
  minalt = altmin_arr[model_index];
  maxalt = altmax_arr[model_index];

  // for (auto val : irec_posarr) {
  //   std::cout << val << std::endl;
  // }

  // get the spherical harmonic coeffiecients for the specific date
  // and either interpolate or extrapolate if we have to
  if (max2_arr[model_index] == 0)  // if model isnt the last model in the file
  {
    getshc(fname, 1, irec_posarr[model_index], max1_arr[model_index], 1);
    getshc(fname, 1, irec_posarr[model_index + 1], max1_arr[model_index + 1],
           2);
    interpsh(sdate_, yrmin_arr[model_index], max1_arr[model_index],
             yrmin_arr[model_index + 1], max1_arr[model_index + 1], 3);
    interpsh(sdate_ + 1, yrmin_arr[model_index], max1_arr[model_index],
             yrmin_arr[model_index + 1], max1_arr[model_index + 1], 4);
  } else {
    getshc(fname, 1, irec_posarr[model_index], max1_arr[model_index], 1);
    getshc(fname, 0, irec_posarr[model_index], max2_arr[model_index], 2);
    extrapsh(sdate_, epoch_arr[model_index], max1_arr[model_index],
             max2_arr[model_index], 3);
    extrapsh(sdate_ + 1, epoch_arr[model_index], max1_arr[model_index],
             max2_arr[model_index], 4);
  }
}

/*
Read the file from fname and store the model information into
the arrays.

Parameters:
-----------
- fname (std::string) :
  path to the .COF file
*/
void IGRF::readmdfile(const std::string& fname) {
  std::ifstream ifs{fname};
  if (!ifs) {
    throw std::runtime_error("Cannot open file!");
  }

  std::string line;

  // some temporary variables for writing the variables to
  // the storage arrays
  int max1, max2, max3;
  std::string model;
  double epoch, yrmin, yrmax, altmin, altmax;

  int modelI = 0;

  while (std::getline(ifs, line)) {  // read each line in file
    // if second entry has no value, then its on model line
    if (line[1] != ' ') {
      continue;
    } else {
      std::istringstream ss{line};        // stringstream of each row
      irec_posarr[modelI] = ifs.tellg();  // position in function
      // same thing as using sscanf, but this is kind of ugly...
      while (ss >> model >> epoch >> max1 >> max2 >> max3 >> yrmin >> yrmax >>
             altmin >> altmax) {
        // std::cout << max1 << std::endl;
        model_arr[modelI] = model;
        epoch_arr[modelI] = epoch;
        max1_arr[modelI] = max1;
        max2_arr[modelI] = max2;
        max3_arr[modelI] = max3;
        yrmin_arr[modelI] = yrmin;
        yrmax_arr[modelI] = yrmax;
        altmin_arr[modelI] = altmin;
        altmax_arr[modelI] = altmax;
      }  // while (ss >> values)
    }    // if (line[1] != ' ')

    /* Compute date range for all models */
    if (modelI == 0) /*If first model */
    {
      minyr = yrmin_arr[0];
      maxyr = yrmax_arr[0];
    } else {
      if (yrmin_arr[modelI] < minyr) {
        minyr = yrmin_arr[modelI];
      }
      if (yrmax_arr[modelI] > maxyr) {
        maxyr = yrmax_arr[modelI];
      }
    } /* if modelI == 0 */
    ++modelI;
  }  // while (std::getline(ifs, line))

  // std::cout << "Completed reading the file." << std::endl;

  // std::cout << "Checking values..." << std::endl
  //           << "model_index : " << model_index << std::endl
  //           << "model_arr : " << model_arr[model_index - 1] << std::endl;

  nmodel = modelI + 1;
}  // readfile

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
void IGRF::getshc(const std::string& fname, int iflag, int strec,
                  int nmax_of_gh, int gh) {
  std::ifstream ifs{fname};
  if (!ifs) {
    throw std::runtime_error("Cannot open file!");
  }

  std::string line;

  // variables that are used to keep track of the coeffs
  // and their indices
  int ii = 0;
  double g, hh, m, n;

  // extra variables to make file reading work
  double tr, tr1;
  int ln;
  std::string ir;

  // some temporary variables for writing the variables to
  // the storage arrays (i.e. irrelavant variables)
  int line_num;
  double gt, hht, mt, nt;
  double trash, trash1;
  std::string irat;

  ifs.seekg(strec);

  // std::cout << strec << std::endl;

  for (int nn = 1; nn <= nmax_of_gh; ++nn) {
    for (int mm = 0; mm <= nn; ++mm) {
      // std::cout << "in inner loop" << std::endl;
      // while (std::getline(ifs, line)) {
      std::getline(ifs, line);
      // std::cout << line << std::endl;
      std::istringstream ss{line};

      if (iflag == 1) {
        // fgets(inbuff, MAXREAD, stream);
        while (ss >> nt >> mt >> gt >> hht >> trash >> trash1 >> irat >>
               line_num) {
          // std::cout << "inside ss reading" << std::endl;
          // std::cout << nt << mt << gt << hht << trash << trash1 << irat
          //           << line_num << std::endl;
          n = nt;
          m = mt;
          g = gt;
          hh = hht;
          tr = trash;
          tr1 = trash1;
          ir = irat;
          ln = line_num;
        }
        // std::cout << n << ' ' << m << std::endl;
      } else {
        while (ss >> nt >> mt >> trash >> trash1 >> gt >> hht >> irat >>
               line_num) {
          // std::cout << nt << mt << std::endl;
          n = nt;
          m = mt;
          g = gt;
          hh = hht;
          tr = trash;
          tr1 = trash1;
          ir = irat;
          ln = line_num;
        }
        // std::cout << n << ' ' << m << std::endl;
      }
      // std::cout << "Read line" << std::endl;
      if ((nn != n) || (mm != m)) {
        // show some error condition
        // std::cout << "\nError in subroutine getshc" << std::endl;
        throw std::runtime_error("\nnn != n or mm != m");
        break;
      }
      ii += 1;
      switch (gh) {
        case 1:
          gh_first[ii] = g;
          // std::cout << gh_first[ii] << std::endl;
          break;
        case 2:
          gh_second[ii] = g;
          break;
        default:
          // std::cout << "\nError in subroutine getshc" << std::endl;
          throw std::runtime_error("\nError in subroutine getshc");
          break;
      }
      if (m != 0) {
        ii += 1;
        switch (gh) {
          case 1:
            gh_first[ii] = hh;
            // std::cout << gh_first[ii] << std::endl;
            break;
          case 2:
            gh_second[ii] = hh;
            break;
          default:
            // std::cout << "\nError in subroutine getshc" << std::endl;
            throw std::runtime_error("\nError in subroutine getshc");
            break;
        }  // switch (gh)
      }    // if (m!=0)
    }      // for (int mm = 0; mm <= nn; ++mm)
  }        //  for (int nn = 1; nn <= nmax_of_gh; ++nn)

  // std::cout << "Finished importing coefficients." << std::endl;

  // std::cout << "Check results: " << std::endl;
  // for (int l = 0; l < igrf_const::MAXCOEFF; ++l) {
  //   std::cout << gh_first[l] << std::endl;
  // }
  // delete[] irat;
}  // getshc

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
void IGRF::interpsh(double date, double dte1, int nmax1, double dte2, int nmax2,
                    int gh) {
  // int nmax;
  int k, l;
  int ii;
  double factor;

  factor = (date - dte1) / (dte2 - dte1);
  if (nmax1 == nmax2) {
    k = nmax1 * (nmax1 + 2);
    nmax_ = nmax1;
  } else {
    if (nmax1 > nmax2) {
      k = nmax2 * (nmax2 + 2);
      l = nmax1 * (nmax1 + 2);
      switch (gh) {
        case 3:
          for (ii = k + 1; ii <= l; ++ii) {
            gh_arr[ii] = gh_first[ii] + factor * (-gh_first[ii]);
          }
          break;
        case 4:
          for (ii = k + 1; ii <= l; ++ii) {
            ghsv_arr[ii] = gh_first[ii] + factor * (-gh_first[ii]);
          }
          break;
        default:
          std::cout << "\nError in subroutine extrapsh" << std::endl;
          break;
      }
      nmax_ = nmax1;
    } else {
      k = nmax1 * (nmax1 + 2);
      l = nmax2 * (nmax2 + 2);
      switch (gh) {
        case 3:
          for (ii = k + 1; ii <= l; ++ii) {
            gh_arr[ii] = factor * gh_second[ii];
          }
          break;
        case 4:
          for (ii = k + 1; ii <= l; ++ii) {
            ghsv_arr[ii] = factor * gh_second[ii];
          }
          break;
        default:
          std::cout << "\nError in subroutine extrapsh" << std::endl;
          break;
      }
      nmax_ = nmax2;
    }
  }
  switch (gh) {
    case 3:
      for (ii = 1; ii <= k; ++ii) {
        gh_arr[ii] = gh_first[ii] + factor * (gh_second[ii] - gh_first[ii]);
      }
      break;
    case 4:
      for (ii = 1; ii <= k; ++ii) {
        ghsv_arr[ii] = gh_first[ii] + factor * (gh_second[ii] - gh_first[ii]);
      }
      break;
    default:
      std::cout << "\nError in subroutine extrapsh" << std::endl;
      break;
  }

  // std::cout << "Finished interpolating between two models" << std::endl;
  // std::cout << "Check results: " << std::endl;
  // for (int l = 0; l < igrf_const::MAXCOEFF; ++l) {
  //   std::cout << gh_arr[l] << std::endl;
  // }
}

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
void IGRF::extrapsh(double date, double dte1, int nmax1, int nmax2, int gh) {
  // int nmax;
  int k, l;
  int ii;
  double factor;

  factor = date - dte1;
  if (nmax1 == nmax2) {
    k = nmax1 * (nmax1 + 2);
    nmax_ = nmax1;
  } else {
    if (nmax1 > nmax2) {
      k = nmax2 * (nmax2 + 2);
      l = nmax1 * (nmax1 + 2);
      switch (gh) {
        case 3:
          for (ii = k + 1; ii <= l; ++ii) {
            gh_arr[ii] = gh_first[ii];
          }
          break;
        case 4:
          for (ii = k + 1; ii <= l; ++ii) {
            ghsv_arr[ii] = gh_first[ii];
          }
          break;
        default:
          std::cout << "\nError in subroutine extrapsh" << std::endl;
          break;
      }
      nmax_ = nmax1;
    } else {
      k = nmax1 * (nmax1 + 2);
      l = nmax2 * (nmax2 + 2);
      switch (gh) {
        case 3:
          for (ii = k + 1; ii <= l; ++ii) {
            gh_arr[ii] = factor * gh_second[ii];
          }
          break;
        case 4:
          for (ii = k + 1; ii <= l; ++ii) {
            ghsv_arr[ii] = factor * gh_second[ii];
          }
          break;
        default:
          std::cout << "\nError in subroutine extrapsh" << std::endl;
          break;
      }
      nmax_ = nmax2;
    }
  }
  switch (gh) {
    case 3:
      for (ii = 1; ii <= k; ++ii) {
        gh_arr[ii] = gh_first[ii] + factor * gh_second[ii];
      }
      break;
    case 4:
      for (ii = 1; ii <= k; ++ii) {
        ghsv_arr[ii] = gh_first[ii] + factor * gh_second[ii];
      }
      break;
    default:
      std::cout << "\nError in subroutine extrapsh" << std::endl;
      break;
  }
  // std::cout << "Finished extrapolating between two models" << std::endl;
  // std::cout << "Check results: " << std::endl;
  // for (int l = 0; l < igrf_const::MAXCOEFF; ++l) {
  //   std::cout << gh_arr[l] << std::endl;
  // }
}

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
void IGRF::shval3(int igdgc, double flat, double flon, double elev, int gh,
                  int iext, int ext1, int ext2, int ext3) {
  double earths_radius = 6371.2;
  double dtr = 0.01745329;
  double slat;
  double clat;
  double ratio;
  double aa, bb, cc, dd;
  double sd;
  double cd;
  double r;
  double a2;
  double b2;
  double rr;
  double fm, fn;
  double sl[14];
  double cl[14];
  double p[119];
  double q[119];
  int ii, j, k, l, m, n;
  int npq;
  int ios;
  double argument;
  double power;
  a2 = 40680631.59; /* WGS84 */
  b2 = 40408299.98; /* WGS84 */
  ios = 0;
  r = elev;
  argument = flat * dtr;
  slat = sin(argument);
  if ((90.0 - flat) < 0.001) {
    aa = 89.999; /*  300 ft. from North pole  */
  } else {
    if ((90.0 + flat) < 0.001) {
      aa = -89.999; /*  300 ft. from South pole  */
    } else {
      aa = flat;
    }
  }
  argument = aa * dtr;
  clat = cos(argument);
  argument = flon * dtr;
  sl[1] = sin(argument);
  cl[1] = cos(argument);
  switch (gh) {
    case 3:
      bfield_.x = 0;
      bfield_.y = 0;
      bfield_.z = 0;
      break;
    case 4:
      bfield_temp_.x = 0;
      bfield_temp_.y = 0;
      bfield_temp_.z = 0;
      break;
    default:
      std::cout << "\nError in subroutine shval3" << std::endl;
      break;
  }
  sd = 0.0;
  cd = 1.0;
  l = 1;
  n = 0;
  m = 1;
  npq = (nmax_ * (nmax_ + 3)) / 2;
  if (igdgc == 1) {
    aa = a2 * clat * clat;
    bb = b2 * slat * slat;
    cc = aa + bb;
    argument = cc;
    dd = sqrt(argument);
    argument = elev * (elev + 2.0 * dd) + (a2 * aa + b2 * bb) / cc;
    r = sqrt(argument);
    cd = (elev + dd) / r;
    sd = (a2 - b2) / dd * slat * clat / r;
    aa = slat;
    slat = slat * cd - clat * sd;
    clat = clat * cd + aa * sd;
  }
  ratio = earths_radius / r;
  argument = 3.0;
  aa = sqrt(argument);
  p[1] = 2.0 * slat;
  p[2] = 2.0 * clat;
  p[3] = 4.5 * slat * slat - 1.5;
  p[4] = 3.0 * aa * clat * slat;
  q[1] = -clat;
  q[2] = slat;
  q[3] = -3.0 * clat * slat;
  q[4] = aa * (slat * slat - clat * clat);
  for (k = 1; k <= npq; ++k) {
    if (n < m) {
      m = 0;
      n = n + 1;
      argument = ratio;
      power = n + 2;
      rr = pow(argument, power);
      fn = n;
    }
    fm = m;
    if (k >= 5) {
      if (m == n) {
        argument = (1.0 - 0.5 / fm);
        aa = sqrt(argument);
        j = k - n - 1;
        p[k] = (1.0 + 1.0 / fm) * aa * clat * p[j];
        q[k] = aa * (clat * q[j] + slat / fm * p[j]);
        sl[m] = sl[m - 1] * cl[1] + cl[m - 1] * sl[1];
        cl[m] = cl[m - 1] * cl[1] - sl[m - 1] * sl[1];
      } else {
        argument = fn * fn - fm * fm;
        aa = sqrt(argument);
        argument = ((fn - 1.0) * (fn - 1.0)) - (fm * fm);
        bb = sqrt(argument) / aa;
        cc = (2.0 * fn - 1.0) / aa;
        ii = k - n;
        j = k - 2 * n + 1;
        p[k] = (fn + 1.0) * (cc * slat / fn * p[ii] - bb / (fn - 1.0) * p[j]);
        q[k] = cc * (slat * q[ii] - clat / fn * p[ii]) - bb * q[j];
      }
    }
    switch (gh) {
      case 3:
        aa = rr * gh_arr[l];
        break;
      case 4:
        aa = rr * ghsv_arr[l];
        break;
      default:
        std::cout << "\nError in subroutine shval3" << std::endl;
        break;
    }
    if (m == 0) {
      switch (gh) {
        case 3:
          bfield_.x += aa * q[k];
          bfield_.z -= aa * p[k];
          break;
        case 4:
          bfield_temp_.x += aa * q[k];
          bfield_temp_.z -= aa * p[k];
          break;
        default:
          std::cout << "\nError in subroutine shval3" << std::endl;
          break;
      }
      l = l + 1;
    } else {
      switch (gh) {
        case 3:
          bb = rr * gh_arr[l + 1];
          cc = aa * cl[m] + bb * sl[m];
          bfield_.x += cc * q[k];
          bfield_.z -= cc * p[k];
          if (clat > 0) {
            bfield_.y +=
                (aa * sl[m] - bb * cl[m]) * fm * p[k] / ((fn + 1.0) * clat);
          } else {
            bfield_.y += (aa * sl[m] - bb * cl[m]) * q[k] * slat;
          }
          l = l + 2;
          break;
        case 4:
          bb = rr * ghsv_arr[l + 1];
          cc = aa * cl[m] + bb * sl[m];
          bfield_temp_.x += cc * q[k];
          bfield_temp_.z -= cc * p[k];
          if (clat > 0) {
            bfield_temp_.y +=
                (aa * sl[m] - bb * cl[m]) * fm * p[k] / ((fn + 1.0) * clat);
          } else {
            bfield_temp_.y += (aa * sl[m] - bb * cl[m]) * q[k] * slat;
          }
          l = l + 2;
          break;
        default:
          std::cout << "\nError in subroutine shval3" << std::endl;
          break;
      }
    }
    m = m + 1;
  }
  if (iext != 0) {
    aa = ext2 * cl[1] + ext3 * sl[1];
    switch (gh) {
      case 3:
        bfield_.x -= ext1 * clat + aa * slat;
        bfield_.y += ext2 * sl[1] - ext3 * cl[1];
        bfield_.z += ext1 * slat + aa * clat;
        break;
      case 4:
        bfield_temp_.x -= ext1 * clat + aa * slat;
        bfield_temp_.y += ext2 * sl[1] - ext3 * cl[1];
        bfield_temp_.z += ext1 * slat + aa * clat;
        break;
      default:
        std::cout << "\nError in subroutine shval3" << std::endl;
        break;
    }
  }
  switch (gh) {
    case 3:
      aa = bfield_.x;
      bfield_.x = bfield_.x * cd + bfield_.z * sd;
      bfield_.z = bfield_.z * cd - aa * sd;
      break;
    case 4:
      aa = bfield_temp_.x;
      bfield_temp_.x = bfield_temp_.x * cd + bfield_temp_.z * sd;
      bfield_temp_.z = bfield_temp_.z * cd - aa * sd;
      break;
    default:
      std::cout << "\nError in subroutine shval3" << std::endl;
      break;
  }
}

std::array<double, 3> IGRF::values(const double& r, const double& theta,
                                   const double& phi) {
  // convert to lat, long, since shval3 is in this format as well
  // latitude is the north latitude in degrees (+ is north)
  // longitude is east latitude in degrees (+ is east)
  // elev is radius from Earth's center
  double elev = r * (1e-3);
  double lat = 90. - (theta * constants::RAD_TO_DEG);
  double lng = phi * constants::RAD_TO_DEG;

  // std::cout << lat << ' ' << lng << ' ' << elev << std::endl;

  int gh = 3;  // gh=3 is for the main field coefficents (4 for sv)

  // get the x, y, z components of the bfield
  // igdgc = 2 for geocentric, which is what we want
  shval3(2, lat, lng, elev, gh);

  // transform bfield results to spherical coordinates
  std::array<double, 3> values;

  values[0] = sqrt((bfield_.x * bfield_.x) + (bfield_.y * bfield_.y) +
                   (bfield_.z * bfield_.z));
  values[1] = acos(bfield_.z / values[0]);
  values[2] = atan2(bfield_.y, bfield_.x);

  // std::cout << "Finished evaluating bfield values." << std::endl;

  // for (auto val : values) {
  //   std::cout << val << std::endl;
  // }

  return values;
}

// get the declination, inclination, horizontal intensity, and
// total intensity from field values X, Y, Z
// std::array<double, 4> dihf(int gh);
// {
//   int ios;
//   int j;
//   double sn;
//   double h2;
//   double hpx;
//   double argument, argument2;

//   ios = gh;
//   sn = 0.0001;

//   switch (gh) {
//     case 3:
//       for (j = 1; j <= 1; ++j) {
//         h2 = x * x + y * y;
//         argument = h2;
//         h = sqrt(argument); /* calculate horizontal intensity */
//         argument = h2 + z * z;
//         f = sqrt(argument); /* calculate total intensity */
//         if (f < sn) {
//           d = NaN; /* If d and i cannot be determined, */
//           i = NaN; /*       set equal to NaN         */
//         } else {
//           argument = z;
//           argument2 = h;
//           i = atan2(argument, argument2);
//           if (h < sn) {
//             d = NaN;
//           } else {
//             hpx = h + x;
//             if (hpx < sn) {
//               d = PI;
//             } else {
//               argument = y;
//               argument2 = hpx;
//               d = 2.0 * atan2(argument, argument2);
//             }
//           }
//         }
//       }
//       break;
//     case 4:
//       for (j = 1; j <= 1; ++j) {
//         h2 = xtemp * xtemp + ytemp * ytemp;
//         argument = h2;
//         htemp = sqrt(argument);
//         argument = h2 + ztemp * ztemp;
//         ftemp = sqrt(argument);
//         if (ftemp < sn) {
//           dtemp = NaN; /* If d and i cannot be determined, */
//           itemp = NaN; /*       set equal to 999.0         */
//         } else {
//           argument = ztemp;
//           argument2 = htemp;
//           itemp = atan2(argument, argument2);
//           if (htemp < sn) {
//             dtemp = NaN;
//           } else {
//             hpx = htemp + xtemp;
//             if (hpx < sn) {
//               dtemp = PI;
//             } else {
//               argument = ytemp;
//               argument2 = hpx;
//               dtemp = 2.0 * atan2(argument, argument2);
//             }
//           }
//         }
//       }
//       break;
//     default:
//       std::cout << "\nError in subroutine dihf" << std::endl;
//       break;
//   }
//   return (ios);
// }