# The Earth's magnetic field

## Introduction

- The Earth's magnetic field is the primary influencer of the behaviour of the particle trajectory
- They can be defined in many different ways:

  - By using geocentric spherical coordinates
  - By using L-shells
  - By using geomagnetic latitude and longitude
  - Most often, we will use spherical coordinates as it is most easily expressed in that form and thus easiest to integrate it with

- All particle trajectories are assumed to be in the source-free region, that is, on Earth's surface and above, which allows us to express the magnetic field as the negative gradient of the potential $\vec{B} = -\nabla V$.

## Approximations

- Ideal Dipole Approximation

  - Here, we approximate the Earth's magnetic field as an ideal dipole.

  - This can be done based on the fact that we assume that the Earth is:

    - A perfect sphere
    - Spherically symmetric

      - This allows azimuthal symmetry to occur with the magnetic field, so that the $\phi$ component is zero.

  - The actual equations are as follows:

    - $B_r = 2\left(\dfrac{R_E}{r}\right)^3 g_1^0\cos\theta$
    - $B_\theta = \left(\dfrac{R_E}{r}\right)^3 g_1^0\sin\theta$
    - $B_\phi = 0$
    - $g_1^0 = -29404.8\times10^{-9}$ is the first parameter value from the IGRF model and the mean value of the magnetic field at the equator
    - $R_E = 6371.2\text{km}$ is the radius of Earth

- The true IGRF model

  - Based on the IGRF 13th generation model announced in 2020 (<https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html>)

  - We can alternatively use the DGRF model, which is the definite model, but this is from 2015 and I think we can trust the IGRF model as many papers use this instead of the DGRF model

  - The order of truncation is $N=13$ (as of 2020).

  - The potential is expressed as the series expansion in spherical coordinates:

    - $V(r, \theta, \phi, t) = R_E \sum\limits_{n=1}^{N}\sum\limits_{m=0}^{n}\left(\dfrac{R_E}{r}\right)^{n+1}\left[g_n^m(t)\cos(m\phi) + h_n^m(t)\sin(m\phi)\right]P_n^m(\cos\theta)$

  - To get the magnetic field components, we take the gradient in spherical coordinates. This yields us with:

    - $B_r = R_E \sum\limits_{n=1}^{N}\sum\limits_{m=0}^{n}(n+1)\left(\dfrac{R_E^{n+1}}{r^{n+2}}\right)\left[g_n^m(t)\cos(m\phi) + h_n^m(t)\sin(m\phi)\right]P_n^m(\cos\theta)$

    - $B_\theta = \dfrac{R_E}{r} \sum\limits_{n=1}^{N}\sum\limits_{m=0}^{n}\left(\dfrac{R_E}{r}\right)^{n+1}\left[g_n^m(t)\cos(m\phi) + h_n^m(t)\sin(m\phi)\right]P_n^{m'}(\cos\theta)(\sin\theta)$

    - $B_\phi = \dfrac{R_E}{r\sin\theta} \sum\limits_{n=1}^{N}\sum\limits_{m=0}^{n}\left(\dfrac{R_E}{r}\right)^{n+1}\left[m g_n^m(t)\sin(m\phi) - mh_n^m(t)\cos(m\phi)\right]P_n^m(\cos\theta)$

  - The Gauss coeffients $g_n^m(t), h_n^m(t)$ have a secular variation in a linear fashion over a 5-year interval:

    - $g_n^m(t) = g_n^m(T_0) + g_n^{m'}(T_0)(t-T0)$
    - $h_n^m(t) = h_n^m(T_0) + h_n^{m'}(T_0)(t-T0)$
    - $T_0$ is the epoch preceding $t$, which is an exact multiple of 5 years.
    - The first time derivatives can be obtained by the secular variation coeffients (for 2015-2020 models only).
    - As we are only dealing with second-like time intervals, we can approximate this time variation of the coefficients as constant, which are provided as the IGRF coefficients.

  - The function $P_n^m(x)$ are the associated Legendre polynomials, but are Schmidt quasi-normalized (source (IGRF 12th generation paper): <https://earth-planets-space.springeropen.com/articles/10.1186/s40623-015-0228-9>)

    - The form of the Schmidt quasi-normalized associated Legendre polynomials are used primarily in geomagnetic calculations and can be found in this link: <https://www.spenvis.eu/help/background/magfield/legendre.html#Schmidt1>.
    - Implementing this is very simple, just a problem of adding a normalization term in front of a certain case (when m > 0).

## Implementation

We will perform the following to implement the IGRF magnetic field:

- import the coefficients (which are stored in a csv file) and store them into a matrix-like data structure
- use this matrix for the double sum of the magnetic field coefficients

Some ideas:

- could we make an entire class and have all of this as members / member functions?
- This can be a derived class from the dipole approximation, since the member functions of the dipole approximation can appropriately be used for the IGRF version.
