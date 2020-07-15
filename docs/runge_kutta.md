# Runge-Kutta Integration

## Initial implementation

We used the differential equations obtained by D.F. Smart from the following source: <https://www.dartmouth.edu/~sshepherd/research/Shielding/docs/Smart_00.pdf>:

- $\dfrac{\mathrm{d} v_{r}}{\mathrm{d} t}=\dfrac{e\left(v_{\theta} B_{\phi}-v_{\phi} B_{\theta}\right)}{m c}+\dfrac{v_{\theta}^{2}}{r}+\dfrac{v_{\phi}^{2}}{r}$
- $\dfrac{\mathrm{d} v_{\theta}}{\mathrm{d} t}=\dfrac{e\left(v_{\phi} B_{r}-v_{r} B_{\phi}\right)}{m c}-\dfrac{v_{r} v_{\theta}}{r}+\dfrac{v_{\phi}^{2}}{r \tan \theta}$
- $\dfrac{\mathrm{d} v_{\phi}}{\mathrm{d} t}=\dfrac{e\left(v_{r} B_{\theta}-v_{\theta} B_{r}\right)}{m c}-\dfrac{v_{r} v_{\phi}}{r}-\dfrac{v_{\theta} v_{\phi}}{r \tan \theta}$

The coordinate components are evaluated as such:

- $\dfrac{dr}{dt} = v_r$
- $\dfrac{d\theta}{dt}=\dfrac{v_\theta}{r}$
- $\dfrac{d\phi}{dt} = \dfrac{v_\phi}{r\sin\theta}$

## Issues

- I was checking the radius values during the RK integration to see why no trajectories were going back to Earth, and I figured out that the integration process kept going even after the stop altitude. This was due to the fact that I thought that the stepsize for the integration process was determined by the formula $h = \dfrac{b-a}{n}$ where a,b are the starting / stopping altitude and n is the maximum number of iterations. It seems that this only applies to Euler's method / numerical integration methods based on Riemann sums.

  - As such, I will make the stepsize be set by the user, and determine the maximum number of steps from this. I think the formula will go like $n = \dfrac{b-a}{h}$, but this may not be the case as this is essentially identical to the formula for step size.

    - ~~This now seems to work. We will use this instead from now on.~~ This doesnt actually work, as doing this will make the array length a lot longer than it actually would be. The way I am doing this now is by making user inputs for both stepsize and maximum number of iterations, and trim the zero values at the end that were not used

      - there should be a way to determine the condition that max number of iterations > max array length in which contents will be appended

- My friend had told me that the extra terms that pop up at the end of the Lorenz force term are due to unit vector conversions from Cartesian to Spherical coordinates. Now that I know this fact, I decided to test this out by hand to see if I do get this equation (the Equation obtained from D.F.Smart and M.A.Shea, Advances in Spce Research, 2004).

  - It seems like they somehow did not consider the velocity term within the Lorenz factor $\gamma$ since their expression only contains a relativistic term in the Lorenz force component. This is a mystery to me, as when I evaluated this myself I got something much more complicated than his version:
  - $\eta = \dfrac{q}{m\gamma c^2}$
  - $\dfrac{dv_r}{dt} = rv_\theta^2 + rv_\phi^2\sin^2\theta + \eta\left[(v_\theta B_\phi - B_\theta v_\phi)(c^2-v_r^2) - rv_rv_\theta (v_rB_\phi - v_\phi B_r) + rv_\theta v_\phi \sin\theta (v_rB_\theta - v_\theta B_r)\right]$
  - $\dfrac{dv_\theta}{dt} = \dfrac{1}{r}\left(- 2v_rv_\theta + rv_\phi^2\sin\theta\cos\theta + \eta\left[rv_rv_\theta(v_\theta B_\phi - B_\theta v_\phi) - (c^2 - (rv_\theta)^2) (v_rB_\phi - v_\phi B_r) + r^2v_\theta v_\phi \sin\theta (v_rB_\theta - v_\theta B_r)\right]\right)$
  - $\dfrac{dv_\phi}{dt} = \dfrac{1}{r\sin\theta}\left(- 2v_rv_\phi\sin\theta -2rv_\theta v_\phi \cos\theta + \eta\left[rv_rv_\phi\sin\theta(v_\theta B_\phi - B_\theta v_\phi) - r^2v_\phi v_\theta\sin\theta (v_rB_\phi - v_\phi B_r) + (c^2 - (rv_\phi\sin\theta)^2) (v_rB_\theta - v_\theta B_r)\right]\right)$

  - This was evaluated using the relativistic acceleration equation from the relativistic version of Newton's 2nd law: $\vec{a} = \dfrac{1}{m\gamma}\left(\vec{F} - \dfrac{(\vec{v}\cdot\vec{F})\vec{v}}{c^2}\right)$ (obtained from wikipedia)

    - not really sure how this is properly derived though...

  - Testing the equations above, this was shown to not be correct, as I get the same trajectory for any particle energy. Switching back to the original equations, we get the trajectory that is identical to Baldini's paper, which is a good relief.

### Revisiting my own implementation of the equations

I realized that the equation for acceleration above $\left(\vec{a} = \dfrac{1}{m\gamma}\left(\vec{F} - \dfrac{(\vec{v}\cdot\vec{F})\vec{v}}{c^2}\right)\right)$ has a dot product with the same variable $\vec{v}$. This means that the second term within the brackets vanishes, leaving us with: $\vec{a} = \dfrac{\vec{F}}{m\gamma} = \dfrac{q}{m\gamma}\left(\vec{v}\times\vec{B}\right)$. Using the expression for acceleration in spherical coordinates, we obtain a much simpler set of ODEs:

- $\dfrac{dv_r}{dt} = \dfrac{q}{m\gamma}(v_\theta B_\phi - B_\theta v_\phi) + rv_\theta^2 + rv_\phi^2\sin^2\theta$
- $\dfrac{dv_\theta}{dt} = \dfrac{q}{m\gamma r}(v_rB_\phi - v_\phi B_r) - \dfrac{2v_r v_\theta}{r}+ v_\phi^2\sin\theta\cos\theta$
- $\dfrac{dv_\phi}{dt} = \dfrac{q}{m\gamma r\sin\theta}(v_rB_\theta - v_\theta B_r) - \dfrac{2v_r v_\phi}{r} - \dfrac{2v_\theta v_\phi}{\tan\theta}$

- These set of ODEs are derived from the relativistic acceleration formula in terms of the relativistic force:

  - $\vec{a} = \dfrac{\vec{F}}{m\gamma} = \dfrac{q}{m\gamma}\left(\vec{v}\times\vec{B}\right)$
  - Simplifications were made from the actual equation: $\vec{a} = \dfrac{1}{m\gamma}\left(\vec{F} - \dfrac{(\vec{v}\cdot\vec{F})\vec{v}}{c^2}\right)$ since $\vec{v}\cdot\vec{F} = q \vec{v}\cdot \left(\vec{v}\times\vec{B}\right) = 0$

  - The extra terms are obtained from using the expression of acceleration defined in spherical coordinates (from Wikipedia):

    - $\vec{a} = \left(\dfrac{dv_r}{dt} - rv_\theta^2 - rv_\phi^2\sin^2\theta\right)\hat{r}$ $+ \left(r\dfrac{dv_\theta}{dt} +2v_r v_\theta - v_\phi^2\sin\theta\cos\theta\right)\hat{\theta}$ $+ \left(r\sin\theta\dfrac{dv_\phi}{dt} + 2v_r v_\phi\sin\theta + 2rv_\theta v_\phi\cos\theta\right)\hat{\phi}$

~~These equations were cross-checked with Baldini's equations, and the trajectories were consistent. As such, I will be using these equations above from now on.~~

### Revisiting the Runge Kutta implementation once again

We found out that the equations obtained by D.F. Smart were indeed correct, and were the same equations as which we derived. The only difference is that:

- In _our equation_ we assumed that the theta and phi component of the velocity are in _radians / second_
- In _their equation_ they assumed that those components are in _meters / second_ by letting:

  - vtheta = r * vtheta
  - vphi = r _sin(theta)_ vphi In this sense, **their version is more correct.** Why? Since each component of the velocity / momentum should have **some sort** of units.

As of such, we now use Smart's version of the coupled ODEs, which provide more physically sound reasoning behind the equations.

### Using momentum instead of velocity

To allow natural units to be played around safely within the code, we will be using momentum instead of velocity. This will allow the equations to make more sense in a conventional view. This requires a few conversions between momentum and velocity via the relation $p = \gamma m v$.

As such, the final equations are as followed:

- $\dfrac{\mathrm{d} p_{r}}{\mathrm{d} t}=\dfrac{1}{m\gamma}\left(e\left(p_{\theta} B_{\phi}-p_{\phi} B_{\theta}\right)+\dfrac{p_{\theta}^{2}}{r}+\dfrac{p_{\phi}^{2}}{r}\right)$
- $\dfrac{\mathrm{d} p_{\theta}}{\mathrm{d} t}=\dfrac{1}{m\gamma}\left(e\left(p_{\phi} B_{r}-p_{r} B_{\phi}\right)-\dfrac{p_{r} p_{\theta}}{r}+\dfrac{p_{\phi}^{2}}{r \tan \theta}\right)$
- $\dfrac{\mathrm{d} p_{\phi}}{\mathrm{d} t}=\dfrac{1}{m\gamma}\left(e\left(p_{r} B_{\theta}-p_{\theta} B_{r}\right)-\dfrac{p_{r} p_{\phi}}{r}-\dfrac{p_{\theta} p_{\phi}}{r \tan \theta}\right)$
- $\dfrac{dr}{dt} = \dfrac{1}{m\gamma}\left(p_r\right)$
- $\dfrac{d\theta}{dt}=\dfrac{1}{m\gamma}\left(\dfrac{p_\theta}{r}\right)$
- $\dfrac{d\phi}{dt} = \dfrac{1}{m\gamma}\left(\dfrac{p_\phi}{r\sin\theta}\right)$
