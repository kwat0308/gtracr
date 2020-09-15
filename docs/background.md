# Background

This code is a sub-project of a much larger project. The original intentions of this project is to aid the evaluation of the three-dimensional atmospheric neutrino flux with high performance.

The evaluation of the atmospheric neutrino flux depends on various factors. The first factor is dependent on the probability of muon production $P_{\mu, prod}$ within some wedge in the Earth's atmosphere. Of course, this depends on the number of particle interactions that occur within this "packet" of atmosphere depending on the energy of such particles. This will also influence the direction of propagation and momentum of the muon.

Such values then affect the muon trajectory within the geomagnetic field. The probability of muon decay $P_{\mu, dec}$ is then evaluated at each step based on the momentum of the muon. The distribution of the momentum magnitude along each direction is a key component, as this will influence the direction in which the neutrinos are produced (i.e. the probability of neutrino production in a certain solid angle relative to the muon $P_{\nu, prod}$) in the muon decay. Lastly, we need to calculate the probability that the detector will detect $P_{\nu, prod}$ number of neutrinos at a certain solid angle $P_{detect}$.

As such, there are a total of 4 probability calculations required, and the atmospheric neutrino flux can be evaluated by the product of each of these probabilities (actually, even before all of these, we need to evaluate the energies in which cosmic rays can come towards Earth's atmosphere!). Additionally, this is done with $N_{part}$ many particles **simultaneously!** Essentially, there are a lot of calculations to perform.

As one can imagine, the evaluation of the 3-D atmospheric neutrino flux has been a real computational challenge that many had a hard time to tackle. The current calculations, made by Honda, uses a Monte-Carlo approach in which he injects $N_{part}$ number of cosmic ray trajectories at random locations within an injection sphere in a Monte-Carlo-like manner. He then performs the calculations above for each cosmic ray. This works (as with most Monte-Carlo codes), but it is terribly slow, taking almost two weeks (as with most Monte-Carlo codes)!

Today, we have a lot of new and modern computation methods available to us with new technologies emerging. As such, we want to use such technologies to upgrade the evaluation of the 3-D atmospheric neutrino flux using multi-threading with GPU cores. GPUs, or Graphical Processing Units, contain more than a thousand cores per unit. Since such trajectory evaluations are relatively simple, these cores can handle such evaluation processes. This means that we can perform such trajectory calculations _simultaneously_ in the time of the calculation of one trajectory!

As such, our goal is to create such a high-performance 3-D atmospheric neutrino flux calculator!

This project will serve three key purposes for the atmospheric neutrino calculation:

- Evaluate the geomagnetic cutoff rigidities at a certain location on Earth (for example at IceCube).
- Evaluate the trajectories for muons, including the probability of muon decay throughout its propagation.
- Create an experimental platform for [CORSIKA 8](https://arxiv.org/pdf/1902.02822.pdf) with GPU acceleration (multi-threading with GPU cores).

We hope to achieve these results in our package, and with the code being open-source, anyone with some GPU can run this very easily (we also hope to integrate CPU parallelization for those without GPUs)!
