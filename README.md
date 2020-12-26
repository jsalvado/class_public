CLASS_lrs: Cosmic Linear Anisotropy Solving System with Long Range Interactions
==============================================

Forked by Ivan Esteban and Jordi Salvado from the CLASS code by Julien
Lesgourgues and Thomas Tram.

This is a fork of CLASS 2.9.2, designed to implement the cosmology of
scalar-mediated long range interactions among fermions. For details,
see the companion paper.

There are two ipython notebooks in the notebooks/ folder, _longrange.ipynb_
and _longrange\_checks.ipynb_, that you can have a look at to familiarize
yourself with this code.

Input parameters
-----------------------------------
* _longrangescalar_ -- If 'y' or 'Y', include scalar-mediated long range
interactions
* _longrangescalar\_phi\_pt_ -- If 'y' or 'Y', include scalar field
perturbations (see caveats below)
* _longrangescalar\_nuggets_ -- If 'y' or 'Y', include fermion _nugget_
formation (see caveats below)
* _lrs\_g\_over\_M_ -- Interaction coupling over scalar mass in
eV<sup>-1</sup>
* _lrs\_M\_phi_ -- Scalar mass in eV. Only relevant if scalar field
perturbations are being taken into account
* _lrs\_m\_F_ -- Fermion mass in eV
* _log10\_lrs_ -- If 'y' or 'Y', the input parameters _lrs\_g\_over\_M_,
_lrs\_M\_phi_ and _lrs\_m\_F_ are interpreted as the base-10 logarithms
of the corresponding physical parameters
* _lrs\_g\_F_ -- Number of internal degrees of freedom of the fermion
* _lrs\_T\_F_ -- Fermion temperature over photon temperature, assumed to
be constant

Caveats
-----------------------------------
The main physical issue is that density perturbations of a system of
attractive long-range interacting fermions are unstable. When fermions become
non-relativistic, they collapse into non-linear _nuggets_ that a linear
code as CLASS cannot properly model.
This can be overcome by assuming that _nugget_ formation takes place
over timescales and distances much smaller than cosmological scales. In
this case, one can enforce an _instantaneous_ transition to a dust-like
behaviour below some pre-computed fermion temperature T<sub>nug</sub>. This is
equivalent to setting

_longrangescalar\_phi\_pt_ = 'N'

_longrangescalar\_nuggets_ = 'y'


It is important to turn off scalar field perturbations, as they may make
perturbations exponentially grow if T<sub>nug</sub> is not 100% precise. For the
details on how we estimate T<sub>nug</sub>, see Appendix B in the companion paper.

The code can also be run with scalar field perturbations, in which case
non-relativistic fermion density perturbations will exponentially grow.
Notice, though, that this growth is stronger at scales much smaller than
those computed by CLASS, what in turn affects the background. Such
simulations are therefore unphysical.

Apart from that, the implementation largely follows the CLASS 2.9.2
implementation of warm relics (named _ncdm_ in the code). The fermion momentum 
distribution function is hard-coded to be Fermi-Dirac of an ultrarelativistic relic.
This can be easily modified in the function _background\_lrs\_distribution_
of the file _source/longrange.c_. Fluid approximations are in principle implemented, but 
their accuracy has not been tested. Tensor perturbations have not been properly implemented.
