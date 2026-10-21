![logo_small](logos/logo_small.png)

This repository contains all the core code for the artificial protomodeller, https://arxiv.org/abs/2012.12246. For results and discussions, see also https://smodels.github.io/protomodels/. This work was also presented to the reinterpretation forum, https://cds.cern.ch/record/2752548.

DISCLAIMER: This is not a finished software product. The code is public, but not published. Your mileage may vary.

## walker

Contains the code for the random walker. The hiscore class is also here.

## builder

Contains the protomodel class -- the data object -- and the manipulator, i.e. the algorithm object.

## tester

All the code around likelihoods, combining likelihoods, statistics.

## user-provided priors

The prior that is computed in `Combiner.computePrior` can be modified by the
user. To do so, write a python file that contains a function `myprior`, which
takes the protomodel and the nll flag and returns a float, and tell the
walker about it:

    slurm_walk.py --userprior /path/to/myprior.py ...

For `nll=True` the returned value is added to the (negative log likelihood
of the) prior, for `nll=False` it is multiplied to the prior.

The user prior is part of the run environment (i.e. of `run.dict`, see
`pbase/runEnviron.py`), where it shows up as the `userprior` entry. So, if
you run a walker yourself, use

    environ = RunEnviron.create ( userprior = "/path/to/myprior.py" )
    walker = RandomWalker ( environ, rvars )

(or set `environ.userprior` directly). See `snippets/exampleUserPrior.py` for
a simple example, and `snippets/userPriorCharginos.py` for a more
physics-inspired one (it only accepts models with an off-shell chargino, i.e.
no on-shell W, Z or h decays of 1000023 or 1000024).

## ptools

Various tools and helpers needed for protomodelling, e.g. history recording, or the code that creates "fake" databases.

## plotting

Plotting scripts

## notebooks

Jupyter notebooks

## unittests

Unit testing.
