## VENUS 1986-1997

    The project starts with the creation of the VENUS model by Klaus WERNER,
    during a postdoc stay at the Brookhaven National Laboratory between 1986
    and 1989. This was one of the first codes able to simulate relativistic
    proton-proton and heavy-ion collisions, the latter ones being realized
    at the time at the SPS accelerator at CERN. VENUS was already using the
    concept of parallel scatterings in the Gribov-Regge framework, but the
    scatterings were considered to be completely "soft", sufficient for the
    moderately high SPS energies.

## NEXUS 1997-2004

    In the early 2000s, the RHIC collider in Brookhaven started operation,
    with Au+Au collisions at an CMS energy of 200 AGeV. This made it
    necessary (in the late 90s) to reconsider the VENUS project in order to
    implement "hard processes". It was decided to fuse the QGSJET code of
    Sergej OSTAPCHENKO and the VENUS code of Klaus WERNER, with the help of
    several PhD students (Hajo DRESCHER, Michael HLADIK, Tanguy PIEROG), to
    have a consistent description of soft and hard scatterings, and
    therefore an up-to-date model for RHIC simuations, named NEXUS.

## EPOS1 2004-2012

    In the NEXUS approach, it was clear from the beginning that a simple
    approach with just elementary Pomerons (I-diagrams) will lead to
    contradictions, and it was decided to explicitely implement as well
    Y-diagrams, which amounts to splitting a parton ladder into two. But
    then one also needs to allow for splitting the two legs again, and so
    on. Unfortunately, such a cascade of splittings is impossible to do in a
    framework with rigorous energy conservation. Two possible solution
    emerged, and were further developed around 2006 as EPOS (keeping full
    energy consevation, but treating Y diagrams in an effective way), and as
    QGSJETII (giving up energy conservation, but summing up the splittings
    to infinit order).

    In EPOS, in addition to the multiple Pomeron exchanges, referred to as
    primary scatterings happening at t=0, secondary scatterings were added,
    treating the reinteractions of the particles produced initially. In
    2007, the core-corona picture was introduced (Klaus WERNER), which
    amounts to identify a core part, which then evolves macroscopically as a
    fluid. Those particles which do not participate to the core are called
    corona particles. In EPOS1, the core expansion was simply parameterized,
    with a flow put in by hand. EPOS1 turned out to b quite successful to
    descibe the first LHC results (with the code being contructed and tuned
    before LHC). Around 2012, EPOS LHC was created (Tanguy PIEROG, K. Werner), 
    which is essentially the last EPOS1 version (1.99) with some finetuning 
    compared to LHC results from run 1.{" "}

## EPOS2 2010-2013

    In EPOS2 (I. KARPENKO, T. PIEROG, and K. WERNER), an efficient code for 
    solving the hydrodynamic equations in 3+1 dimensions was implemented,
    inluding the conservation of baryon number, strangeness, and electric
    charge, employing a realistic equation-of-state, compatible with lattice
    gauge results, using a complete hadron resonance table, making our
    calculations compatible with the results from statistical models.
    Nonlinear effects (parton ladder fusions) were treated by simply adding
    by hand a term x^epsilon to the Pomeron amplitudes, with x being an
    energy fraction (epsilon method).

## EPOS3 2013-2017

    In EPOS3 (B. GUIOT, I. KARPENKO, T. PIEROG, and K. WERNER), a 3D+1 viscous
    hydrodynamical evolution was implemented, starting from flux tube initial 
    conditions, being generated in the Gribov-Regge multiple scattering framework.
    An individual scattering is referred to as Pomeron, identified with a
    parton ladder, eventually showing up as flux tubes (or strings). Each
    parton ladder is composed of a pQCD hard process, plus initial and final
    state linear parton emission. Klaus WERNER introduced a new way to deal
    with nonlinear effects, namely by associating a saturation scale Q_s to
    each parton ladder (Q_s method). Tanguy PIEROG introduced a smart way to
    combine the epsilon and the Q_s method, by first using epsilon to define
    an effective Pomeron amplitude G_eff, and then making the link to the
    QCD ladder via G_eff = k*G_QCD(Q_s), where the momentum dependence of k
    should be chosen such that factorization is assured (which turned out to
    be difficult). Finally, Benjamin GUIOT made major contributions
    concerning the implementation of heavy flavor.

