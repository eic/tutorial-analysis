---
title: "Full Chain Analysis"
teaching: 15
exercises: 10
questions:
- "How do I bring all of this together if I'm starting from scratch?"
objectives:
- "Become familiar with the full analysis chain"
keypoints:
- "There are a few steps to go through before we get to the file we analysed previously."
- "Good for testing, but use simulation campaign output where possible!."
---

In this short session, we'll go through a brief run through of how we actually ended up with a file like the one we ran our script on before. There are X basic steps which we'll look at individually, and then combine together:

1. Generate an input file (typically hepmc, other formats are useable). This is usually from some external event generator.
2. Afterburn the file and apply beam effects (might be skipped in some cases).
3. Process the input through the simulation, DD4HEP.
4. Reconstruct the DD4HEP output with EICrecon.
5. Analyse the EICrecon output with analysis script.

Note that for low level analyses, you could also directly analyse the DD4HEP output from step 3.

## Event Generator Input Files

I won't say too much on this since this strongly depends upon the channel you want to simulate and analyse. The files here will likely come from an external event generator, for example -

- PYTHIA
- BeAGLE
- DJANGOH
- MILOU
- eSTARlight
- LaGER
- DEMPgen
- elSPectro

... and may others. However, regardless of what you use, the output is likely some form of .hepmc file with event by event particle/vertex info. For example -

> Example HEPMC Event:
> An example event from a HEPMC file is shown below. In this example event, we have an input 5 GeV electron on a 41 GeV proton. We have one vertex and three outgoing particles, a scattered electron, a pion, and a neutron. In our header, we also have an event weight included.
> 
> E 1 1 5
> U GEV MM
> A 0 weight 4.813926604168258e-07
> P 1 0 11 6.123233963758798e-16 0.000000000000000e+00 -4.999999973888007e+00 5.000000000000000e+00 5.109989488070365e-04 4
> P 2 0 2212 -0.000000000000000e+00 -0.000000000000000e+00 4.100000000000000e+01 4.101073462535657e+01 9.382720881600054e-01 4
> V -1 0 [1,2]
> P 3 -1 11 -6.872312444834133e-01 1.924351128807063e+00 -4.281657822517654e+00 4.744260534644128e+00 5.109989488070365e-04 1
> P 4 -1 211 1.042011265882083e+00 -1.600831989262599e+00 1.404460452649878e+00 2.374960954263115e+00 1.395701800000037e-01 1
> P 5 -1 2112 -3.547800213986697e-01 -3.235191395444645e-01 3.887719739597977e+01 3.889151313644933e+01 9.395654204998098e-01 1
{: .callout}

Typically, we also need to incorporate beam effects. This is done via the use of the afterburner.

## Applying Beam Effects - Afterburner

Afterburner applies beam effects to an existing hepmc file. These include effects due to the crabbing of the beam bunches and the crossing angle. Afterburner is pre-installed in eic-shell. We can run it via -

```console
abconv
```

However, we'll need an input file to do anything, we can also check other options quickly with -

```console
abconv -h
```

Note that when we run Afterburner, it will try to pick up the input beam energies and apply the relevant configuration. We can force a different configuration if we want (see the options from the help printout). We could for example though run -

```console
abconv $File -o $OutputFilename
```

where $File is our input hepmc file from our generator, and $OutputFilename is whatever we want our output to be called.

Regardless of whether we want or need to afterburn the file, we can feed in our hepmc file to DD4HEP and process our events through the detector simulation.

## Simulation



## Reconstruction

## Reconstruction Output Analysis

## Combining Everything

## Farming