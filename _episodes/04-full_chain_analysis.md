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

In this short session, we'll go through a brief run through of how we actually ended up with a file like the one we ran our script on before. There are 5 basic steps which we'll look at individually, and then combine together:

1. Generate an input file (typically hepmc, other formats are useable). This is usually from some external event generator.
2. Afterburn the file and apply beam effects (might be skipped in some cases).
3. Process the input through the simulation, DD4HEP.
4. Reconstruct the DD4HEP output with EICrecon.
5. Analyse the EICrecon output with analysis script.

Note that for low level analyses, you could also directly analyse the DD4HEP output from step 3. You may also wish to consult [Holly's slides from the April 2024 software meeting for an overview](https://indico.cern.ch/event/1343984/contributions/5927492/attachments/2843633/4971409/tutorial_overview.pdf) of each of these steps and how they fit into this production chain.

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

To process our events through the simulation, we need to get the detector geometry. The simplest way is simply to source the nightly build within eic-shell -

```console
./eic-shell
source /opt/detector/epic-main/bin/thisepic.sh
```

We can check this worked as intended by checking that the DETECTOR_PATH variable is now set. Do so via -

```console
ls $DETECTOR_PATH
```

If we do this without sourcing thisepic.sh, we should get an error. Now, we should see a range of .xml files (which outline various detector configurations). We could also compile our own version of the detector within eic-shell. You might want to do so if you are actively iterating on the design of specific detector for example. See the [GitHub page](https://github.com/eic/epic) for instructions.

We can now process a simulation. Be aware that this may take some time, so to test it, try processing a small number of events first. Check the options we can provide via -

```console
npsim -h
```

A typical simulation command might look something like -

```console
npsim --compactFile $DETECTOR_PATH/epic_craterlake.xml --numberOfEvents 1000 --inputFiles input.hepmc --outputFile output.edm4hep.root
```

Most of the arguments are pretty self explanatory. As a quick demo, I'll run -

```console
npsim --compactFile $DETECTOR_PATH/epic_craterlake_5x41.xml --numberOfEvents 10 --inputFiles eic_DEMPgen_5on41_ip6_pi+_1B_1.hepmc --outputFile DEMPgen_5on41_pi+_10_TestOutput.edm4hep.root
```

When we run this, we'll get lots of printouts to screen, we can supress this by adding the -v5 argument too (only errors will be printed). Once we have our simulation output, we can now reconstruct our events.

## Reconstruction

We can run eicrecon pretty straightforwardly, within eicshell, try -

```console
eicrecon -h
```

which should again, print out the various options we have available. An example command to run the reconstruction on a file might look like this -

```console
eicrecon -Ppodio:output_file=eicrecon_out.root -Pjana:nevents=1000 -Pdd4hep:xml_files=epic_craterlake.xml sim_output.edm4hep.root
```

Again, this might take a long time. So test a small sample of events first. Following up on my simulation demo, I'll run -

```console
eicrecon -Ppodio:output_file=DEMPgen_5on41_pi+_10_TestReconOutput.edm4hep.root -Pjana:nevents=10 -Pdd4hep:xml_files=epic_craterlake_5x41.xml DEMPgen_5on41_pi+_10_TestOutput.edm4hep.root
```

eicrecon will look for the detector .xml file in $DETECTOR_PATH, so make sure the detector geometry is sourced before running eicrecon.

## Combining Everything

Ok, great. We now have a file we could run our earlier analysis script on. But what if we wanted to do all of this from scratch? Well, the easiest way might be to put all of this in a shell script. So, pulling all of our commands together -

```console
#! /bin/bash

source /opt/detector/epic-main/bin/thisepic.sh
eval npsim --compactFile $DETECTOR_PATH/epic_craterlake.xml --numberOfEvents 1000 --inputFiles input.hepmc --outputFile output.edm4hep.root
sleep 3
eval eicrecon -Ppodio:output_file=eicrecon_out.root -Pjana:nevents=1000 -Pdd4hep:xml_files=epic_craterlake.xml sim_output.edm4hep.root

exit 0
```

I've premade a version of this with the commands I ran earlier, so we can run it and see what happens.

## Farming

Ok, great. We can do (almost) all of the processes we need in one command. But as we've seen, the processing can take a while. Realistically, we're probably going to want to parallelise this in some way. With access to the JLab iFarm or the BNL systems (Condor). We can create and submit compute jobs for this purpose. This is getting a bit beyond the scope of this tutorial, but some things to consider -

- Our job needs to either access the container, or process eic-shell within the job (more on this below)
- The job itself should be as simple as possible, just exectuing a command with some arugments. Our script above is a good candidate (with some work)
- As is, our script is fairly inflexible. We should probably make things like the input and output file names variables that are set based upon arguments we provide.
- We need to consider the resource usage of our job carefully.
- Pathing can be tricky, we need to make sure the farm/compute node picks up the correct paths such as $DETECTOR_PATH (this is a common job error).
- As always, TEST first. Run a small job that runs quickly interactively, THEN submit it is a small compute job. Compare the outputs.

For our first point, one easy (and not recommended for Condor jobs!) way to do this is via an EOF line -

```console
#! /bin/bash

cat <<EOF | eic-shell
./Basic_Bash.sh
EOF

exit 0
```
This just starts eic-shell and runs our earlier script. We can run this WITHOUT running eic-shell first. Note that this is a bit of a cheat to address a pathing issue. $DETECTOR_PATH will be interpreted by the script BEFORE the EOF script so our variable will be mis-set. We can get around this by running a script. Ideally, for our compute job, we should probably also explicitly set our paths to eic-shell and the bash script in some way.

With changes like this made, we could then make a quick job and submit it. This is a bit beyond the tutorials, but for some farm examples, see [this job script](https://github.com/sjdkay/ePIC_PairSpec_Sim/blob/main/Farm_Bash_Scripts/PairSpec_Sim_Job.sh) and this [job submission script](https://github.com/sjdkay/ePIC_PairSpec_Sim/blob/main/Farm_Bash_Scripts/PairSpec_Sim_Job.sh) I use as a template. Feel free to use these as a template for your own jobs, but please thoroughly read through and understand them before submitting a huge number of jobs. Keep the comments above in mind too.

Also a quick disclaimer, my experience in running jobs is limited to systems I know (which does not include the BNL systems). As such, I can't advise much beyond general slurm job style questions on BNL/Condor. I'm also aware that using EOF in scripts was not encouraged in BNL jobs, see [this discussion on mattermost](https://chat.epic-eic.org/main/pl/fo9954siwigyjckasrnd4xufxw) for more.

## Warnings

Finally, a major disclaimer. A lot of the time, you should NOT be starting from scratch and processing through the simulation and reconstruction yourself. There are numerous reasons -

- Computing time intensive
- Versioning errors/mismatch
- Not as reproducible (if you find an error, people will need to try and reproduce it from your environment)

Where possible, use files from official simulation campaigns (bringing us full circle, see the first lesson for using a simulation campaign file in a script!). That being said, for testing and iterating rapidly on a design change, running small jobs yourself may be the way to go. It may also help you to understand the full process by seeing the steps involved.
