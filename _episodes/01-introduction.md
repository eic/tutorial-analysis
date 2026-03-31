---
title: "Introduction"
teaching: 5
exercises: 5
questions:
- "How do I locate and access the simulation output?"
objectives:
- "Understand how the simulation output is organised."
- "Download a file for the next step of the tutorial."
keypoints:
- "Use `xrdcp` from within eic-shell to copy files to your local environment."
---

More detailed information on the simulation productions, including the information presented below, can be found on the [Simulation Production Campaign Website](https://eic.github.io/epic-prod/). 

**Note that as of March 2026, Rucio will soon become the default and preferred method to browse and find files and datasets. A tutorial on using Rucio for this purpose will be presented soon, please see [here](https://eic.github.io/tutorial-file-access/) for the latest version of this tutorial.**

## Simulation Files Organization

There are three broad classes of files stored on xrootd, each in their own directory:
- EVGEN: The input hepmc3 datasets
    - E.g. some files that have been supplied by a physics event generator
- FULL: The full GEANT4 output root files (usually only saved for a fraction of runs)
    - If running a simulation yourself, this would be your output from processing npsim
- RECO: The output root files from the reconstruction
    - And again, if running yourself, this would be your output from EICrecon (after you've used your awesome new reconstruction algorithm from the later tutorial of course)

Most users will interact with the files in the RECO directory and that is what we will focus on in this tutorial. Within the RECO directory, files are organized by campaign (26.02.0 for the February 2026 campaign, for example), detector configuration and then physics process. Each physics process will have different sub directories, for example generator version, energy, or Q2. The directory structure and number of reconstructed files for each campaign can be found on the Simulation Website [here](https://eic.github.io/epic-prod/campaigns/campaigns_reco.html).

> **Note that campaigns more than ~6 months old will not directly be accessible.**
> **If you are running this tutorial and encounter a file access error, check the campaign you are trying to access.**
> **Where possible, use the latest campaign available.**
{: .callout}

## Download a file for the next step!

We will need a file to analyse going forward, if you have not done so, download a file now!

Grab a file from -

```console
epic:/RECO/26.02.0/epic_craterlake/DIS/BeAGLE1.03.02-1.2/eHe3/10x110/q2_2to10/
```

> Reminder, you can check the *content* of files within this dataset via:
> ```bash
> rucio did content list --short epic:/RECO/26.02.0/epic_craterlake/DIS/BeAGLE1.03.02-1.2/eHe3/10x110/q2_2to10
> ```
> and check the location of files in the dataset via:
> ```bash
> rucio replica list file --protocols root --pfns --rses isopenaccess epic:/RECO/26.02.0/epic_craterlake/DIS/BeAGLE1.03.02-1.2/eHe3/10x110/q2_2to10
> ```
{: .callout}

For example -

```console
xrdcp root://dtn-eic.jlab.org:1094//volatile/eic/EPIC//RECO/26.02.0/epic_craterlake/DIS/BeAGLE1.03.02-1.2/eHe3/10x110/q2_2to10/BeAGLE1.03.02-1.2_DIS_eHe3_10x110_q2_2to10_ab.0001.eicrecon.edm4eic.root
```
Note that the ./ at the end is the target location to copy to. Change this as desired.

> Note that we can also specify a different filename to copy to as we could with a normal cp command. You might want to do this as the filename is a little cumbersome.
> I called mine `3He_10x110_Feb26Campaign.root`, just replace ./ with your file name of choice.
{: .callout}

You can also stream the file if you prefer, just copy the path of the file above. You will need to modify the scripts later in the tutorial accordingly to account for this. Check the [File Access Tutorial](https://eic.github.io/epic-prod/) for information and examples on how to do this.

> Typically, if you are processing more than a handful of files, it is probably best to stream files from the server rather than downloading a local copy of all files.
{: .callout}
<!---
## Advanced Use Case - Grabbing a whole bunch of files

I won't go through this in the tutorial, but this may be something you want to come back to as you get deeper into writing and using your own analysis code. This advanced use case involves copying/using a large number of processed files. Something you might want to do once your analysis is out of the testing phase and onto the "Let's process ALL of the data!" stage.

If you're moving a lot of files around, you might normally resort to using a wildcard -

```console
cp File* My_Folder/
```

or similar. However, with xrdcp, this isn't so trivial. Some methods to test and try are include below. 

where here we're finding things in the given path that match the name pattern provided, and copying them to our current directory.

Alternatively, you could grab a list of the files you want and pipe them to a file -

```console
xrdfs root://dtn-eic.jlab.org ls /volatile/eic/EPIC/RECO/26.02.0/epic_craterlake/DIS/NC/18x275/minQ2=10 | sed 's|^|root://dtn-eic.jlab.org/|g' > list.txt
```

In this case, we're listing all files on the server in that path, piping them to sed and inserting "root://dtn-eic.jlab.org/" at the front and then feeding the output to the file "list.txt".

```console
more list.txt
root://dtn-eic.jlab.org//volatile/eic/EPIC/RECO/26.02.0/epic_craterlake/DIS/NC/18x275/minQ2=10/pythia8NCDIS_18x275_minQ2=10_beamEffects_xAngle=-0.025_hiDiv_1.0000.eicrecon.tree.edm4eic.root
root://dtn-eic.jlab.org//volatile/eic/EPIC/RECO/26.02.0/epic_craterlake/DIS/NC/18x275/minQ2=10/pythia8NCDIS_18x275_minQ2=10_beamEffects_xAngle=-0.025_hiDiv_1.0001.eicrecon.tree.edm4eic.root
...
```
We could then, for example, feed this list to a TChain -

```console
TChain events("events")
std::ifstream in("list.txt")
std::string file("")
while (in >> file) events.Add(file.data())
events.Scan("@MCParticles.size()","","",10)
```
Where in the final line we're only going to skim over the first 10 events.

It should be noted that the best solution may just be to run the files from the server, rather than copying them to somewhere else and running them there.
-->
