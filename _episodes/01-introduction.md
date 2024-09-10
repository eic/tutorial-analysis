---
title: "Introduction"
teaching: 15
exercises: 5
questions:
- "How do I locate and access the simulation output?"
objectives:
- "Understand how the simulation output is organized"
- "Know how to access the simulation output using Jefferson Lab xrootd"
- "Know how to access the simulation output using BNL S3"
keypoints:
- "Use `xrdfs` from within the eic-shell to browse available files from simulations campaigns."
- "Use `xrdcp` from within eic-shell to copy files to your local environment."
- "Alternatively, use the minio client to access simulation files."
- "Within eic-shell, you can also stream files directly in your root macros."
---

More detailed information on the simulation productions, including the information presented below, can be found on the [Simulation Production Campaign Website](https://eic.github.io/epic-prod/). 

## Simulation Files Organization

There are three broad classes of files stored on xrootd/S3, each in their own directory:
- EVGEN: The input hepmc3 datasets
    - E.g. some files that have been supplied by a physics event generator
- FULL: The full GEANT4 output root files (usually only saved for a fraction of runs)
    - If running a simulation yourself, this would be your output from processing npsim
- RECO: The output root files from the reconstruction
    - And again, if running yourself, this would be your output from EICrecon (after you've used your awesome new reconstruction algorithm from the later tutorial of course)

Most users will interact with the files in the RECO directory and that is what we will focus on in this tutorial. Within the RECO directory, files are organized by campaign (24.04.0 for the April 2024 campaign, for example), detector configuration and then physics process. Each physics process will have different sub directories, for example generator version, energy, or Q2. The directory structure and number of reconstructed files for each campaign can be found on the Simulation Website [here](https://eic.github.io/epic-prod/campaigns/campaigns_reco.html).

## Access Simulation from Jefferson Lab xrootd

The prefered method for browsing the simulation output is to use xrootd from within the eic-shell. To browse the directory structure and exit, one can run the commands:
```console
./eic-shell
xrdfs root://dtn-eic.jlab.org
ls /work/eic2/EPIC/RECO/24.04.0
exit
```
It is also possible to copy a file and open it locally using the `xrdcp` command:
```console
./eic-shell
xrdcp root://dtn-eic.jlab.org//work/eic2/EPIC/RECO/24.04.0/path-to-file .
exit
```
## Access Simulation from BNL S3

The simulation files can also be accessed from S3 storage at BNL using the MinIO client for S3 storage. It is included in eic-shell. To install it natively, you can issue the following commands to install minio:
```console
mkdir --parent ~/bin
curl https://dl.min.io/client/mc/release/linux-amd64/mc --create-dirs -o ~/bin/mc
chmod +x ~/bin/mc
```
From here on out, we assume `mc` is in your PATH variable, otherwise you can use the full path, in the above example `~/bin/mc`.
After the client is installed, it needs to be configured for read access:
```console
export S3_ACCESS_KEY=<credential>; export S3_SECRET_KEY=<credential>
mc config host add S3 https://eics3.sdcc.bnl.gov:9000 $S3_ACCESS_KEY $S3_SECRET_KEY
```
The <credential> for read access values can be obtained by asking on Mattermost. Assuming the minio client is installed and configured as above, one can browse the file structure using the minio `ls` command:
```console
mc ls S3/eictest/EPIC/RECO
```

Files can also be coppied locally by replacing `ls` with `cp`.

## Streaming Files

It is also possible to open a file directly in ROOT. Note that the following command should be executed after opening root and `TFile::Open()` should be used:
```console
auto f = TFile::Open("root://dtn-eic.jlab.org//work/eic2/EPIC/RECO/path-to-file")
```
or alternatively
```console
auto f = TFile::Open("s3https://eics3.sdcc.bnl.gov:9000/eictest/EPIC/RECO/path-to-file");
```

## Reminder - Download a file for the next step!

We will need a file to analyse going forward, if you have not done so, download a file now!

Grab a file from -

```console
/work/eic2/EPIC/RECO/24.07.0/epic_craterlake/DIS/NC/18x275/minQ2=10/
```
For example -

```console
xrdcp root://dtn-eic.jlab.org//work/eic2/EPIC/RECO/24.04.0/epic_craterlake/DIS/NC/18x275/minQ2=10/pythia8NCDIS_18x275_minQ2=10_beamEffects_xAngle=-0.025_hiDiv_5.0001.eicrecon.tree.edm4eic.root ./
```
Note that the ./ at the end is the target location to copy to. Change this as desired.

You can also stream the file if you prefer, just copy the path of the file above. You will need to modify the scripts later in the tutorial accordingly to account for this.

## Advanced Use Case - Grabbing a whole bunch of files

I won't go through this in the tutorial, but this may be something you want to come back to as you get deeper into writing and using your own analysis code. This advanced use case involves copying/using a large number of processed files. Something you might want to do once your analysis is out of the testing phase and onto the "Let's process ALL of the data!" stage.

If you're moving a lot of files around, you might normally resort to using a wildcard -

```console
cp File* My_Folder/
```

or similar. However, with the mc or xrdcp, this isn't so trivial. Some methods to test and try are include below. 

```console
mc find S3/eictest/EPIC/RECO/main/CI/ --name "*0001*.root" --exec "mc cp {} ."
```
where here we're finding things in the given path that match the name pattern provided, and copying them to our current directory.

Alternatively, you could grab a list of the files you want and pipe them to a file -

```console
xrdfs root://dtn-eic.jlab.org ls /work/eic2/EPIC/RECO/24.04.0/epic_craterlake/DIS/NC/18x275/minQ2=10 | sed 's|^|root://dtn-eic.jlab.org/|g' > list.txt
```

In this case, we're listing all files on the server in that path, piping them to sed and inserting "root://dtn-eic.jlab.org/" at the front and then feeding the output to the file "list.txt".

```console
more list.txt
root://dtn-eic.jlab.org//work/eic2/EPIC/RECO/24.04.0/epic_craterlake/DIS/NC/18x275/minQ2=10/pythia8NCDIS_18x275_minQ2=10_beamEffects_xAngle=-0.025_hiDiv_1.0000.eicrecon.tree.edm4eic.root
root://dtn-eic.jlab.org//work/eic2/EPIC/RECO/24.04.0/epic_craterlake/DIS/NC/18x275/minQ2=10/pythia8NCDIS_18x275_minQ2=10_beamEffects_xAngle=-0.025_hiDiv_1.0001.eicrecon.tree.edm4eic.root
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
