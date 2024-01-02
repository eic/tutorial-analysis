---
title: "Introduction"
teaching: 20
exercises: 0
questions:
- "How do I locate and access the simulation output?"
objectives:
- "Understand how the simulation output is organized"
- "Know how to access the simulation output using Jefferson Lab xrootd"
- "Know how to access the simulation output using BNL S3"
keypoints:
- "Use `xrdfs` from within the eic shell or the minio client to access simulation files"
---

More detailed information on the simulation productions, including the information presented below, can be found on the [Simulation Production Campaign Website](https://eic.github.io/epic-prod/). 

## Simulation Files Organization

There are three broad classes of files stored on xrootd/S3, each in their own directory:
- EVGEN: The input hepmc3 datasets
- FULL: The full GEANT4 output root files (usually only saved for a fration of runs)
- RECO: The output root files from the reconstruction

Most users will interact with the files in the RECO directory and that is what we will focus on in this tutorial. Within the RECO directory, files are organized by campaign (23.12.0 for the December 2023 campaign, for example), detector configuration, physics process, energy, and Q2. The directory structure and number of reconstructed files for each campaign can be found on the Simulation Website [here](https://eic.github.io/epic-prod/campaigns/campaigns_reco.html).

## Access Simulation from Jefferson Lab xrootd

The prefered method for browsing the simulation output is to use xrootd from within the eic-shell. To browse the directory structure and exit, one can run the commands:
```console
./eic-shell
xrdfs root://dtn-eic.jlab.org
ls /work/eic2/EPIC/RECO/23.12.0
exit
```

It is alos possible to copy a file and open it locally using the `xrdcp` command:
```console
./eic-shell
xrdcp root://dtn-eic.jlab.org//work/eic2/EPIC/RECO/23.12.0/path-to-file .
exit
```

## Access Simulation from BNL S3

The simulation files can also be accessed from S3 storage at BNL using the client. Issue the following commands to install minio:
```console
mkdir --parent ~/bin
curl https://dl.min.io/client/mc/release/linux-amd64/mc --create-dirs -o ~/bin/mc
chmod +x ~/bin/mc
```

After the client is installed, it needs to be configured for read access:
```console
~/bin/mc config host add S3 https://eics3.sdcc.bnl.gov:9000 <credential> <credential>
```

The <credential> values can be obtained by asking on Mattermost. Assuming the minio client is installed and configured as above, one can browse the file structure using the minio `ls` command:
```console
~/bin/mc ls S3/eictest/EPIC/RECO
```

Files can also be coppied locally by replacing `ls` with `cp`.

## Streaming Files

It is also possible to open a file directly in ROOT. Note that the following command should be executed after opening root and `TFile::Open()` should be used:
```console
auto f = TFile::Open("root://dtn-eic.jlab.org//work/eic2/EPIC/RECO/path-to-file")
```

{% include links.md %}

