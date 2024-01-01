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
More detailed information on the simulation productions can be found on the [Simulation Production Campaign Website](https://eic.github.io/epic-prod/). 

# Simulation Files Organization

There are three broad classes of files stored on xrootd/S3, each in their own directory:
- EVGEN: The input hepmc3 datasets
- FULL: The full GEANT4 output root files (usually only saved for a fration of runs)
- RECO: The output root files from the reconstruction

Most users will interact with the files in the RECO directory and that is what we will focus on in this tutorial. Within the RECO directory, files are organized by campaign (23.12.0 for the December 2023 campaign, for example), detector configuration, physics process, energy, and Q2. 

{% include links.md %}

