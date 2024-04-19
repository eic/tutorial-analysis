---
title: Setup
---
If you have not done so already, please follow the instructions [here](https://eic.github.io/tutorial-setting-up-environment/setup.html) well before the start of the tutorial to ensure your system is ready.

This tutorial will go over how to analyze the reconstructed simulation, so you will need to download a file to work with locally. The files are on the order of 35MB each. For consistency, we will use neutral current DIS events from the December campaign (23.12.0) with minimum Q2 = 10 GeV2 and top beam energy (if you wish to make an energy comparison, you can download additional files). To browse the available files, you can run the following commands from within the eic-shell environment:

```console
xrdfs root://dtn-eic.jlab.org
ls /work/eic2/EPIC/RECO/23.12.0/epic_craterlake/DIS/NC/18x275/minQ2=10
exit
```

You can download any of the files you want in here. You can do this by (still within eic-shell environment) navigating to the directory you will store your file(s) and run the command:

```console
xrdcp root://dtn-eic.jlab.org//work/eic2/EPIC/RECO/23.12.0/epic_craterlake/DIS/NC/18x275/minQ2=10/pythia8NCDIS_18x275_minQ2=10_beamEffects_xAngle=-0.025_hiDiv_5.0001.eicrecon.tree.edm4eic.root ./
```

Do not forget the trailing ./ (or just . works too) as this tells the progam to put the file in your current dir.

This command will download the file (0001) specified. You can of course, download a different file in the same directory if you want.

## Additional Comments for this Tutorial - April 2024 Version

Note that this tutorial is a little odd in that, for the most part, we don't rely on eic-shell for the majority of the lesson. We will need a working ROOT install though. This ROOT install must also be one that we can work with interactively, with minimal lag. There are two straightforward options for this -

1. eic-shell running locally on your own local machine. You should be able to run ROOT interactively from within eic-shell.
2. A working version (and relatively recent, 6.22 or above) of ROOT on your local machine.

If you use option 2, note that you will not be able to "stream" files to your ROOT script, you will need them available locally. I will be using option 2 for this tutorial.

{% include links.md %}
