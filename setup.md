---
title: Setup
---
If you have not done so already, please follow the instructions [here](https://eic.github.io/tutorial-setting-up-environment/setup.html) well before the start of the tutorial to ensure your system is ready.

This tutorial will go over how to analyze the reconstructed simulation, so you will need to download a few files to work with locally. The files are on the order of 35MB each. For consistency, we will use neutral current DIS events from the December campaign (23.12.0) with minimum Q2 = 10 GeV2 and top beam energy (if you wish to make an energy comparison, you can download additional files). To browse the available files, you can run the following commands from within the eic-shell environment:

```console
xrdfs root://dtn-eic.jlab.org
ls /work/eic2/EPIC/RECO/23.12.0/epic_craterlake/DIS/NC/18x275/minQ2=10
exit
```

Once you have identified the file(s) you wish to copy, (still within eic-shell environment) navigate to the directory you will store your file(s) and run the command:

```console
xrdcp root://dtn-eic.jlab.org//work/eic2/EPIC/RECO/23.12.0/epic_craterlake/DIS/NC/18x275/minQ2=10/filename.eicrecon.tree.edm4eic.root ./
```

Do not forget the trailing ./ (or just . works too) as this tells the progam to put the file in your current dir. 

{% include links.md %}
