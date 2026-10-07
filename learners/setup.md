---
title: Setup
---

If you have not done so already, please follow the
[Setting Up Your Environment tutorial setup instructions](https://eic.github.io/tutorial-setting-up-environment/#setup)
well before the start of the tutorial to ensure your system is ready.

This tutorial will go over how to analyze the reconstructed simulation, so you will need to download a file to work with locally. The files are on the order of 200-350MB each. For consistency, we will use neutral current DIS events from the April 2026 campaign (26.04.1) with minimum Q2 = 10 GeV2 and at the highest electron-proton beam energy combination (if you wish to make an energy comparison, you can download additional files). To browse the available files, you can run the following [Rucio](https://eic.github.io/tutorial-file-access/) command from within the eic-shell environment:

```bash
rucio did content list --short epic:/RECO/26.04.1/epic_craterlake/DIS/NC/18x275/minQ2=10
```

You can download any of the files you want in here. You can do this by (still within eic-shell environment) navigating to the directory you will store your file(s) and run the command:

```bash
xrdcp $(rucio replica list file --protocols root --pfns --rses isopenaccess epic:/RECO/26.04.1/epic_craterlake/DIS/NC/18x275/minQ2=10/pythia8NCDIS_18x275_minQ2=10_beamEffects_xAngle=-0.025_hiDiv_1.0001.eicrecon.edm4eic.root | head -1) ./
```

Do not forget the trailing ./ (or just . works too) as this tells the progam to put the file in your current dir.

::::::::::::::::::::::::::::::::::::::::::::: callout

Note that we can also specify a different filename to copy to as we could with a normal cp command. You might want to do this as the filename is a little cumbersome.
I called mine NC_DIS_18x275_Apr26Campaign.root, just replace ./ with your file name of choice.

:::::::::::::::::::::::::::::::::::::::::::::

This command will download the file (0001) specified. You can of course, download a different file in the same directory if you want.

## Additional Comments for this Tutorial

Note that this tutorial is a little odd in that, for the most part, we don't rely on eic-shell for the majority of the lesson. We will need a working ROOT install though. This ROOT install must also be one that we can work with interactively, with minimal lag. There are two straightforward options for this -

1. eic-shell running locally on your own local machine. You should be able to run ROOT interactively from within eic-shell.
2. A working version (and relatively recent, 6.30 or above, ideally 6.34.02) of ROOT on your local machine.

If you use option 2, note that you will not be able to "stream" files to your ROOT script unless you have xrootd installed too, you will need them available locally. I will be using option 2 for this tutorial.
