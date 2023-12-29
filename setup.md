---
title: Setup
---
In advance of the training session, please ensure that
- You have a GitHub account ([sign up here](https://github.com/signup))
- Your GitHub account is a member of the [EIC organization](https://github.com/eic) on GitHub
  - Email [the EICUG SWG conveners](mailto:eicug-software-conveners@eicug.org) with your GitHub account to be added
- You have some basic familiarity with Unix shell and Git
  - Software Carpentry [Unix Shell](https://swcarpentry.github.io/shell-novice/) training (required)
  - Software Carpentry [Basic Git](https://swcarpentry.github.io/git-novice/) training (recommended)
- You have singularity, apptainer, or docker (on Mac) installed and working
  - `module load singularity` (on most JLab or BNL systems)
  - you should at a minimum be able to run either of the following commands and open a shell:
    - `singularity run docker://alpine`
    - `docker run --rm -it alpine`
  - Software Carpentry [Incubator Singularity](https://carpentries-incubator.github.io/singularity-introduction/) training (optional)
  - HEP Software Foundation [Docker](https://hsf-training.github.io/hsf-training-docker/index.html) training (optional)

{% include links.md %}
