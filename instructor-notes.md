---
title: "Instructor Notes"
---

This tutorial teaches learners how to analyze the reconstructed simulation output produced by
EICrecon. It follows the [Setting Up Your Environment](https://eic.github.io/tutorial-setting-up-environment/)
tutorial, which establishes the `eic-shell` environment learners rely on here.

## Before the session

- Ask learners to complete the [Setup](../learners/setup.md) page in advance, in particular
  downloading a reconstruction output file (50-80 MB) from the current simulation campaign, since
  this depends on network access to the JLab xrootd server.
- Confirm learners have a working, interactive ROOT install (either inside `eic-shell` or ROOT
  6.30+ locally) as described on the Setup page.

## Timing

Most of the lesson is hands-on analysis. Learners should work the exercises in whichever language
(ROOT C++, Python/uproot, RDataFrame, or PODIO) they are most comfortable with; you do not need to
cover every approach live.
