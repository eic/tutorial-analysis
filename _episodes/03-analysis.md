---
title: "Analyzing the Reconstruction Output"
teaching: 30
exercises: 0
questions:
- "How does one utilize the reconstruction output trees to do an analysis?"
objectives:
- "Become familiar with methods for reading the trees"
- "Understand how to access truth/particle information"
- "Perform some basic analyses"
keypoints:
- "FIXME"
---

Discussion of steps needed to perform and analysis ...

## Reading the Output Trees

The simulation output trees are "flat" in the sense there is no event class structure embedded within the tree and no additional libraries are needed to handle the output. Therefore, the end user can simply read the values stored in each branch using whatever method/workflow they are most comfortable with. Examples of several common methods for reading the trees are provided below.

### ROOT TTreeArrayReader

### ROOT RDataFrames

### PYTHON

> Note: Section to be filled.
{: .callout}

## The MCParticles Record

Nearly every analysis will include some comparison to the truth level, so it is important to understand how to access generator level information

## Sample Analysis: Track Efficiency and Resolution

