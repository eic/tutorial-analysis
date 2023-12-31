---
title: "The Reconstruction Output Tree"
teaching: 30
exercises: 0
questions:
- "How is the reconstruction output tree populated?"
objectives:
- "Become familiar with the tree branches"
- "Locate the EICrecon factory/algorithm used to fill a specific branch"
- "Become familiar with the edm4hep and edm4eic data models"
- "Understand associations and relations"
keypoints:
- "FIXME"
---

What we generally call the simulation output are root trees generated by the reconstruction software, [EICrecon](https://github.com/eic/EICrecon/tree/main). The branches which appear in the output trees and their content are determined by which EICrecon factories and algorithms are run. 

## EICrecon Output Tree Structure

The output tree contains various branches appropriate for the individual detector subsystems. For example, the calorimeter subsystems will have branches for reconstructed hits and clusters. In addition to individual subsystem information, there are also branches for reconstructed quantities (e.g. reconstructed particles, inclusive kinematics, and jets) which may combine information from several subsystems. There are also branches encoding relationships between different reconstructed quantities as well as reconstructed and truth quantities. 

> Exercise
> - Stream a simulation output tree from within root (see previous lesson) and browse the structure by calling `new TBrowser()`
> - Take some time to explore the branches of the tree. What information is included for various subsystems? What are some of the reconstructed quantities?
> - Try plotting some basic quantities. Plot the cluster energy of the negative endcap ECal - do you see the peak from the scattered electron?
{: .challenge}

## EICrecon Factories and Algorithms

## Data Models

## Associations and Relations

