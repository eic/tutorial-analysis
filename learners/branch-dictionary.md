---
title: "Branch Dictionary"
---

Data files from EICrecon can sometimes look a little overwhelming with the large number of branches included. However, for the exercises in this tutorial, we end up needing only a small subset of the branches and the information they contain. A brief overview of some commonly used branches is included below. You might find it useful to refer to these as you're constructing your analysis.

## Detector Branches

The naming convention of many of the branches is based around the detector they correspond to. For example, take the branch -

```bash
EcalEndcapNClusters
```

Breaking this down piece by piece, the branch tells us at the start which detector it corresponds to -

```bash
EcalEndcap
```

ECal -> Electromagnetic calorimeter, Endcap -> The endcap (as opposed to the barrel). The second part tells us a little bit more about which information pertaining to that detector is in this branch -

```bash
...NClusters
```

Clusters tells us this is probably something to do with the clusters of hits in this calorimeter. The N is telling is that this is for negatively charged clusters in this detector. If we open this branch, we see the sort of information we might expect to be associated with calorimeter clusters, for example -

```bash
EcalEndcapNClusters.energy
EcalEndcapNClusters.nhits
EcalEndcapNClusters.position.x
```

So stored in these branches is event by event information on the energy of clusters in this detector, the number of hits per cluster and the x position of those clusters. Remember, each event may have more than one cluster so this branch is actually an array of values for each event.

## MCParticles

The MCParticles branch contains truth level information on the input events into the simulation and reconstruction. Many simulation studies will involve a comparison to this "true" information. This branch contains event level information for many quantities  -

```bash
MCParticles.PDG
MCParticles.charge
MCParticles.vertex.x
MCParticles.momentum.x
...
```

The PDG value tells us what the input particle actually was, we can cross-reference the output of this with the [PDG code](https://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf) to see what our inputs actually were. We can also easily extract the charge, vertex (x component) and momentum (x component) of the particles in our event using the leaves above.

## ReconstructedParticles

This branch contains information on all of the particles that EICrecon reconstructed in our simulation. This branch includes all reconstructed particles, charged or not. It contains very similar information to the MCParticles branch, but also includes some other key information. For example, since these particles are reconstructed, their PDG is assigned based upon what EICrecon thinks they are. As such, we also get information on quantities such as GoodnessOfPID in this branch.

## ReconstructedChargedParticles

Similar to reconstructed particles. In this case though, we only get information on reconstructed particles that are charged.

## ReconstructedParticlesLinks

This branch contains information on links between MC truth information and reconstructed particles. In this case, this is for all reconstructed particles. We can find the link index matching the index of the reconstructed particle and compare what we have reconstructed to its actual "true" information. See the examples in [episode 3](../episodes/03-analysis.md) for some example usage of this.

## ReconstructedChargedParticlesLinks

Similar to ReconstructedParticlesLinks, but for reconstructed charged particles.
