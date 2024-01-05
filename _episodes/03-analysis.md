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

The simulation output trees are "flat" in the sense that there is no event class structure embedded within the tree and no additional libraries are needed to handle the output. Therefore, the end user can simply read the values stored in each branch using whatever method/workflow they are most comfortable with. Examples of several common methods for reading the trees are provided below.

### ROOT TTreeReaderArray

TTreeReader and the associated TTreeReaderArray is a simple interface for reading data from a TTree. The class description and examples can be seen [here](https://root.cern/doc/v630/classTTreeReader.html). To instantiate the reader and access values from a given branch (e.g. the MCParticles branch), one would use the following calls:

```c++
// Set up input file chain
TChain *mychain = new TChain("events");
mychain->Add(input_file_name);

// Initialize reader
TTreeReader tree_reader(mychain);

// Access whatever data-members you need
TTreeReaderArray<int> partGenStat(tree_reader, "MCParticles.generatorStatus");
TTreeReaderArray<float> partMomX(tree_reader, "MCParticles.momentum.x");
...
```

Once you have defined the `TTreeReaderArray` objects for the data-members you want to look at, you can loop over the events and the members within that event:

```c++
while(tree_reader.Next()) { // Loop over events
  for(unsigned int i=0; i<partGenStat.GetSize(); i++) // Loop through particles in the event
    {
      int particleStatus = partGenStat[i]; // Access data-members as you would an array
      float particleXMomentum = partMomX[i]; // partMomX should have same number of entries as partGenStat because they are in the same branch
      ...
    }
}
```
All members of the same branch should have the same number of entries, so it is sufficient to use any member of the branch to set the limit of your loop.


### ROOT RDataFrames

### PYTHON

> Note: Section to be filled.
{: .callout}

## The MCParticles Record

Nearly every analysis will include some comparison to the truth level, so it is important to understand how to access generator level information

## Sample Analysis: Track Efficiency and Resolution

