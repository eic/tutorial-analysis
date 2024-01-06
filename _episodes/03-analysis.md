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

The branches and their members can be viewed by opening a file with TBrowser (`new TBrowser()`) from within ROOT. Once you have defined the `TTreeReaderArray` objects for the data-members you want to look at, you can loop over the events and the members within that event:

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

Nearly every analysis will include some comparison to the truth level, so it is important to understand how to access generator level information. Particles produced by the Monte Carlo Event Generator and during the interaction of these primary particles with the detector material as modeled by GEANT are stored in the MCParticles branch, whoes structure is defined by the datatype [edm4hep::MCParticle](https://github.com/key4hep/EDM4hep/blob/main/edm4hep.yaml#L140). The particle's [PDG](https://pdg.lbl.gov/2020/reviews/rpp2020-rev-monte-carlo-numbering.pdf) code, charge, production time, mass, production vertex, endpoint, momentum at the production vertex, and momentum at the endpoint are all available. In addition, the status of the particle as defined by the event generator and the detector simulation are stored. For example, if one wanted to look at stable particles from the event generator, they would require `MCParticles.generatorStatus == 1`. The field `MCParticles.simulatorStatus` is a bit-field which encodes some information on how the particle propagated through the detector. The detailed definition of the bit assignments can be found in the [edm4hep yaml file](https://github.com/key4hep/EDM4hep/blob/main/edm4hep.yaml#L140).

## Sample Analysis: Track Efficiency and Resolution

As a sample exercise to become familiar with the simulation output and how to use it in a realistic setting, we will find the tracking eficiency and resolution. We will need to access the reconstructed track information and the truth particle information and we will have to associate the individual tracks and particles to one another. 

Before we begin, we should create a skeleton macro to handle file I/O. For the `TTreeReader` and `RDataFrames` examples, we will use a simple ROOT macro. Using your favorite editor, create a file with a name like `trackAnalysis.C` or something similar and copy in the following code: 

```c++
void trackAnalysis(TString infile="path_to_your_simu_file")
  {
    // Set output file for the histograms
    TFile *ofile = TFile::Open("out.hist.root","RECREATE");

    // Analysis code will go here

    ofile->Write(); // Write histograms to file
    ofile->Close(); // Close output file
  }
```

Next, we need to access the appropriate branches, we saw how to do this in the "Reading the Output Trees" section. We will need momentum, generator status, and possibly particle species information for the truth particles and momentum information for the reconstructed tracks. The reconstructed track information can be accessed from two different branches: CentralCKFTrackParameters and ReconstructedChargedParticles. 

