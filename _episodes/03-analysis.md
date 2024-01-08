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

> Note: Section to be filled.
{: .callout}

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

Next, we need to access the appropriate branches, we saw how to do this in the "Reading the Output Trees" section. We will need momentum, generator status, and particle species information for the truth particles and momentum information for the reconstructed tracks. The reconstructed track information can be accessed from two different branches: CentralCKFTrackParameters and ReconstructedChargedParticles. We will proceed using the ReconstructedChargedParticles branch as this will give us a chance to practice using associations, copy the following lines into your analysis macro.

```c++
// Set up input file chain
TChain *mychain = new TChain("events");
mychain->Add(infile);

// Initialize reader
TTreeReader tree_reader(mychain);

// Get Particle Information
TTreeReaderArray<int> partGenStat(tree_reader, "MCParticles.generatorStatus");
TTreeReaderArray<float> partMomX(tree_reader, "MCParticles.momentum.x");
TTreeReaderArray<float> partMomY(tree_reader, "MCParticles.momentum.y");
TTreeReaderArray<float> partMomZ(tree_reader, "MCParticles.momentum.z");
TTreeReaderArray<int> partPdg(tree_reader, "MCParticles.PDG");

// Get Reconstructed Track Information
TTreeReaderArray<float> trackMomX(tree_reader, "ReconstructedChargedParticles.momentum.x");
TTreeReaderArray<float> trackMomY(tree_reader, "ReconstructedChargedParticles.momentum.y");
TTreeReaderArray<float> trackMomZ(tree_reader, "ReconstructedChargedParticles.momentum.z");

// Get Associations Between MCParticles and ReconstructedChargedParticles
TTreeReaderArray<unsigned int> recoAssoc(tree_reader, "ReconstructedChargedParticleAssociations.recID");
TTreeReaderArray<unsigned int> simuAssoc(tree_reader, "ReconstructedChargedParticleAssociations.simID");
```

The last two lines encode the association between a ReconstructedChargedParticle and a MCParticle where the matching is determined in the [ParticlesWithPID](https://github.com/eic/EICrecon/blob/main/src/algorithms/pid/ParticlesWithPID.cc) algorithm which generates the ReconstructedChargedParticle objects. 

### Efficiency Analysis

Now that we have access to the data we need we will begin constructing our efficiency plots, starting with efficiency as a function of the true particle pseudorapidity. The basic strategy is outlined below:

1. Loop over all events in the file
2. Within each event, loop over all stable charged particles
3. Identify the ReconstructedChargedParticle (if any) associated with the truth particle we are looking at
4. Create and fill the necessary histograms

Here is the code to implement these steps:

```c++
// Define Histograms
TH1D *partEta = new TH1D("partEta","Eta of Thrown Charged Particles;Eta",100,-5.,5.);
TH1D *matchedPartEta = new TH1D("matchedPartEta","Eta of Thrown Charged Particles That Have Matching Track",100,-5.,5.);

TH1D *matchedPartTrackDeltaR = new TH1D("matchedPartTrackDeltaR","Delta R Between Matching Thrown and Reconstructed Charged Particle",5000,0.,5.);

while(tree_reader.Next()) { // Loop over events

  for(unsigned int i=0; i<partGenStat.GetSize(); i++) // Loop over thrown particles
    {
	    if(partGenStat[i] == 1) // Select stable thrown particles
	      {
	        int pdg = TMath::Abs(partPdg[i]);

	        if(pdg == 11 || pdg == 13 || pdg == 211 || pdg == 321 || pdg == 2212) // Look at charged particles (electrons, muons, pions, kaons, protons)
	          {
		          TVector3 trueMom(partMomX[i],partMomY[i],partMomZ[i]);

		          float trueEta = trueMom.PseudoRapidity();
		          float truePhi = trueMom.Phi();
	    
		          partEta->Fill(trueEta);

              // Loop over associations to find matching ReconstructedChargedParticle
		          for(unsigned int j=0; j<simuAssoc.GetSize(); j++)
		            {
		              if(simuAssoc[j] == i) // Find association index matching the index of the thrown particle we are looking at
		                {
			                TVector3 recMom(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]); // recoAssoc[j] is the index of the matched ReconstructedChargedParticle

                      // Check the distance between the thrown and reconstructed particle
			                float deltaEta = trueEta - recMom.PseudoRapidity();
			                float deltaPhi = TVector2::Phi_mpi_pi(truePhi - recMom.Phi());
			                float deltaR = TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);

			                matchedPartTrackDeltaR->Fill(deltaR);

			                matchedPartEta->Fill(trueEta); // Plot the thrown eta if a matched ReconstructedChargedParticle was found
                    }
                }
            }            
        }
    }
}
```

We should now have everything we need to find the track efficiency as a function of pseudorapidity. To run the macro and produce an output file containing the histograms we defined, simply type `root -l -q trackAnalysis.C`. After the macro runs, you can open the root file to inspect the histograms. The efficiency can be found by taking the ratio of matchedPartEta over partEta.

> Question:
> - We plot the distance between thrown and reconstructed charged partices, does this distribution look reasonable?
> - When filling the matchedPartEta histogram (the numerator in our efficiency), why do we use again the true thrown eta instead of the associated reconstructed eta?
{: .callout}

> Exercise:
> - Find the efficiency as a function of particle momentum. Are there cuts on any other quantities you should place to get a sensible result?
> - Find the efficiency for some 2-D correlations: momentum vs eta; phi vs eta
> - Plot some kinematic distributions (momentum, eta, etc) for all ReconstructedChargedParticles, not just those that are associated with a thrown particle
{: .challenge}

### Resolution Analysis

Next, we will look at track momentum resolution, that is, how well the momentum of the reconstructed track matches that of the thrown particle. We should have all of the "infrastructure" we need in place to do the analysis, we just need to define the appropriate quantities and make the histograms. It only makes sense to define the resolution for tracks and particles which are associated with one another, so we will work within the loop over associations. Define the resolution expression and fill a simple histogram:

```c++
TH1D *trackMomentumRes = new TH1D("trackMomentumRes","Track Momentum Resolution",2000,=10.,10.);
...
// Loop over associations to find matching ReconstructedChargedParticle
for(unsigned int j=0; j<simuAssoc.GetSize(); j++)
  {
    if(simuAssoc[j] == i) // Find association index matching the index of the thrown particle we are looking at
      {
        ...
        double momRes = (recMom.Mag() - trueMom.Mag())/trueMom.Mag();

        trackMomentumRes->Fill(momRes);
      }
  }  
```

While this plot will give us a sense of what the tracking resolution is, we don't expect the resolution to be constant for all momenta or eta.

> Exercise:
> - Make 2-D plots of resolution vs true momentum and vs true pseudorapidity
> - Don't forget to place appropriate cuts on kinematic quantities you are not explicitly plotting
{: .challenge}



