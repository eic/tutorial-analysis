---
title: "Analyzing the Reconstruction Output"
teaching: 20
exercises: 40
questions:
- "How does one utilize the reconstruction output trees to do an analysis?"
objectives:
- "Become familiar with methods for reading the trees"
- "Understand how to access truth/particle information"
- "Find track efficiency and resolution"
keypoints:
- "Flat tree structure provides flexibility in analysis."
- "The ReconstructedChargedParticles branch holds information on reconstructed tracks."
---

So far, we have only looked at (and plotted) some information from our file interactively. This is very useful and can help us identify the variables we want to deal with. However, we can't really use these techniques to conduct a full analysis of the data. To do so, we typically use a script or macro. In this part of the tutorial, we will create a script that we can use to do a relatively straightforward analysis of our file.

## Reading the Output Trees

The simulation output trees are "flat" in the sense that there is no event class structure embedded within the tree and no additional libraries are needed to handle the output. Therefore, the end user can simply read the values stored in each branch using whatever method/workflow they are most comfortable with. Examples of several common methods for reading the trees are provided below. We will see a ROOT TTreeReader based example using a ROOT macro and a python/uproot based version. During the tutorial, you should try the exercise using whichever language you feel most comfortable with.

## Sample Analysis with ROOT TTreeReader: Track Efficiency and Resolution

As a sample exercise to become familiar with the simulation output and how to use it in a realistic setting, we will find the tracking eficiency and resolution. We will need to access the reconstructed track information and the truth particle information and we will have to associate the individual tracks and particles to one another. 

Before we begin, we should create a skeleton macro to handle file I/O. For the `TTreeReader` example, we will use a simple ROOT macro. Using your favorite editor, create a file with a name like `trackAnalysis.C` or something similar and copy in the following code: 

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

We will need momentum, generator status, and particle species information for the truth particles and momentum information for the reconstructed tracks. The reconstructed track information can be accessed from two different branches: CentralCKFTrackParameters and ReconstructedChargedParticles. We can access these branches using a TTreeReaderArray.

> ROOT TTreeReaderArrays:
>
>TTreeReader and the associated TTreeReaderArray is a simple interface for reading data from a TTree. The class description and examples can be seen [here](https://root.cern/doc/v630/classTTreeReader.html). To instantiate the reader and access values from a given branch (e.g. the MCParticles branch), one would use the following calls:
>
> ```c++
> // Set up input file chain
> TChain *mychain = new TChain("events");
> mychain->Add(input_file_name);
>
> // Initialize reader
> TTreeReader tree_reader(mychain);
>
> // Access whatever data-members you need
> TTreeReaderArray<int> partGenStat(tree_reader, "MCParticles.generatorStatus");
> TTreeReaderArray<float> partMomX(tree_reader, "MCParticles.momentum.x");
> ...
> ```
>
> The branches and their members can be viewed by opening a file with TBrowser (`new TBrowser()`) from within ROOT. Once you have defined the `TTreeReaderArray` objects for the data-members you want to look at, you can loop over the events and the members within that event:
>
> ```c++
> while(tree_reader.Next()) { // Loop over events
>  for(unsigned int i=0; i<partGenStat.GetSize(); i++) // Loop through particles in the event
>    {
>      int particleStatus = partGenStat[i]; // Access data-members as you would an array
>      float particleXMomentum = partMomX[i]; // partMomX should have same number of entries as partGenStat because they are in the same branch
>      ...
>    }
> }
> ```
> All members of the same branch should have the same number of entries, so it is sufficient to use any member of the branch to set the limit of your loop.
{: .callout}

We will proceed using the ReconstructedChargedParticles branch as this will give us a chance to practice using associations, copy the following lines into your analysis macro.

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

> Hint:
> Refer to [the script template](https://eic.github.io/tutorial-analysis/exercise_scripts/index.html#efficiencyanalysisc) if you're having trouble putting things in the right place.
{: .callout}

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
> - Do the histogram ranges make sense?
> - We plot the distance between thrown and reconstructed charged partices, does this distribution look reasonable?
> - When filling the matchedPartEta histogram (the numerator in our efficiency), why do we use again the true thrown eta instead of the associated reconstructed eta?
{: .callout}

> Exercise:
> - Find the efficiency as a function of particle momentum. Are there cuts on any other quantities you should place to get a sensible result?
> - Find the efficiency for some 2-D correlations: momentum vs eta; phi vs eta
> - Plot some kinematic distributions (momentum, eta, etc) for all ReconstructedChargedParticles, not just those that are associated with a thrown particle
{: .challenge}

### Resolution Analysis

> Hint:
> Refer to [the script template](https://eic.github.io/tutorial-analysis/exercise_scripts/index.html#resolutionanalysisc) if you're having trouble putting things in the right place.
{: .callout}

Next, we will look at track momentum resolution, that is, how well the momentum of the reconstructed track matches that of the thrown particle. We should have all of the "infrastructure" we need in place to do the analysis, we just need to define the appropriate quantities and make the histograms. It only makes sense to define the resolution for tracks and particles which are associated with one another, so we will work within the loop over associations. Define the resolution expression and fill a simple histogram:

```c++
TH1D *trackMomentumRes = new TH1D("trackMomentumRes","Track Momentum Resolution",2000,-10.,10.);
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

While this plot will give us a sense of what the tracking resolution is, we don't expect the resolution to be constant for all momenta or eta. We can get a more complete picture by plotting the resolution as a function of different kinematic quantities. 

> Exercise:
> - Make 2-D plots of resolution vs true momentum and vs true pseudorapidity.
> - Break resolution plots down by particle species.
{: .challenge}

> Question:
> - Will the histogram ranges for each particle species be the same?
> - Could we present the resolution values in a more understandable way?
{: .callout}

## Sample Analysis with Python/uproot: Track Efficiency and Resolution

If you are more familiar with python than you are with C/C++, you might find that using a python based root macro is easier for you. Outlined below are sample blocks of code for creating and running a python based analysis script.

With python, some tasks become easier, e.g. string manipulation and writing to (non ROOT) files.

Before we begin, we should create a skeleton macro to handle file I/O. For this example, we will make a simple python script. Using your favorite editor, create a file with a name like `trackAnalysis.py` or something similar and copy in the following code: 

```python
#! /usr/bin/python
         
#Import relevant packages
import ROOT, math, array                                
from ROOT import TH1F, TH2F, TMath, TTree, TVector3, TVector2
import uproot as up

# Define and open files
infile="PATH_TO_FILE" 
ofile=ROOT.TFile.Open("TrackAnalysis_OutPy.root", "RECREATE")

# Open input file and define branches we want to look at with uproot
events_tree = up.open(infile)["events"]

# Define histograms below

# Add main analysis loop(s) below

# Write output histograms to file below                

# Close files
ofile.Close()                    
```
Note that we are using the module uproot to access the data here. See [further documentation here](https://masonproffitt.github.io/uproot-tutorial/03-trees/index.html). You may also need some of the other included packages too.

> We will use uproot a little bit like we use the TTreeReader in the other example. We can define the branches we want and assign them to arrays with uproot.
> We can do this via:
>  ```python
> # Open input file and define branches we want to look at with uproot
> events_tree = up.open(infile)["events"]
> # Get particle information# Get particle information
> partGenStat = events_tree["MCParticles.generatorStatus"].array()
> partMomX = events_tree["MCP articles.momentum.x"].array() 
> partMomY = events_tree["MCParticles.momentum.y"].array()
> partMomZ = events_tree["MCParticles.momentum.z"].array()
> partPdg = events_tree["MCParticles.PDG"].array()
>
> # Get reconstructed track information
> trackMomX = events_tree["ReconstructedChargedParticles.momentum.x"].array()
> trackMomY = events_tree["ReconstructedChargedParticles.momentum.y"].array()
> trackMomZ = events_tree["ReconstructedChargedParticles.momentum.z"].array()
>  ...
> ```
>  We can then access them as an array in a loop -
> ```python
> # Add main analysis loop(s) below
> for i in range(0, len(events_tree)): # Loop over all events
>     for j in range(0, len(partGenStat[i])): # Loop over all thrown particles
>         if partGenStat[i][j] == 1: # Select stable particles
>             pdg = abs(partPdg[i][j]) # Get PDG for each stable particle
>             ...
> ```
> Uproot effectively takes the information in the tree, and turns it into an array. We can then acces and manipulate this array in the same way that we can with any array in python.
>
> Note that if you are using an older version of uproot (v2.x.x), you will need to access the branches slightly differently via -
> ```python
> partGenStat = events_tree.array("MCParticles.generatorStatus")
> ```
{: .callout}

You can run this file with ``python3 trackAnalysis.py``. It should open your file and create an empty output root file as specified. We will add histograms to this script and fill them in the next step.

Note that depending upon your setup, ``python trackAnalysis.py`` may work too.

### Efficiency Analysis

> Hint:
> Refer to [the script template](https://eic.github.io/tutorial-analysis/exercise_scripts/index.html#efficiencyanalysispy) if you're having trouble putting things in the right place.
{: .callout}

As with the ROOT TTreeReader example, we will find the tracking eficiency and resolution. We will need to access the reconstructed track information and the truth particle information and we will have to associate the individual tracks and particles to one another.

The basic strategy is the same:

1. Loop over all events in the file
2. Within each event, loop over all stable charged particles
3. Identify the ReconstructedChargedParticle (if any) associated with the truth particle we are looking at
4. Create and fill the necessary histograms

Here is the sample code to implement these steps:

```python
# Get assocations between MCParticles and ReconstructedChargedParticles
recoAssoc = events_tree["ReconstructedChargedParticleAssociations.recID"].array()
simuAssoc = events_tree["ReconstructedChargedParticleAssociations.simID"].array()

# Define histograms below
partEta = ROOT.TH1D("partEta","Eta of Thrown Charged Particles;Eta",100, -5 ,5 )
matchedPartEta = ROOT.TH1D("matchedPartEta","Eta of Thrown Charged Particles That Have Matching Track", 100, -5 ,5);
matchedPartTrackDeltaR = ROOT.TH1D("matchedPartTrackDeltaR","Delta R Between Matching Thrown and Reconstructed Charge Particle", 5000, 0, 5);

# Add main analysis loop(s) below
for i in range(0, len(events_tree)): # Loop over all events
    for j in range(0, len(partGenStat[i])): # Loop over all thrown particles
        if partGenStat[i][j] == 1: # Select stable particles
            pdg = abs(partPdg[i][j]) # Get PDG for each stable particle
            if(pdg == 11 or pdg == 13 or pdg == 211 or pdg == 321 or pdg == 2212):
                trueMom = ROOT.TVector3(partMomX[i][j], partMomY[i][j], partMomZ[i][j])
                trueEta = trueMom.PseudoRapidity()
                truePhi = trueMom.Phi()
                
                partEta.Fill(trueEta)
                for k in range(0,len(simuAssoc[i])): # Loop over associations to find matching ReconstructedChargedParticle
                    if (simuAssoc[i][k] == j):
                        recMom = ROOT.TVector3(trackMomX[i][recoAssoc[i][k]], trackMomY[i][recoAssoc[i][k]], trackMomZ[i][recoAssoc[i][k]])
                        deltaEta = trueEta - recMom.PseudoRapidity()
                        deltaPhi = TVector2. Phi_mpi_pi(truePhi - recMom.Phi())
                        deltaR = math.sqrt((deltaEta*deltaEta) + (deltaPhi*deltaPhi))

                        matchedPartEta.Fill(trueEta)
                        matchedPartTrackDeltaR.Fill(deltaR)
                        
# Write output histograms to file below
partEta.Write()
matchedPartEta.Write()
matchedPartTrackDeltaR.Write()

# Close files
ofile.Close()
```
Insert this block of code appropriately. We should now have everything we need to find the track efficiency as a function of pseudorapidity. Run the script with `python3 trackAnalysis.py``. This should produce a root file with a few histograms in place. The efficiency can be found by taking the ratio of matchedPartEta over partEta.

> Question:
> - Do the hisotgram ranges make sense?
> - We plot the distance between thrown and reconstructed charged partices, does this distribution look reasonable?
> - When filling the matchedPartEta histogram (the numerator in our efficiency), why do we use again the true thrown eta instead of the associated reconstructed eta?
{: .callout}

> Exercise:
> - Find the efficiency as a function of particle momentum. Are there cuts on any other quantities you should place to get a sensible result?
> - Find the efficiency for some 2-D correlations: momentum vs eta; phi vs eta
> - Plot some kinematic distributions (momentum, eta, etc) for all ReconstructedChargedParticles, not just those that are associated with a thrown particle
{: .challenge}

### Resolution Analysis

> Hint:
> Refer to [the script template](https://eic.github.io/tutorial-analysis/exercise_scripts/index.html#resolutionanalysispy) if you're having trouble putting things in the right place.
{: .callout}

Next, we will look at track momentum resolution, that is, how well the momentum of the reconstructed track matches that of the thrown particle. We should have all of the "infrastructure" we need in place to do the analysis, we just need to define the appropriate quantities and make the histograms. It only makes sense to define the resolution for tracks and particles which are associated with one another, so we will work within the loop over associations. Define the resolution expression and fill a simple histogram by inserting this block of code appropriately:

```python
trackMomentumRes = ROOT.TH1D("trackMomentumRes","Track Momentum Resolution",2000,-10.,10.);
...
                for k in range(0,len(simuAssoc[i])): # Loop over associations to find matching ReconstructedChargedParticle
                    if (simuAssoc[i][k] == j):
                        recMom = ROOT.TVector3(trackMomX[i][recoAssoc[i][k]], trackMomY[i][recoAssoc[i][k]], trackMomZ[i][recoAssoc[i][k]])
                        momRes = (recMom.Mag() - trueMom.Mag())/trueMom.Mag()

                        trackMomentumRes.Fill(momRes)
```

Remember to write this histogram to the output file too! While this plot will give us a sense of what the tracking resolution is, we don't expect the resolution to be constant for all momenta or eta. We can get a more complete picture by plotting the resolution as a function of different kinematic quantities. 

> Exercise:
> - Make 2-D plots of resolution vs true momentum and vs true pseudorapidity.
> - Break resolution plots down by particle species.
{: .challenge}

> Question:
> - Will the histogram ranges for each particle species be the same?
> - Could we present the resolution values in a more understandable way?
{: .callout}

## ROOT RDataFrames

Newer versions of root, such as the version in eic-shell, have access to a relatively new class, [RDataFrames](https://root.cern/doc/master/classROOT_1_1RDataFrame.html). These are similar to pythonic data frame style strucutres that you may be familiar with. Some people are moving towards utilising RDataFrames in their analysis. If you are more familiar with working with data frames, you may wish to investigate these further.

Included below is a quick script from [Simon Gardner](https://github.com/simonge/EIC_Analysis/blob/main/Analysis-Tutorial/EfficiencyAnalysisRDF.C) that utilises RDataFrames to analyse a data file. Copy the following into a new file called `EfficiencyAnalysisRDF.C` -

```c++
// Define the function to perform the efficiency analysis
void EfficiencyAnalysisRDF(TString infile="PATH_TO_FILE"){
   
  // Set up input file 
  ROOT::RDataFrame df("events", infile);

  // Define new dataframe node with additional columns
  auto df1 =  df.Define("statusFilter",  "MCParticles.generatorStatus == 1"    )
                .Define("absPDG",        "abs(MCParticles.PDG)"                )
                .Define("pdgFilter",     "absPDG == 11 || absPDG == 13 || absPDG == 211 || absPDG == 321 || absPDG == 2212")
                .Define("particleFilter","statusFilter && pdgFilter"           )
                .Define("filtMCParts",   "MCParticles[particleFilter]"         )
                .Define("assoFilter",    "Take(particleFilter,ReconstructedChargedParticleAssociations.simID)") // Incase any of the associated particles happen to not be charged
                .Define("assoMCParts",   "Take(MCParticles,ReconstructedChargedParticleAssociations.simID)[assoFilter]")
                .Define("assoRecParts",  "Take(ReconstructedChargedParticles,ReconstructedChargedParticleAssociations.recID)[assoFilter]")
                .Define("filtMCEta",     getEta<MCP>   , {"filtMCParts"} )
                .Define("filtMCPhi",     getPhi<MCP>   , {"filtMCParts"} )
                .Define("accoMCEta",     getEta<MCP>   , {"assoMCParts"} )
                .Define("accoMCPhi",     getPhi<MCP>   , {"assoMCParts"} )
                .Define("assoRecEta",    getEta<RecoP> , {"assoRecParts"})
                .Define("assoRecPhi",    getPhi<RecoP> , {"assoRecParts"})
                .Define("deltaR",        "ROOT::VecOps::DeltaR(assoRecEta, accoMCEta, assoRecPhi, accoMCPhi)");

  // Define histograms
  auto partEta                = df1.Histo1D({"partEta","Eta of Thrown Charged Particles;Eta",100,-5.,5.},"filtMCEta");
  auto matchedPartEta         = df1.Histo1D({"matchedPartEta","Eta of Thrown Charged Particles That Have Matching Track",100,-5.,5.},"accoMCEta");
  auto matchedPartTrackDeltaR = df1.Histo1D({"matchedPartTrackDeltaR","Delta R Between Matching Thrown and Reconstructed Charged Particle",5000,0.,5.},"deltaR");

  // Write histograms to file
  TFile *ofile = TFile::Open("EfficiencyAnalysis_Out_RDF.root","RECREATE");

  // Booked Define and Histo1D lazy actions are only performed here
  partEta->Write();
  matchedPartEta->Write();
  matchedPartTrackDeltaR->Write();
      
  ofile->Close(); // Close output file
}
```

> Note:
> You will need to run this script with the command 'root -l -q EfficiencyAnalysisRDF.C++', within eic-shell (or somewhere else with the correct EDM4hep/EDM4eic libraries installed).
> Remember to put in the correct file path.
> {: .callout}

If you like, you can try completing the exercises using this example to start from.
