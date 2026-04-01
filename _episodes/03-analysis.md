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

> Note:
> - The [branch dictionary]({{ page.root }}{% link _extras/branch_dictionary.md %}) outlines all of the branches we will need to utilise in this section.
> - If you want, you can prune the branches you don't need from the input file using the [TreePrune.C script]({{ page.root }}{% link _extras/tree_pruning_script.md %}).
{: .callout}

## Reading the Output Trees

The simulation output trees are "flat" in the sense that there is no event class structure embedded within the tree and no additional libraries are needed to handle the output. Therefore, the end user can simply read the values stored in each branch using whatever method/workflow they are most comfortable with. Examples of several common methods for reading the trees are provided below. We will see a ROOT TTreeReader based example using a ROOT macro and a python/uproot based version. There is also an example using the (relatively) new RDataFrame class of ROOT. During the tutorial, you should try the exercise using whichever language you feel most comfortable with. Five different approaches are currently provieded:

- [TTreeReaders](https://eic.github.io/tutorial-analysis/03-analysis/index.html#sample-analysis-with-root-ttreereader-track-efficiency-and-resolution) - ROOT/C based
- [Python/Uproot - Pythonic](https://eic.github.io/tutorial-analysis/03-analysis/index.html#sample-analysis-with-pythonuproot---pythonic-method-track-efficiency-and-resolution) - A pythonic based appoach using arrays directly
- [Python/Uproot - ROOT/Pyroot](https://eic.github.io/tutorial-analysis/03-analysis/index.html#sample-analysis-with-pythonuproot---rootpyroot-style-track-efficiency-and-resolution) - An approach using Pyroot, halfway house between C and python
- [ROOT RDataFrames](https://eic.github.io/tutorial-analysis/03-analysis/index.html#root-rdataframes) - An approach using RDataFrames
- [PODIO](https://eic.github.io/tutorial-analysis/03-analysis/index.html#podio---direct-analysis) - An approach using the Plane Old Data IO (PODIO) approach. Use the flat datastructure directly

## Sample Analysis with ROOT TTreeReader: Track Efficiency and Resolution

As a sample exercise to become familiar with the simulation output and how to use it in a realistic setting, we will find the tracking eficiency and resolution. We will need to access the reconstructed track information and the truth particle information and we will have to associate the individual tracks and particles to one another. 

Before we begin, we should create a skeleton macro to handle file I/O. For the `TTreeReader` example, we will use a simple ROOT macro. Using your favorite text editor, create a file with a name like `trackAnalysis.C` or something similar and copy in the following code: 

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
> mychain->Add(infile);
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
TTreeReaderArray<double> partMomX(tree_reader, "MCParticles.momentum.x");
TTreeReaderArray<double> partMomY(tree_reader, "MCParticles.momentum.y");
TTreeReaderArray<double> partMomZ(tree_reader, "MCParticles.momentum.z");
TTreeReaderArray<int> partPdg(tree_reader, "MCParticles.PDG");

// Get Reconstructed Track Information
TTreeReaderArray<float> trackMomX(tree_reader, "ReconstructedChargedParticles.momentum.x");
TTreeReaderArray<float> trackMomY(tree_reader, "ReconstructedChargedParticles.momentum.y");
TTreeReaderArray<float> trackMomZ(tree_reader, "ReconstructedChargedParticles.momentum.z");

// Get Associations Between MCParticles and ReconstructedChargedParticles
TTreeReaderArray<int> recoAssoc(tree_reader, "_ReconstructedChargedParticleAssociations_rec.index");
TTreeReaderArray<int> simuAssoc(tree_reader, "_ReconstructedChargedParticleAssociations_sim.index");
```

The last two lines encode the association between a ReconstructedChargedParticle and an MCParticle where the matching is determined in the [ParticlesWithPID](https://github.com/eic/EICrecon/blob/main/src/algorithms/pid/ParticlesWithPID.cc) algorithm which generates the ReconstructedChargedParticle objects.

> Compiling ROOT Macros:
> - If you are analysing a large number of events, you may wish to compile your macro to increase throughput. An example of how you can create and compile a root macro is included in the [Exercise Scripts section](https://eic.github.io/tutorial-analysis/exercise_scripts/index.html#compiled-root-scripts)
{: .callout}

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
> For all **scattered electrons**, **charged pions** and **protons** in our events:
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
> For all **scattered electrons**, **charged pions** and **protons** in our events:
> - Make 2-D plots of resolution vs true momentum and vs true pseudorapidity.
{: .challenge}

> Question:
> - Will the histogram ranges for each particle species be the same?
> - Could we present the resolution values in a more understandable way?
{: .callout}


## Sample Analysis with Python/uproot - Pythonic Method: Track Efficiency and Resolution

> For some examples of using uproot to access information in .root files, please consult [this notebook](https://github.com/eic/HSF-India/blob/main/Working_With_Uproot/Working_With_Uproot_Standalone.ipynb) which can be run in Google Collab.
{: .callout}

If you are more familiar with python than you are with C/C++, you might find that using a python based root macro is easier for you. Outlined below are sample blocks of code for creating and running a python based analysis script.

With python, some tasks become easier, e.g. string manipulation and writing to (non ROOT) files.

Before we begin, we should create a skeleton macro to handle file I/O. For this example, we will make a simple python script. Using your favourite editor, create a file with a name like `trackAnalysis.py` or something similar and copy in the following code: 

```python
#! /usr/bin/python
# Import some relevant packages
import uproot as up
import awkward as ak
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.ticker as ticker
import matplotlib.cm as cm
import matplotlib.pylab as plt
import scipy, vector, os
from XRootD import client
from scipy import stats
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import colors as colours

# Set some matplot lib features
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['xaxis.labellocation'] = 'right'
plt.rcParams['yaxis.labellocation'] = 'top'
plt.rcParams["figure.figsize"] = (16,9)
kP6 = ['#5790fc','#f89c20','#e42536','#964a8b','#9c9ca1','#7a21dd'] # Set ROOT kP6 colours - see https://root.cern.ch/doc/v636/classTColor.html

# Open our file
fname = "INPUT_FILE.root"
if os.path.isfile(fname):
    file=up.open(fname)
else:
    print("Error opening file - ", fname, " check your fname variable!")

# Open the tree
tree = file['events']

# Convert relevant branches to arrays
MCPartBr = tree["MCParticles"].arrays()

# Define some filters

# Use filters to manipulate data

# Create and write some plots
```

Make sure you change the input file name to match whatever you saved your file as earlier.

Note that we are using the module uproot to access the data here. See [further documentation here](https://masonproffitt.github.io/uproot-tutorial/03-trees/index.html). You may also need some of the other included packages too.

We will use uproot a little bit like we use the TTreeReader in the other example. We can define the branches we want and assign them to arrays with uproot.
We can do this via:
```python
# Open input file and define branches we want to look at with uproot
fname = "INPUT_FILE.root" # Your file name
file=up.open(fname)
tree = file['events']
# Convert relevant branches to arrays
MCPartBr = tree["MCParticles"].arrays()
# If we want, convert a specific branch element to an array
partPdg = tree["MCParticles.PDG"].array()
```
Uproot effectively takes the information in the tree, and turns it into an array. We can then access and manipulate this array in the same way that we can with any array in python.

> Warning: Note that if you are using an older version of uproot (v2.x.x), you will need to access the branches slightly differently via -
> ```python
> partGenStat = tree.array("MCParticles.generatorStatus")
> ```
{: .callout}
Once you create your script and add the template code in, you can try running it with``python3 trackAnalysis.py`` or ``python trackAnalysis.py``. At the moment, it shouldn't *do* anything, but we can change that!

Try assigning a quantity to an array, such as the MC particles PDG values above and printing some values of that array. Or, perhaps try printing the length of that array. We could also quickly make a plot of the values with -

```python
plt.hist(ak.flatten(partPdg),alpha=0.75, color=kP6[0])
plt.savefig("TestOut.png", dpi = (160))
plt.clf # Clear figure
```
The script should now write out a figure, ``TestOut.png`` when you run it, showing the MC PDG values of entries in the file.

> We did not specify a number of bins or a range, so our plot looks a bit odd. What might be a useful range and number of bins to use here?
> We can specify our number of bins and our range with -
> ```python
> plt.hist(ak.flatten(partPdg),bins=NBins, range=(X,Y),alpha=0.75, color=kP6[0])
> ```
> With NBins being our number of bins and X and Y being our min/max range - think carefully about these numbers!
{: .callout}

Note that we don't really need to define individual arrays either, we can just directly access the array we want once we've converted the branch to a series of arrays -

```python
MCPartBr = tree["MCParticles"].arrays()
plt.hist(ak.flatten(MCPartBr["MCParticles.PDG"]),alpha=0.75, color=kP6[0])
plt.savefig("TestOut.png", dpi = (160))
```

> Note:
> Remember to call:
> ```python
> plt.clf() # Clear figure
> ```
> After a figure to avoid drawing on the same plot.
{: .callout}

We can also define and apply filters to our arrays as we plot or print from them -

```python
# We will filter on the particle status. Generator status == 1 corresponds to a stable particle (as opposed to a beam or intermediate particle) that we could detect in our detector. 4 is for beam particles
BoolStable=(MCPartBr["MCParticles.generatorStatus"]==1) # This filter is actually an array of booleans. Any where the generator status for a particle == 1 will return true
plt.hist(ak.flatten(MCPartBr["MCParticles.PDG"][BoolStable]),alpha=0.75, color=kP6[0]) # Apply filter to PDG array as we plot it. Only stable particles will now be plotted
plt.savefig("TestOut2.png", dpi = (160))
```

And we can also combine conditions in our filters -

```python
BoolStablePos=((MCPartBr["MCParticles.generatorStatus"]==1) & (MCPartBr["MCParticles.charge"]>0)) # Create a filter to select out stable, positively chrarged particles 
plt.hist(ak.flatten(MCPartBr["MCParticles.PDG"][BoolStablePos]),alpha=0.75, color=kP6[0]) # Apply filter to PDG array as we plot it. Only stable particles will now be plotted
plt.savefig("TestOut3.png", dpi = (160))
```

> We did not specify a number of bins or a range, so our plot looks a bit odd. What might be a useful range and number of bins to use here?
{: .callout}

### Efficiency Analysis

> Hint:
> Refer to [the script template](https://eic.github.io/tutorial-analysis/exercise_scripts/index.html#pythonic_efficiencyanalysispy) if you're having trouble putting things in the right place.
{: .callout}

Our approach here is a bit different to the TTreeReader example, but we will still need to utilise our simulation and reconstruction association IDs. We can access them via -

```python
RecoAssocRec = tree['_ReconstructedChargedParticleAssociations_rec'].arrays()
RecoAssocSim = tree['_ReconstructedChargedParticleAssociations_sim'].arrays()
RecID=RecoAssocRec['_ReconstructedChargedParticleAssociations_rec.index'] # Array of reconstructed IDs
SimID=RecoAssocSim['_ReconstructedChargedParticleAssociations_sim.index'] # Array of simulated IDs
```

These are arrays of indices which map our simulated data to events which have actually been reconstructed in our simulation. We can use these to index our arrays and pick out *only* events where a simulated particle has a matching reconstructed particle (OR vice versa, where a reconstructed particle has a matching simulated particle). Note that when we index by these particles, we also need to make sure any filters are also indexed appropriately. E.g.

```python
BoolChargeTrack = ((abs(MCPartBr["MCParticles.charge"])!=0) & (MCPartBr["MCParticles.generatorStatus"]==1)) # A filter to select out charged stable particles in the MC data
```

If we apply this filter to an array which has been indexed by the Simulation ID, we will run into problems -

```python
BoolChargeTrack = ((abs(MCPartBr["MCParticles.charge"])!=0) & (MCPartBr["MCParticles.generatorStatus"]==1)) # A filter to select out charged stable particles in the MC data
MatchPDG = MCPartBr["MCParticles.PDG"][SimID] # An array of the MC PDG values for all MC particles with a matching reconstructed particle.
# print(MatchPDG[BoolChargeTrack]) # Will return an error, arrays don't match sizing!
```

We could either -

1. Create a new filter which is explicitly indexed by the SimID's
2. Index our previous filter by the SimID as we use it

Both should return the same answer -

```python
print(MatchPDG[BoolChargeTrack[SimID]]) # Explicitly index our filter by SimID first
BoolChargeTrackMatch = ((abs(MCPartBr["MCParticles.charge"][SimID])!=0) & (MCPartBr["MCParticles.generatorStatus"][SimID]==1)) # A filter to select out charged stable particles in the MC data
print(MatchPDG[BoolChargeTrackMatch]) # Use a newly defined filter which has been indexed by the SimID
```
To select out reconstructed particles that have a matching simulated particle, we can index by ``RecID`` in a similar way. **Note that we need to be careful what we conclude from our analysis if we match particles in this way. We are "cheating" in the sense that we are directly matching the truth to what we reconstruct. We cannot do this in a real experiment.**

We may want to group some of our quantities together as vectors for easy manipulation. For example, we could create an array of vectors corresponding to our charged particles -

```python
MC_Parts = vector.zip({'px': MCPartBr["MCParticles.momentum.x"], 'py': MCPartBr["MCParticles.momentum.y"], 'pz': MCPartBr["MCParticles.momentum.z"]})
```
We can filter and index this like any other array. Creating 3 or 4-vectors in this way is useful as we can then use various functions to extract information from our vectors. For example, we can easily get -

  - *Pseudorapidity* - Eta
  - *Polar angle* - Theta they make wrt the origin (in the lab frame, our bunch crossing point) - *Note, this is in radians by default*
  - *Transverse momentum - PT
  - ...

```python
print(MC_Parts.eta)
print(MC_Parts.theta)
print(MC_Parts.pt)
```

Now, we can use what we know to determine and plot some simple efficiency graphs for our reconstructed data. Efficiency is a measure of the probability that we will detect an incident particle. We could calculate our efficiency by straightforwardly counting how many particles we detect vs how many we "threw". For example, if we detect 9 particles and 10 were generated, our efficiency is -

9/10

i.e.

90%

This might be a useful figure. However, if we are evaluating the performance of a detector, it might be useful to consider the efficiency as a function of another quantity, for example *eta* or *p*. What this will tell us is how likely we are to detect particles incident on certain areas of the detector (or with a certain momentum for instance). We should not really expect these distributions to be completely flat.

In terms of our code, we can straightforwardly determine this by dividing some histograms. We can divide histograms via -

```python
RecPartBr = tree["ReconstructedChargedParticles"].arrays()
BoolElec=((MCPartBr["MCParticles.PDG"][SimID]==11) & (MCPartBr["MCParticles.generatorStatus"][SimID]==1)) # Define a filter to select out electrons which reconstruct in our detector
MCHist = np.histogram(ak.flatten(MCPartBr["MCParticles.momentum.x"][BoolElec]), bins=100, range=(0,25)) # Create a hisogram of x-momenta for MC particles that are electrons
RecHist = np.histogram(ak.flatten(RecPartBr["ReconstructedChargedParticles.momentum.x"][RecID][BoolElec]), bins=100, range=(0,25)) # Create a histogram of reconstructed x-momenta values for particles that are actual MC electrons that have reconstructed
with np.errstate(divide='ignore'):
    Division = RecHist[0] / MCHist[0]
Division = np.nan_to_num(Division,nan=0, posinf = 0) # Convert any nan or pos inf from /0 (empty bins) to 0
Bin_Edges=MCHist[1]
Bars = 0.5 * (Bin_Edges[1:] + Bin_Edges[:-1])
BarWidth=Bars[1]-Bars[0]
plt.bar(Bars, Division, width=BarWidth, alpha=0.75, color=kP6[2])
plt.savefig("TestOut4.png", dpi = (160))
```

**Important** - What we have create here is __not__ our efficiency! We have simply divided the **reconstructed** electron P_{X} by its true value for MC electrons that have reconstructed in our detector. We used the PDG code for electrons, 11, to pick out electrons from our MC particles branch. There are a few caveats to actually calculating our efficiency. This just demonstrates how we can divide histograms in python.

For our efficiency. We need to compare our thrown particles of a given type to our detected particles of a given type.

- When we do our division, we should do this for the same quantity in each case (i.e. compare the true and values to each other).
  - How can you select the thrown MC Particles of a specific type?
  - How can you select the particles we detected of a specific type?
    - Note, this does not mean we need our reconstructed values.

> Exercise:
> For all **scattered electrons**, **charged pions** and **protons** in our events:
> - Find the efficiency as a function of particle momentum. Are there cuts on any other quantities you should place to get a sensible result?
> - Find the efficiency for some 2-D correlations: momentum vs eta; phi vs eta
> - Plot some kinematic distributions (momentum, eta, etc) for all ReconstructedChargedParticles, not just those that are associated with a thrown particle
{: .challenge}

> Hint:
> Getting the right arrays here is a bit tricky. We want three different things -
> - Our MC particles (truth information), regardless of whether we have a matching track or not. This is just:
>   - "MCPartBr['MCParticles.QUANTITY']"[SELECTION_CUTS] - We do not need to index this by the SimID
> - Our MC particles (truth information) that do have a matching reconstructed track, we just need to index these by the SimID:
>   - "MCPartBr['MCParticles.QUANTITY'][SimID]" - We can then apply selection criteria
> - The Reconstructed particle information for events which correspond to a real MC track, we just need to index these by our RecID:
>   - "ReconChPartBr['ReconstructedChargedParticles.QUANTITY'][RecID]"
{: .callout}

> 2D Histograms: We can make 2D histograms in python via -
> ```python
> plt.hist2d(np.asarray(ak.flatten(Quantity1)), np.asarray(ak.flatten(Quantity2)), bins=[NBinsX,NBinsY], range=[[XLow,XHigh],[YLow,YHigh]], cmin=1)
> cb = plt.colorbar()
> cb.set_label('Counts/bin')
> ```
> We can set titles etc as usual. Simply swap on the bin values and ranges, as well as the quantities as needed. Make sure your arrays contain equal numbers of entries.
{: .callout}

### Resolution Analysis

> Hint:
> Refer to [the script template](https://eic.github.io/tutorial-analysis/exercise_scripts/index.html#pythonic_resolutionanalysispy) if you're having trouble putting things in the right place.
{: .callout}

Next, we will look at track momentum resolution. The resolution tells us how well we can reconstruct our "true" value. For example. we might want to know how well we can determine the energy of our particles. As such, we could calculate the energy resolution. Our resolution is simply -

- (Reconstructed - True)/True

This is often expressed as a percentage. So, for example, say we detect a particle and determine its energy to be 0.95 GeV. In reality, the energy was 1 GeV. As such, our energy resolution for this particle is -

- 0.95-1/1 = -5%

Let's see a quick example of this calculation for a a quantity -

```python
ElecMomXRes = ((RecPartBr["ReconstructedChargedParticles.momentum.x"][RecID][BoolElec]-MCPartBr["MCParticles.momentum.x"][SimID][BoolElec])/MCPartBr["MCParticles.momentum.x"][SimID][BoolElec])*100 # Multiply by 100 to get this as a %
plt.hist(ak.flatten(ElecMomXRes), bins=50, range=(-25,25),alpha=0.5, color=kP6[2])
plt.savefig("TestOut5.png", dpi = (160))
```

Here we've calculated and plotted the X momentum resolution for our charged tracks that correspond to true electrons in our sample. Whilst this plot will give us a sense of what the tracking resolution is, we don't expect the resolution to be constant for all momenta or eta. We can get a more complete picture by plotting the resolution as a function of different kinematic quantities. 

> Exercise:
> For all **scattered electrons**, **charged pions** and **protons** in our events:
> - Make 2-D plots of resolution vs true momentum and vs true pseudorapidity.
{: .challenge}

> Question:
> - Will the histogram ranges for each particle species be the same?
> - Could we present the resolution values in a more understandable way?
{: .callout}

## Sample Analysis with Python/uproot - ROOT/Pyroot Style: Track Efficiency and Resolution

> Comment:
> Despite using python/uproot, I have written these in a very "ROOT"/C way. Uproot converts our branches to arrays, so you can manipulate them in various fun ways using more pythonic methods if you want.
> See the previous method for an example of this approach.
{: .callout}

If you are more familiar with python than you are with C/C++, you might find that using a python based root macro is easier for you. Outlined below are sample blocks of code for creating and running a python based analysis script.

With python, some tasks become easier, e.g. string manipulation and writing to (non ROOT) files.

Before we begin, we should create a skeleton macro to handle file I/O. For this example, we will make a simple python script. Using your favourite editor, create a file with a name like `trackAnalysis.py` or something similar and copy in the following code: 

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
> for i in range(0, len(partGenStat)): # Loop over all events
>     for j in range(0, len(partGenStat[i])): # Loop over all thrown particles
>         if partGenStat[i][j] == 1: # Select stable particles
>             pdg = abs(partPdg[i][j]) # Get PDG for each stable particle
>             ...
> ```
> Uproot effectively takes the information in the tree, and turns it into an array. We can then access and manipulate this array in the same way that we can with any array in python.
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
recoAssoc = events_tree["_ReconstructedChargedParticleAssociations_rec.index"].array()
simuAssoc = events_tree["_ReconstructedChargedParticleAssociations_sim.index"].array()

# Define histograms below
partEta = ROOT.TH1D("partEta","Eta of Thrown Charged Particles;Eta",100, -5 ,5 )
matchedPartEta = ROOT.TH1D("matchedPartEta","Eta of Thrown Charged Particles That Have Matching Track", 100, -5 ,5);
matchedPartTrackDeltaR = ROOT.TH1D("matchedPartTrackDeltaR","Delta R Between Matching Thrown and Reconstructed Charge Particle", 5000, 0, 5);

# Add main analysis loop(s) below
for i in range(0, len(partGenStat)): # Loop over all events
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
> For all **scattered electrons**, **charged pions** and **protons** in our events:
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
> For all **scattered electrons**, **charged pions** and **protons** in our events:
> - Make 2-D plots of resolution vs true momentum and vs true pseudorapidity.
{: .challenge}

> Question:
> - Will the histogram ranges for each particle species be the same?
> - Could we present the resolution values in a more understandable way?
{: .callout}

## ROOT RDataFrames

> Note:
> - This method does actually need you to be within eic-shell (or somewhere else with the correct EDM4hep/EDM4eic libraries installed).
{: .callout}

Newer versions of root, such as the version in eic-shell, have access to a relatively new class, [RDataFrames](https://root.cern/doc/master/classROOT_1_1RDataFrame.html). These are similar to pythonic data frame style structures that you may be familiar with. Some people are moving towards utilising RDataFrames in their analysis. If you are more familiar with working with data frames, you may wish to investigate these further.

Included below is a quick script from [Simon Gardner](https://github.com/simonge/EIC_Analysis/blob/main/Analysis-Tutorial/EfficiencyAnalysisRDF.C) that utilises RDataFrames to analyse a data file. Copy the following into a new file called `EfficiencyAnalysisRDF.C` -

```c++
#include <edm4hep/utils/vector_utils.h>
#include <edm4hep/MCParticle.h>
#include <edm4eic/ReconstructedParticle.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TFile.h>

// Define aliases for the data types 
using MCP = edm4hep::MCParticleData;
using RecoP = edm4eic::ReconstructedParticleData;

// Define function to vectorize the edm4hep::utils methods
template <typename T>
auto getEta = [](ROOT::VecOps::RVec<T> momenta) {
  return ROOT::VecOps::Map(momenta, [](const T& p) { return edm4hep::utils::eta(p.momentum); });
};

template <typename T>
auto getPhi = [](ROOT::VecOps::RVec<T> momenta) {
  return ROOT::VecOps::Map(momenta, [](const T& p) { return edm4hep::utils::angleAzimuthal(p.momentum); });
};

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
                .Define("assoFilter",    "Take(particleFilter,_ReconstructedChargedParticleAssociations_simID.index)") // Incase any of the associated particles happen to not be charged
                .Define("assoMCParts",   "Take(MCParticles,_ReconstructedChargedParticleAssociations_simID.index)[assoFilter]")
                .Define("assoRecParts",  "Take(ReconstructedChargedParticles,_ReconstructedChargedParticleAssociations_recID.index)[assoFilter]")
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
> - You will need to run this script with the command `root -l -q EfficiencyAnalysisRDF.C++`, within eic-shell (or somewhere else with the correct EDM4hep/EDM4eic libraries installed).
> - Remember to put in the correct file path.
{: .callout}

If you like, you can try completing the exercises using this example to start from.

## PODIO - Direct Analysis

If you want to avoid ROOT entirely, you can analyse the PODIO files directly in a variety of ways.

See [Wouter's example use cases](https://indico.cern.ch/event/1343984/contributions/5908856/attachments/2842958/4970156/2024-04-23%20-%20Examples%20for%20Data%20Model%20Usage.pdf) from 23/04/24. Wouter shows a few ways in which the PODIO file can be accessed and analysed directly.

**As of March 2026, a full example and version of this method will be provided in the near future.**
