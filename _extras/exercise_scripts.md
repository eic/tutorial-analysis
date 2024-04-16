---
title: "Exercise Scripts"
---

Included below is a selection of scripts for the exercises in part 3 of this tutorial. Shortly after the tutorial, I will also include "complete" examples for future reference.

You should be able to copy the code text directly into a new file. The name of the file is included as the title of each script section and in the accompanying descriptive text.

## ROOT TTreeReader Scripts

Note that in copying these scripts, the tabbing seemed to get completely screwed up. I *highly* recommend that you select all of the code and hit `tab` once you copy it into a file. Most editors should then correctly tab each line.

### EfficiencyAnalysis.C

Create a file called `EfficiencyAnalysis.C` and copy in the code below to get started on the efficiency analysis exercise. Note that you will need to correctly specifiy your input file path in the first line.

```c++
void EfficiencyAnalysis(TString infile="PATH_TO_INPUT_FILE"){
  // Set output file for the histograms
  TFile *ofile = TFile::Open("EfficiencyAnalysis_Out.root","RECREATE");
  
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
      
  // Define Histograms
  TH1D *partEta = new TH1D("partEta","Eta of Thrown Charged Particles;Eta",100,-5.,5.);
  TH1D *matchedPartEta = new TH1D("matchedPartEta","Eta of Thrown Charged Particles That Have Matching Track",100,-5.,5.);
  TH1D *matchedPartTrackDeltaR = new TH1D("matchedPartTrackDeltaR","Delta R Between Matching Thrown and Reconstructed Charged Particle",5000,0.,5.);

  while(tree_reader.Next()) { // Loop over events
    for(unsigned int i=0; i<partGenStat.GetSize(); i++){ // Loop over thrown particles
	if(partGenStat[i] == 1){ // Select stable thrown particles
	    int pdg = TMath::Abs(partPdg[i]);
	    if(pdg == 11 || pdg == 13 || pdg == 211 || pdg == 321 || pdg == 2212){ // Look at charged particles (electrons, muons, pions, kaons, protons)
		TVector3 trueMom(partMomX[i],partMomY[i],partMomZ[i]);

		float trueEta = trueMom.PseudoRapidity();
		float truePhi = trueMom.Phi();
	    
		partEta->Fill(trueEta);

		for(unsigned int j=0; j<simuAssoc.GetSize(); j++){ // Loop over associations to find matching ReconstructedChargedParticle
		    if(simuAssoc[j] == i){ // Find association index matching the index of the thrown particle we are looking at
			TVector3 recMom(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]); // recoAssoc[j] is the index of the matched ReconstructedChargedParticle

			// Check the distance between the thrown and reconstructed particle
			float deltaEta = trueEta - recMom.PseudoRapidity();
			float deltaPhi = TVector2::Phi_mpi_pi(truePhi - recMom.Phi());
			float deltaR = TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);

			matchedPartTrackDeltaR->Fill(deltaR);

			matchedPartEta->Fill(trueEta); // Plot the thrown eta if a matched ReconstructedChargedParticle was found
                    }
                } // End loop over associations 
            } // End PDG check          
        } // End stable particles condition  
    } // End loop over thrown particles
  } // End loop over events 
  ofile->Write(); // Write histograms to file
  ofile->Close(); // Close output file
}
```

### ResolutionAnalysis.C

Create a file called `ResolutionAnalysis.C` and copy in the code below to get started on the resolution analysis exercise. Note that you will need to correctly specifiy your input file path in the first line.

```c++
void ResolutionAnalysis(TString infile="PATH_TO_INPUT_FILE"){
  // Set output file for the histograms
  TFile *ofile = TFile::Open("ResolutionAnalysis_Out.root","RECREATE");

  // Analysis code will go here
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
    
  // Define Histograms
  TH1D *trackMomentumRes = new TH1D("trackMomentumRes","Track Momentum Resolution", 400, -2, 2);
 
  TH1D *matchedPartTrackDeltaEta = new TH1D("matchedPartTrackDeltaEta","#Delta#eta Between Matching Thrown and Reconstructed Charged Particle; #Delta#eta", 100, -0.25, 0.25);
  TH1D *matchedPartTrackDeltaPhi = new TH1D("matchedPartTrackDeltaPhi","#Detla #phi Between Matching Thrown and Reconstructed Charged Particle; #Delta#phi", 200, -0.2, 0.2);
  TH1D *matchedPartTrackDeltaR = new TH1D("matchedPartTrackDeltaR","#Delta R Between Matching Thrown and Reconstructed Charged Particle; #Delta R", 300, 0, 0.3);
  TH1D *matchedPartTrackDeltaMom = new TH1D("matchedPartTrackDeltaMom","#Delta P Between Matching Thrown and Reconstructed Charged Particle; #Delta P", 200, -10, 10);
  while(tree_reader.Next()) { // Loop over events
    for(unsigned int i=0; i<partGenStat.GetSize(); i++){ // Loop over thrown particles
	if(partGenStat[i] == 1){ // Select stable thrown particles
	    int pdg = TMath::Abs(partPdg[i]);
	    if(pdg == 11 || pdg == 13 || pdg == 211 || pdg == 321 || pdg == 2212){ // Look at charged particles (electrons, muons, pions, kaons, protons)
		TVector3 trueMom(partMomX[i],partMomY[i],partMomZ[i]);

		float trueEta = trueMom.PseudoRapidity();
		float truePhi = trueMom.Phi();
	    
		for(unsigned int j=0; j<simuAssoc.GetSize(); j++){ // Loop over associations to find matching ReconstructedChargedParticle
		    if(simuAssoc[j] == i){ // Find association index matching the index of the thrown particle we are looking at
			TVector3 recMom(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]); // recoAssoc[j] is the index of the matched ReconstructedChargedParticle

			// Check the distance between the thrown and reconstructed particle
			float deltaEta = trueEta - recMom.PseudoRapidity();
			float deltaPhi = TVector2::Phi_mpi_pi(truePhi - recMom.Phi());
			float deltaR = TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
			float deltaMom = ((trueMom.Mag()) - (recMom.Mag()));
	    double momRes = (recMom.Mag() - trueMom.Mag())/trueMom.Mag();
	
			trackMomentumRes->Fill(momRes);

			matchedPartTrackDeltaEta->Fill(deltaEta);
			matchedPartTrackDeltaPhi->Fill(deltaPhi);
			matchedPartTrackDeltaR->Fill(deltaR);
			matchedPartTrackDeltaMom->Fill(deltaMom);
                    }
                } // End loop over associations 
            } // End PDG check          
        } // End stable particles condition  
    } // End loop over thrown particles
  } // End loop over events 
  ofile->Write(); // Write histograms to file
  ofile->Close(); // Close output file
}
```
## Python Uproot Scripts

### EfficiencyAnalysis.py

Create a file called `EfficiencyAnalysis.py` and copy in the code below to get started on the resolution analysis exercise. Note that you will need to correctly specifiy your input file path in the variable `infile`.

```python
#! /usr/bin/python

#Import relevant packages
import ROOT, math, array
from ROOT import TCanvas, TColor, TGaxis, TH1F, TH2F, TPad, TStyle, gStyle, gPad, TGaxis, TLine, TMath, TPaveText, TTree, TVector3, TVector2
import uproot as up

#Define and open files
infile="PATH_TO_INPUT_FILE""
ofile=ROOT.TFile.Open("EfficiencyAnalysis_OutPy.root", "RECREATE")

# Open input file and define branches we want to look at with uproot
events_tree = up.open(infile)["events"]

# Get particle information
partGenStat = events_tree.array("MCParticles.generatorStatus")
partMomX = events_tree.array("MCParticles.momentum.x")
partMomY = events_tree.array("MCParticles.momentum.y")
partMomZ = events_tree.array("MCParticles.momentum.z")
partPdg = events_tree.array("MCParticles.PDG")

# Get reconstructed track information
trackMomX = events_tree.array("ReconstructedChargedParticles.momentum.x")
trackMomY = events_tree.array("ReconstructedChargedParticles.momentum.y")
trackMomZ = events_tree.array("ReconstructedChargedParticles.momentum.z")

# Get assocations between MCParticles and ReconstructedChargedParticles
recoAssoc = events_tree.array("ReconstructedChargedParticleAssociations.recID")
simuAssoc = events_tree.array("ReconstructedChargedParticleAssociations.simID")

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

### ResolutionAnalysis.py

Create a file called `ResolutionAnalysis.py` and copy in the code below to get started on the resolution analysis exercise. Note that you will need to correctly specifiy your input file path in the variable `infile`.

```python
#! /usr/bin/python

#Import relevant packages
import ROOT, math, array
from ROOT import TCanvas, TColor, TGaxis, TH1F, TH2F, TPad, TStyle, gStyle, gPad, TGaxis, TLine, TMath, TPaveText, TTree, TVector3, TVector2
import uproot as up

#Define and open files
infile="PATH_TO_INPUT_FILE""
ofile=ROOT.TFile.Open("ResolutionAnalysis_OutPy.root", "RECREATE")

# Open input file and define branches we want to look at with uproot
events_tree = up.open(infile)["events"]

# Get particle information
partGenStat = events_tree.array("MCParticles.generatorStatus")
partMomX = events_tree.array("MCParticles.momentum.x")
partMomY = events_tree.array("MCParticles.momentum.y")
partMomZ = events_tree.array("MCParticles.momentum.z")
partPdg = events_tree.array("MCParticles.PDG")

# Get reconstructed track information
trackMomX = events_tree.array("ReconstructedChargedParticles.momentum.x")
trackMomY = events_tree.array("ReconstructedChargedParticles.momentum.y")
trackMomZ = events_tree.array("ReconstructedChargedParticles.momentum.z")

# Get assocations between MCParticles and ReconstructedChargedParticles
recoAssoc = events_tree.array("ReconstructedChargedParticleAssociations.recID")
simuAssoc = events_tree.array("ReconstructedChargedParticleAssociations.simID")

# Define histograms below
trackMomentumRes = ROOT.TH1D("trackMomentumRes","Track Momentum Resolution", 400, -2, 2)

matchedPartTrackDeltaEta = ROOT.TH1D("matchedPartTrackDeltaEta","#Delta#eta Between Matching Thrown and Reconstructe Charged Particle; #Delta#eta", 100, -0.25, 0.25)
matchedPartTrackDeltaPhi = ROOT.TH1D("matchedPartTrackDeltaPhi","#Detla #phi Between Matching Thrown and Reconstructed Charged Particle; #Delta#phi", 200, -0.2, 0.2)
matchedPartTrackDeltaR = ROOT.TH1D("matchedPartTrackDeltaR","#Delta R Between Matching Thrown and Reconstructed Charged Particle; #Delta R", 300, 0, 0.3)
matchedPartTrackDeltaMom = ROOT.TH1D("matchedPartTrackDeltaMom","#Delta P Between Matching Thrown and Reconstructed Charged Particle; #Delta P", 200, -10, 10)

# Add main analysis loop(s) below
for i in range(0, len(events_tree)): # Loop over all events
    for j in range(0, len(partGenStat[i])): # Loop over all thrown particles
        if partGenStat[i][j] == 1: # Select stable particles
            pdg = abs(partPdg[i][j]) # Get PDG for each stable particle
            if(pdg == 11 or pdg == 13 or pdg == 211 or pdg == 321 or pdg == 2212):
                trueMom = ROOT.TVector3(partMomX[i][j], partMomY[i][j], partMomZ[i][j])
                trueEta = trueMom.PseudoRapidity()
                truePhi = trueMom.Phi()
                
                for k in range(0,len(simuAssoc[i])): # Loop over associations to find matching ReconstructedChargedParticle
                    if (simuAssoc[i][k] == j):
                        recMom = ROOT.TVector3(trackMomX[i][recoAssoc[i][k]], trackMomY[i][recoAssoc[i][k]], trackMomZ[i][recoAssoc[i][k]])
                        deltaEta = trueEta - recMom.PseudoRapidity()
                        deltaPhi = TVector2. Phi_mpi_pi(truePhi - recMom.Phi())
                        deltaR = math.sqrt((deltaEta*deltaEta) + (deltaPhi*deltaPhi))
                        deltaMom = ((trueMom.Mag()) - (recMom.Mag()))

                        momRes = (recMom.Mag() - trueMom.Mag())/trueMom.Mag()
                        
                        matchedPartTrackDeltaEta.Fill(deltaEta)
                        matchedPartTrackDeltaPhi.Fill(deltaPhi)
                        matchedPartTrackDeltaR.Fill(deltaR)
                        matchedPartTrackDeltaMom.Fill(deltaMom)                        

                        trackMomentumRes.Fill(momRes)
                        
# Write output histograms to file below
trackMomentumRes.Write()
matchedPartTrackDeltaEta.Write()
matchedPartTrackDeltaPhi.Write()
matchedPartTrackDeltaR.Write()
matchedPartTrackDeltaMom.Write()

# Close files
ofile.Close()

```

{% include links.md %}
