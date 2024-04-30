---
title: "Exercise Scripts"
---

Included below is a selection of scripts for the exercises in part 3 of this tutorial. Shortly after the tutorial, I will also include "complete" examples for future reference.

You should be able to copy the code text directly into a new file. The name of the file is included as the title of each script section and in the accompanying descriptive text.

## ROOT TTreeReader Scripts

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

A "solution" version of the script for the exercise is included below -

```c++
void EfficiencyAnalysis_Exercise(TString infile="PATH_TO_FILE"){
  
  // Set output file for the histograms
  TFile *ofile = TFile::Open("EfficiencyAnalysis_Exercise_Out.root","RECREATE");

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
  TH1D *partEta = new TH1D("partEta","#eta of Thrown Charged Particles; #eta", 120, -6, 6);
  TH1D *matchedPartEta = new TH1D("matchedPartEta","#eta of Thrown Charged Particles That Have Matching Track; #eta", 120, -6, 6);
  TH1D* partMom = new TH1D("partMom", "Momentum of Thrown Charged Particles (truth); P(GeV/c)", 150, 0, 150);
  TH1D* matchedPartMom = new TH1D("matchedPartMom", "Momentum of Thrown Charged Particles (truth), with matching track; P(GeV/c)", 150, 0, 150);
  TH1D* partPhi = new TH1D("partPhi", "#phi of Thrown Charged Particles (truth); #phi(rad)", 320, -3.2, 3.2);
  TH1D* matchedPartPhi = new TH1D("matchedPartPhi", "#phi of Thrown Charged Particles (truth), with matching track; #phi(rad)", 320, -3.2, 3.2);

  TH2D* partPEta = new TH2D("partPEta", "P vs #eta of Thrown Charged Particles; P(GeV/c); #eta", 150, 0, 150, 120, -6, 6);
  TH2D* matchedPartPEta = new TH2D("matchedPartPEta", "P vs #eta of Thrown Charged Particles, with matching track; P(GeV/c); #eta", 150, 0, 150, 120, -6, 6);
  TH2D* partPhiEta = new TH2D("partPhiEta", "#phi vs #eta of Thrown Charged Particles; #phi(rad); #eta", 160, -3.2, 3.2, 120, -6, 6);
  TH2D* matchedPartPhiEta = new TH2D("matchedPartPhiEta", "#phi vs #eta of Thrown Charged Particles; #phi(rad); #eta", 160, -3.2, 3.2, 120, -6, 6);
    
  TH1D *matchedPartTrackDeltaEta = new TH1D("matchedPartTrackDeltaEta","#Delta#eta Between Matching Thrown and Reconstructed Charged Particle; #Delta#eta", 100, -0.25, 0.25);
  TH1D *matchedPartTrackDeltaPhi = new TH1D("matchedPartTrackDeltaPhi","#Detla #phi Between Matching Thrown and Reconstructed Charged Particle; #Delta#phi", 200, -0.2, 0.2);
  TH1D *matchedPartTrackDeltaR = new TH1D("matchedPartTrackDeltaR","#Delta R Between Matching Thrown and Reconstructed Charged Particle; #Delta R", 300, 0, 0.3);
  TH1D *matchedPartTrackDeltaMom = new TH1D("matchedPartTrackDeltaMom","#Delta P Between Matching Thrown and Reconstructed Charged Particle; #Delta P", 200, -10, 10);
    
  // Define some histograms for our efficiencies
  TH1D *TrackEff_Eta = new TH1D("TrackEff_Eta", "Tracking efficiency as fn of #eta; #eta; Eff(%)", 120, -6, 6); 
  TH1D *TrackEff_Mom = new TH1D("TrackEff_Mom", "Tracking efficiency as fn of P; P(GeV/c); Eff(%)", 150, 0, 150); 
  TH1D *TrackEff_Phi = new TH1D("TrackEff_Phi", "Tracking efficiency as fn of #phi; #phi(rad); Eff(%)", 320, -3.2, 3.2);

  // 2D Efficiencies
  TH2D* TrackEff_PEta = new TH2D("TrackEff_PEta", "Tracking efficiency as fn of P and #eta; P(GeV/c); #eta", 150, 0, 150, 120, -6, 6);
  TH2D* TrackEff_PhiEta = new TH2D("TrackEff_PhiEta", "Tracking efficiency as fn of #phi and #eta; #phi(rad); #eta", 160, -3.2, 3.2, 120, -6, 6);

  // All charged particle histos
  TH1D *ChargedEta = new TH1D("ChargedEta", "#eta of all charged particles; #eta", 120, -6, 6);
  TH1D *ChargedPhi = new TH1D("ChargedPhi", "#phi of all charged particles; #phi (rad)", 120, -3.2, 3.2);
  TH1D *ChargedP = new TH1D("ChargedP", "P of all charged particles; P(GeV/c)", 150, 0, 150);
  
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
		partPhi->Fill(truePhi);
		partMom->Fill(trueMom.Mag());
		partPEta->Fill(trueMom.Mag(), trueEta);
		partPhiEta->Fill(truePhi, trueEta);

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
			float deltaMom = ((trueMom.Mag()) - (recMom.Mag()));

			matchedPartTrackDeltaEta->Fill(deltaEta);
			matchedPartTrackDeltaPhi->Fill(deltaPhi);
			matchedPartTrackDeltaR->Fill(deltaR);
			matchedPartTrackDeltaMom->Fill(deltaMom);

			matchedPartEta->Fill(trueEta); // Plot the thrown eta if a matched ReconstructedChargedParticle was found
			matchedPartPhi->Fill(truePhi);
			matchedPartMom->Fill(trueMom.Mag());

			matchedPartPEta->Fill(trueMom.Mag(), trueEta);
			matchedPartPhiEta->Fill(truePhi, trueEta);
	
		      }
		  }// End loop over associations
	      } // End PDG check
	  } // End stable particles condition
      } // End loop over thrown particles
    // Loop over all charged particles and fill some histograms of kinematics quantities
    for(unsigned int k=0; k<trackMomX.GetSize(); k++){ // Loop over all charged particles, thrown or not
      
      TVector3 CPartMom(trackMomX[k], trackMomY[k], trackMomZ[k]);

      float CPartEta = CPartMom.PseudoRapidity();
      float CPartPhi = CPartMom.Phi();

      ChargedEta->Fill(CPartEta);
      ChargedPhi->Fill(CPartPhi);
      ChargedP->Fill(CPartMom.Mag());
      
    } // End loop over all charged particles
  } // End loop over events

  // Take the ratio of the histograms above to get our efficiency plots
  TrackEff_Eta->Divide(matchedPartEta, partEta, 1, 1, "b");
  TrackEff_Mom->Divide(matchedPartMom, partMom, 1, 1, "b");
  TrackEff_Phi->Divide(matchedPartPhi, partPhi, 1, 1, "b");
  TrackEff_PEta->Divide(matchedPartPEta, partPEta, 1, 1, "b");
  TrackEff_PhiEta->Divide(matchedPartPhiEta, partPhiEta, 1, 1, "b");
  
  ofile->Write(); // Write histograms to file
  ofile->Close(); // Close output file
}
```
Insert your input file path and execute as the example code above.

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
			double momRes = (recMom.Mag()- trueMom.Mag())/trueMom.Mag();
      
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

A "solution" version of the script for the exercise is included below -

```c++
void ResolutionAnalysis_Exercise(TString infile="PATH_TO_FILE"){
  // Set output file for the histograms
  TFile *ofile = TFile::Open("ResolutionAnalysis_Exercise_Out.root","RECREATE");

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
  TH1D *trackMomentumRes = new TH1D("trackMomentumRes","Track Momentum Resolution; (P_{rec} - P_{MC})/P_{MC}", 400, -2, 2);
  TH2D* trackMomResP = new TH2D("trackMomResP", "Track Momentum Resolution vs P; (P_{rec} - P_{MC})/P_{MC}; P_{MC}(GeV/c)", 400, -2, 2, 150, 0, 150);
  TH2D* trackMomResEta = new TH2D("trackMomResEta", "Track Momentum Resolution vs #eta; (P_{rec} - P_{MC})/P_{MC}; #eta_{MC}", 400, -2, 2, 120, -6, 6);

  TH1D *trackMomentumRes_e = new TH1D("trackMomentumRes_e","e^{#pm} Track Momentum Resolution; (P_{rec} - P_{MC})/P_{MC}", 400, -2, 2);
  TH2D* trackMomResP_e = new TH2D("trackMomResP_e", "e^{#pm} Track Momentum Resolution vs P; (P_{rec} - P_{MC})/P_{MC}; P_{MC}(GeV/c)", 400, -2, 2, 150, 0, 25);
  TH2D* trackMomResEta_e = new TH2D("trackMomResEta_e", "e^{#pm} Track Momentum Resolution vs #eta; (P_{rec} - P_{MC})/P_{MC}; #eta_{MC}", 400, -2, 2, 120, -6, 6);

  TH1D *trackMomentumRes_mu = new TH1D("trackMomentumRes_mu","#mu^{#pm} Track Momentum Resolution; (P_{rec} - P_{MC})/P_{MC}", 400, -2, 2);
  TH2D* trackMomResP_mu = new TH2D("trackMomResP_mu", "#mu^{#pm} Track Momentum Resolution vs P; (P_{rec} - P_{MC})/P_{MC}; P_{MC}(GeV/c)", 400, -2, 2, 150, 0, 25);
  TH2D* trackMomResEta_mu = new TH2D("trackMomResEta_mu", "#mu^{#pm} Track Momentum Resolution vs #eta; (P_{rec} - P_{MC})/P_{MC}; #eta_{MC}", 400, -2, 2, 120, -6, 6);

  TH1D *trackMomentumRes_pi = new TH1D("trackMomentumRes_pi","#pi^{#pm} Track Momentum Resolution; (P_{rec} - P_{MC})/P_{MC}", 400, -2, 2);
  TH2D* trackMomResP_pi = new TH2D("trackMomResP_pi", "#pi^{#pm} Track Momentum Resolution vs P; (P_{rec} - P_{MC})/P_{MC}; P_{MC}(GeV/c)", 400, -2, 2, 150, 0, 150);
  TH2D* trackMomResEta_pi = new TH2D("trackMomResEta_pi", "#pi^{#pm} Track Momentum Resolution vs #eta; (P_{rec} - P_{MC})/P_{MC}; #eta_{MC}", 400, -2, 2, 120, -6, 6);

  TH1D *trackMomentumRes_K = new TH1D("trackMomentumRes_K","K^{#pm} Track Momentum Resolution; (P_{rec} - P_{MC})/P_{MC}", 400, -2, 2);
  TH2D* trackMomResP_K = new TH2D("trackMomResP_K", "K^{#pm} Track Momentum Resolution vs P; (P_{rec} - P_{MC})/P_{MC}; P_{MC}(GeV/c)", 400, -2, 2, 150, 0, 150);
  TH2D* trackMomResEta_K = new TH2D("trackMomResEta_K", "K^{#pm} Track Momentum Resolution vs #eta; (P_{rec} - P_{MC})/P_{MC}; #eta_{MC}", 400, -2, 2, 120, -6, 6);

  TH1D *trackMomentumRes_p = new TH1D("trackMomentumRes_p","p Track Momentum Resolution; (P_{rec} - P_{MC})/P_{MC}", 400, -2, 2);
  TH2D* trackMomResP_p = new TH2D("trackMomResP_p", "p Track Momentum Resolution vs P; (P_{rec} - P_{MC})/P_{MC}; P_{MC}(GeV/c)", 400, -2, 2, 150, 0, 150);
  TH2D* trackMomResEta_p = new TH2D("trackMomResEta_p", "p Track Momentum Resolution vs #eta; (P_{rec} - P_{MC})/P_{MC}; #eta_{MC}", 400, -2, 2, 120, -6, 6);
  
  TH1D *matchedPartTrackDeltaEta = new TH1D("matchedPartTrackDeltaEta","#Delta#eta Between Matching Thrown and Reconstructed Charged Particle; #Delta#eta", 100, -0.25, 0.25);
  TH1D *matchedPartTrackDeltaPhi = new TH1D("matchedPartTrackDeltaPhi","#Detla #phi Between Matching Thrown and Reconstructed Charged Particle; #Delta#phi", 200, -0.2, 0.2);
  TH1D *matchedPartTrackDeltaR = new TH1D("matchedPartTrackDeltaR","#Delta R Between Matching Thrown and Reconstructed Charged Particle; #Delta R", 300, 0, 0.3);
  TH1D *matchedPartTrackDeltaMom = new TH1D("matchedPartTrackDeltaMom","#Delta P Between Matching Thrown and Reconstructed Charged Particle; #Delta P", 200, -10, 10);

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
			float deltaMom = ((trueMom.Mag()) - (recMom.Mag()));

			double momRes = (recMom.Mag() - trueMom.Mag())/trueMom.Mag();
	
			trackMomentumRes->Fill(momRes); // Could also multiply by 100 and express as a percentage instead
			trackMomResP->Fill(momRes, trueMom.Mag());
			trackMomResEta->Fill(momRes, trueEta);

			if( pdg == 11){
			  trackMomentumRes_e->Fill(momRes);
			  trackMomResP_e->Fill(momRes, trueMom.Mag());
			  trackMomResEta_e->Fill(momRes, trueEta);
			}
			else if( pdg == 13){
			  trackMomentumRes_mu->Fill(momRes);
			  trackMomResP_mu->Fill(momRes, trueMom.Mag());
			  trackMomResEta_mu->Fill(momRes, trueEta);
			}
			else if( pdg == 211){
			  trackMomentumRes_pi->Fill(momRes);
			  trackMomResP_pi->Fill(momRes, trueMom.Mag());
			  trackMomResEta_pi->Fill(momRes, trueEta);
			}
			else if( pdg == 321){
			  trackMomentumRes_K->Fill(momRes);
			  trackMomResP_K->Fill(momRes, trueMom.Mag());
			  trackMomResEta_K->Fill(momRes, trueEta);
			}
			else if( pdg == 2212){
			  trackMomentumRes_p->Fill(momRes);
			  trackMomResP_p->Fill(momRes, trueMom.Mag());
			  trackMomResEta_p->Fill(momRes, trueEta);
			}
			  
			matchedPartTrackDeltaEta->Fill(deltaEta);
			matchedPartTrackDeltaPhi->Fill(deltaPhi);
			matchedPartTrackDeltaR->Fill(deltaR);
			matchedPartTrackDeltaMom->Fill(deltaMom);
			
		      }
		  }// End loop over associations
	      } // End PDG check
	  } // End stable particles condition
      } // End loop over thrown particles
  } // End loop over events

  ofile->Write(); // Write histograms to file
  ofile->Close(); // Close output file
}
```

Insert your input file path and execute as the example code above.

### Compiled ROOT Scripts 

As brought up in the tutorial, you may wish to compile your ROOT based scripts for faster processing. Included below are some scripts and a short example of a compiled ROOT macro provided by Kolja Kauder.

Each file is uploaded invidually, but your directory should be structured as follows -

- helloroot
    - README.md
    - CMakeLists.txt
    - include
      - helloroot
        - helloroot.hh
    - src
        - helloexec.cxx
        - helloroot.cxx
     
Note that any entry in the above without a file extension is a directory.

The contents of README.md are -

``c++
To build using cmake, create a build directory, navigate to it and run cmake. e.g.:

```
mkdir build
cd build
cmake .. 
make 
```
You can specify a number of parallel build threads with the -j flag, e.g.
```
make -j4
```

You can specify an install directory to cmake with
-DCMAKE_INSTALL_PREFIX=<path>
then, after building, 
```
make install
```
to install the headers and libraries under that location.
There is no "make uninstall" but (on Unix-like systems)
you can do
xargs rm < install_manifest.txt
from the cmake build directory.
```

The contents of CMakeLists.txt are -

```c++
# CMakeLists.txt for helloroot.
# More complicated than needed but demonstrates making and linking your own libraries
# cf. https://cliutils.gitlab.io/modern-cmake/chapters/packages/ROOT.html
# https://root.cern/manual/integrate_root_into_my_cmake_project/

cmake_minimum_required(VERSION 3.10)
project(helloroot VERSION 1.0 LANGUAGES CXX ) # not needed

 # Find ROOT. Use at least 6.20 for smoother cmake support
 find_package(ROOT 6.20 REQUIRED )

 message ( " ROOT Libraries = " ${ROOT_LIBRARIES} )

 ##############################################################################################################

# Main target is the libhelloroot library
add_library(
  # You can use wildcards but it's cleaner to list the files explicitly
   helloroot
   SHARED
   src/helloroot.cxx
   )
## The particular syntax here is a bit annoying because you have to list all the sub-modules you need
## but it picks up automatically all the compile options needed for root, e.g. the c++ std version
## Find all available ROOT modules with `root-config --libs`
target_link_libraries(helloroot PUBLIC ROOT::Core ROOT::RIO ROOT::Rint ROOT::Tree ROOT::EG ROOT::Physics )

## The above _should_ be true, and it is on most systems. If it's not, uncoment one of the following lines
# target_compile_features(helloroot PUBLIC cxx_std_17)
# target_compile_features(helloroot PUBLIC cxx_std_20)

# include directories - this is also overkill but useful if you want to create dictionaries
# Contact kkauder@gmail.com for that - it's too much for this example
target_include_directories(helloroot 
PUBLIC 
$<INSTALL_INTERFACE:include>
$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
 )

# Can add addtional options here
target_compile_options(helloroot PRIVATE -Wall -Wextra -pedantic -g)  

##############################################################################################################

## Build executables
add_executable(helloexec src/helloexec.cxx)
# target_compile_options(helloexec PRIVATE -Wall -Wextra -pedantic -g)
target_link_libraries(helloexec helloroot )
target_include_directories(helloexec
  PRIVATE
  ${ROOT_INCLUDE_DIRS}
  )

install(TARGETS helloexec DESTINATION bin)


##############################################################################################################

## Install library
# Could also use include(GNUInstallDirs)
# and then destinations of the form ${CMAKE_INSTALL_INCLUDEDIR}
install(TARGETS helloroot
  EXPORT helloroot-export
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  )

## Install headers
install (DIRECTORY ${CMAKE_SOURCE_DIR}/include/helloroot
  DESTINATION  include/helloroot
  )

## Generate configuration file - this allows you to use cmake in another project 
## to find and link the installed helloroot library
install(EXPORT helloroot-export
  FILE
  hellorootConfig.cmake
  NAMESPACE
    helloroot::
  DESTINATION
  cmake
  )

## Final message
message( " Done!")
```

The contents of helloroot.hh are -

```c++
#ifndef HELLO_ROOT_H
#define HELLO_ROOT_H

void HelloRoot();

#endif // HELLO_ROOT_H
```

The contents of helloexec.cxx are -

```c++
#include<helloroot/helloroot.hh>

#include<iostream>
#include<string>

int main()
{
    std::cout << "Hello from main " << std::endl;
    HelloRoot(); 
    
    return 0;
}
```

And finally, the contents of helloroot.cxx are -

```c++
#include<helloroot/helloroot.hh>

#include<iostream>
#include<string>

#include<TH1D.h>
#include<TPad.h>

void HelloRoot()
{
  std::cout << "Hello from HelloRoot" << std::endl;

  // do something with root
  TH1D h("h", "h", 100, -5, 5);
  h.FillRandom("gaus", 1000);
  h.Draw();
  gPad->SaveAs("hello.png");

  return;
}
```
Please consult the README and script comments for further instructions.

## Python Uproot Scripts

### EfficiencyAnalysis.py

Create a file called `EfficiencyAnalysis.py` and copy in the code below to get started on the resolution analysis exercise. Note that you will need to correctly specifiy your input file path in the variable `infile`.

```python
#! /usr/bin/python

#Import relevant packages
import ROOT, math, array
from ROOT import TH1F, TH2F, TMath, TTree, TVector3, TVector2
import uproot as up

#Define and open files
infile="PATH_TO_INPUT_FILE"
ofile=ROOT.TFile.Open("EfficiencyAnalysis_OutPy.root", "RECREATE")

# Open input file and define branches we want to look at with uproot
events_tree = up.open(infile)["events"]

# Get particle information
partGenStat = events_tree["MCParticles.generatorStatus"].array()
partMomX = events_tree["MCParticles.momentum.x"].array()
partMomY = events_tree["MCParticles.momentum.y"].array()
partMomZ = events_tree["MCParticles.momentum.z"].array()
partPdg = events_tree["MCParticles.PDG"].array()

# Get reconstructed track information
trackMomX = events_tree["ReconstructedChargedParticles.momentum.x"].array()
trackMomY = events_tree["ReconstructedChargedParticles.momentum.y"].array()
trackMomZ = events_tree["ReconstructedChargedParticles.momentum.z"].array()

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

A "solution" version of the script for the exercise is included below -

```python
#! /usr/bin/python

#Import relevant packages
import ROOT, math, array
from ROOT import TCanvas, TColor, TGaxis, TH1F, TH2F, TPad, TStyle, gStyle, gPad, TGaxis, TLine, TMath, TPaveText, TTree, TVector3, TVector2
import uproot as up

#Define and open files
infile="PATH_TO_FILE"
ofile=ROOT.TFile.Open("EfficiencyAnalysis_Exercise_OutPy.root", "RECREATE")

# Open input file and define branches we want to look at with uproot
events_tree = up.open(infile)["events"]

# Get particle information
partGenStat = events_tree["MCParticles.generatorStatus"].array()
partMomX = events_tree["MCParticles.momentum.x"].array()
partMomY = events_tree["MCParticles.momentum.y"].array()
partMomZ = events_tree["MCParticles.momentum.z"].array()
partPdg = events_tree["MCParticles.PDG"].array()

# Get reconstructed track information
trackMomX = events_tree["ReconstructedChargedParticles.momentum.x"].array()
trackMomY = events_tree["ReconstructedChargedParticles.momentum.y"].array()
trackMomZ = events_tree["ReconstructedChargedParticles.momentum.z"].array()

# Get assocations between MCParticles and ReconstructedChargedParticles
recoAssoc = events_tree["ReconstructedChargedParticleAssociations.recID"].array()
simuAssoc = events_tree["ReconstructedChargedParticleAssociations.simID"].array()

# Define histograms below
partEta = ROOT.TH1D("partEta","#eta of Thrown Charged Particles; #eta", 120, -6, 6)
matchedPartEta = ROOT.TH1D("matchedPartEta","#eta of Thrown Charged Particles That Have Matching Track; #eta", 120, -6, 6)
partMom = ROOT.TH1D("partMom", "Momentum of Thrown Charged Particles (truth); P(GeV/c)", 150, 0, 150)
matchedPartMom = ROOT.TH1D("matchedPartMom", "Momentum of Thrown Charged Particles (truth), with matching track; P(GeV/c)", 150, 0, 150)
partPhi = ROOT.TH1D("partPhi", "#phi of Thrown Charged Particles (truth); #phi(rad)", 320, -3.2, 3.2)
matchedPartPhi = ROOT.TH1D("matchedPartPhi", "#phi of Thrown Charged Particles (truth), with matching track; #phi(rad)", 320, -3.2, 3.2)

partPEta = ROOT.TH2D("partPEta", "P vs #eta of Thrown Charged Particles; P(GeV/c); #eta", 150, 0, 150, 120, -6, 6)
matchedPartPEta = ROOT.TH2D("matchedPartPEta", "P vs #eta of Thrown Charged Particles, with matching track; P(GeV/C); #eta", 150, 0, 150, 120, -6, 6)
partPhiEta = ROOT.TH2D("partPhiEta", "#phi vs #eta of Thrown Charged Particles; #phi(rad); #eta", 160, -3.2, 3.2, 120, -6, 6)
matchedPartPhiEta = ROOT.TH2D("matchedPartPhiEta", "#phi vs #eta of Thrown Charged Particles; #phi(rad); #eta", 160, -3.2, 3.2, 120, -6, 6)

matchedPartTrackDeltaEta = ROOT.TH1D("matchedPartTrackDeltaEta","#Delta#eta Between Matching Thrown and Reconstructe Charged Particle; #Delta#eta", 100, -0.25, 0.25)
matchedPartTrackDeltaPhi = ROOT.TH1D("matchedPartTrackDeltaPhi","#Detla #phi Between Matching Thrown and Reconstructed Charged Particle; #Delta#phi", 200, -0.2, 0.2)
matchedPartTrackDeltaR = ROOT.TH1D("matchedPartTrackDeltaR","#Delta R Between Matching Thrown and Reconstructed Charged Particle; #Delta R", 300, 0, 0.3)
matchedPartTrackDeltaMom = ROOT.TH1D("matchedPartTrackDeltaMom","#Delta P Between Matching Thrown and Reconstructed Charged Particle; #Delta P", 200, -10, 10)

# Define some histograms for our efficiencies
TrackEff_Eta = ROOT.TH1D("TrackEff_Eta", "Tracking efficiency as fn of #eta; #eta; Eff(%)", 120, -6, 6)
TrackEffMom = ROOT.TH1D("TrackEff_Mom", "Tracking efficiency as fn of P; P(GeV/c); Eff(%)", 150, 0, 150)
TrackEffPhi = ROOT.TH1D("TrackEff_Phi", "Tracking efficiency as fn of #phi; #phi(rad); Eff(%)", 320, -3.2, 3.2)

# 2D Efficiencies
TrackEff_PEta = ROOT.TH2D("TrackEff_PEta", "Tracking efficiency as fn of P and #eta; P(GeV/c); #eta", 150, 0, 150, 120, -6, 6)
TrackEff_PhiEta = ROOT.TH2D("TrackEff_PhiEta", "Tracking efficiency as fn of #phi and #eta; #phi(rad); #eta", 160, -3.2, 3.2, 120, -6, 6)

# All charged particle histos
ChargedEta = ROOT.TH1D("ChargedEta", "#eta of all charged particles; #eta", 120, -6, 6)
ChargedPhi = ROOT.TH1D("ChargedPhi", "#phi of all charged particles; #phi (rad)", 120, -3.2, 3.2)
ChargedP = ROOT.TH1D("ChargedP", "P of all charged particles; P(GeV/c)", 150, 0, 150)

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
                partPhi.Fill(truePhi)
                partMom.Fill(trueMom.Mag())
                partPEta.Fill(trueMom.Mag(), trueEta)
                partPhiEta.Fill(truePhi, trueEta)
                for k in range(0,len(simuAssoc[i])): # Loop over associations to find matching ReconstructedChargedParticle
                    if (simuAssoc[i][k] == j):
                        recMom = ROOT.TVector3(trackMomX[i][recoAssoc[i][k]], trackMomY[i][recoAssoc[i][k]], trackMomZ[i][recoAssoc[i][k]])
                        deltaEta = trueEta - recMom.PseudoRapidity()
                        deltaPhi = TVector2. Phi_mpi_pi(truePhi - recMom.Phi())
                        deltaR = math.sqrt((deltaEta*deltaEta) + (deltaPhi*deltaPhi))
                        deltaMom = ((trueMom.Mag()) - (recMom.Mag()))

                        matchedPartTrackDeltaEta.Fill(deltaEta)
                        matchedPartTrackDeltaPhi.Fill(deltaPhi)
                        matchedPartTrackDeltaR.Fill(deltaR)
                        matchedPartTrackDeltaMom.Fill(deltaMom)
                        
                        matchedPartEta.Fill(trueEta)
                        matchedPartPhi.Fill(truePhi)
                        matchedPartMom.Fill(trueMom.Mag())
                        matchedPartPEta.Fill(trueMom.Mag(), trueEta)
                        matchedPartPhiEta.Fill(truePhi, trueEta)
    for x in range (0, len(trackMomX[i])): # Loop over all charged particles, thrown or not
        CPartMom = ROOT.TVector3(trackMomX[i][x], trackMomY[i][x], trackMomZ[i][x])
        CPartEta = CPartMom.PseudoRapidity()
        CPartPhi = CPartMom.Phi()

        ChargedEta.Fill(CPartEta)
        ChargedPhi.Fill(CPartPhi)
        ChargedP.Fill(CPartMom.Mag())
        
# Write output histograms to file below
partEta.Write()
matchedPartEta.Write()
partMom.Write()
matchedPartMom.Write()
partPhi.Write()
matchedPartPhi.Write()
partPEta.Write()
matchedPartPEta.Write()
partPhiEta.Write()
matchedPartPhiEta.Write()
matchedPartTrackDeltaEta.Write()
matchedPartTrackDeltaPhi.Write()
matchedPartTrackDeltaR.Write()
matchedPartTrackDeltaMom.Write()
ChargedEta.Write()
ChargedPhi.Write()
ChargedP.Write()

TrackEff_Eta.Divide(matchedPartEta, partEta, 1, 1, "b")
TrackEffMom.Divide(matchedPartMom, partMom, 1, 1, "b")
TrackEffPhi.Divide(matchedPartPhi, partPhi, 1, 1, "b")
TrackEff_PEta.Divide(matchedPartPEta, partPEta, 1, 1, "b")
TrackEff_PhiEta.Divide(matchedPartPhiEta, partPhiEta, 1, 1, "b")

TrackEff_Eta.Write()
TrackEffMom.Write()
TrackEffPhi.Write()
TrackEff_PEta.Write()
TrackEff_PhiEta.Write()

# Close files
ofile.Close()
```
Insert your input file path and execute as the example code above.

### ResolutionAnalysis.py

Create a file called `ResolutionAnalysis.py` and copy in the code below to get started on the resolution analysis exercise. Note that you will need to correctly specifiy your input file path in the variable `infile`.

```python
#! /usr/bin/python

#Import relevant packages
import ROOT, math, array
from ROOT import TH1F, TH2F, TMath, TTree, TVector3, TVector2
import uproot as up

#Define and open files
infile="PATH_TO_INPUT_FILE"
ofile=ROOT.TFile.Open("ResolutionAnalysis_OutPy.root", "RECREATE")

# Open input file and define branches we want to look at with uproot
events_tree = up.open(infile)["events"]

# Get particle information
partGenStat = events_tree["MCParticles.generatorStatus"].array()
partMomX = events_tree["MCParticles.momentum.x"].array()
partMomY = events_tree["MCParticles.momentum.y"].array()
partMomZ = events_tree["MCParticles.momentum.z"].array()
partPdg = events_tree["MCParticles.PDG"].array()

# Get reconstructed track information
trackMomX = events_tree["ReconstructedChargedParticles.momentum.x"].array()
trackMomY = events_tree["ReconstructedChargedParticles.momentum.y"].array()
trackMomZ = events_tree["ReconstructedChargedParticles.momentum.z"].array()

# Get assocations between MCParticles and ReconstructedChargedParticles
recoAssoc = events_tree["ReconstructedChargedParticleAssociations.recID"].array()
simuAssoc = events_tree["ReconstructedChargedParticleAssociations.simID"].array()

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

A "solution" version of the script for the exercise is included below -

```python
#! /usr/bin/python

#Import relevant packages
import ROOT, math, array
from ROOT import TCanvas, TColor, TGaxis, TH1F, TH2F, TPad, TStyle, gStyle, gPad, TGaxis, TLine, TMath, TPaveText, TTree, TVector3, TVector2
import uproot as up

#Define and open files
infile="PATH_TO_FILE"
ofile=ROOT.TFile.Open("ResolutionAnalysis_Exercise_OutPy.root", "RECREATE")

# Open input file and define branches we want to look at with uproot
events_tree = up.open(infile)["events"]

# Get particle information
partGenStat = events_tree["MCParticles.generatorStatus"].array()
partMomX = events_tree["MCParticles.momentum.x"].array()
partMomY = events_tree["MCParticles.momentum.y"].array()
partMomZ = events_tree["MCParticles.momentum.z"].array()
partPdg = events_tree["MCParticles.PDG"].array()

# Get reconstructed track information
trackMomX = events_tree["ReconstructedChargedParticles.momentum.x"].array()
trackMomY = events_tree["ReconstructedChargedParticles.momentum.y"].array()
trackMomZ = events_tree["ReconstructedChargedParticles.momentum.z"].array()

# Get assocations between MCParticles and ReconstructedChargedParticles
recoAssoc = events_tree["ReconstructedChargedParticleAssociations.recID"].array()
simuAssoc = events_tree.["ReconstructedChargedParticleAssociations.simID"].array()

# Define histograms below
trackMomentumRes = ROOT.TH1D("trackMomentumRes","Track Momentum Resolution", 400, -2, 2)
trackMomResP = ROOT.TH2D("trackMomResP", "Track Momentum Resolution vs P; (P_{rec} - P_{MC})/P_{MC}; P_{MC}(GeV/c)", 400, -2, 2, 150, 0, 150);
trackMomResEta = ROOT.TH2D("trackMomResEta", "Track Momentum Resolution vs #eta; (P_{rec} - P_{MC})/P_{MC}; #eta_{MC}", 400, -2, 2, 120, -6, 6);

trackMomentumRes_e = ROOT.TH1D("trackMomentumRes_e","e^{#pm} Track Momentum Resolution; (P_{rec} - P_{MC})/P_{MC}", 400, -2, 2);
trackMomResP_e = ROOT.TH2D("trackMomResP_e", "e^{#pm} Track Momentum Resolution vs P; (P_{rec} - P_{MC})/P_{MC}; P_{MC}(GeV/c)", 400, -2, 2, 150, 0, 25);
trackMomResEta_e = ROOT.TH2D("trackMomResEta_e", "e^{#pm} Track Momentum Resolution vs #eta; (P_{rec} - P_{MC})/P_{MC}; #eta_{MC}", 400, -2, 2, 120, -6, 6);

trackMomentumRes_mu = ROOT.TH1D("trackMomentumRes_mu","#mu^{#pm} Track Momentum Resolution; (P_{rec} - P_{MC})/P_{MC}", 400, -2, 2);
trackMomResP_mu = ROOT.TH2D("trackMomResP_mu", "#mu^{#pm} Track Momentum Resolution vs P; (P_{rec} - P_{MC})/P_{MC}; P_{MC}(GeV/c)", 400, -2, 2, 150, 0, 25);
trackMomResEta_mu = ROOT.TH2D("trackMomResEta_mu", "#mu^{#pm} Track Momentum Resolution vs #eta; (P_{rec} - P_{MC})/P_{MC}; #eta_{MC}", 400, -2, 2, 120, -6, 6);

trackMomentumRes_pi = ROOT.TH1D("trackMomentumRes_pi","#pi^{#pm} Track Momentum Resolution; (P_{rec} - P_{MC})/P_{MC}", 400, -2, 2);
trackMomResP_pi = ROOT.TH2D("trackMomResP_pi", "#pi^{#pm} Track Momentum Resolution vs P; (P_{rec} - P_{MC})/P_{MC}; P_{MC}(GeV/c)", 400, -2, 2, 150, 0, 150);
trackMomResEta_pi = ROOT.TH2D("trackMomResEta_pi", "#pi^{#pm} Track Momentum Resolution vs #eta; (P_{rec} - P_{MC})/P_{MC}; #eta_{MC}", 400, -2, 2, 120, -6, 6);

trackMomentumRes_K = ROOT.TH1D("trackMomentumRes_K","K^{#pm} Track Momentum Resolution; (P_{rec} - P_{MC})/P_{MC}", 400, -2, 2);
trackMomResP_K = ROOT.TH2D("trackMomResP_K", "K^{#pm} Track Momentum Resolution vs P; (P_{rec} - P_{MC})/P_{MC}; P_{MC}(GeV/c)", 400, -2, 2, 150, 0, 150);
trackMomResEta_K = ROOT.TH2D("trackMomResEta_K", "K^{#pm} Track Momentum Resolution vs #eta; (P_{rec} - P_{MC})/P_{MC}; #eta_{MC}", 400, -2, 2, 120, -6, 6);

trackMomentumRes_p = ROOT.TH1D("trackMomentumRes_p","p Track Momentum Resolution; (P_{rec} - P_{MC})/P_{MC}", 400, -2, 2);
trackMomResP_p = ROOT.TH2D("trackMomResP_p", "p Track Momentum Resolution vs P; (P_{rec} - P_{MC})/P_{MC}; P_{MC}(GeV/c)", 400, -2, 2, 150, 0, 150);
trackMomResEta_p = ROOT.TH2D("trackMomResEta_p", "p Track Momentum Resolution vs #eta; (P_{rec} - P_{MC})/P_{MC}; #eta_{MC}", 400, -2, 2, 120, -6, 6);

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

                        trackMomentumRes.Fill(momRes)
                        trackMomResP.Fill(momRes, trueMom.Mag())
                        trackMomResEta.Fill(momRes, trueEta)
                        
                        if( pdg == 11):
                            trackMomentumRes_e.Fill(momRes)
                            trackMomResP_e.Fill(momRes, trueMom.Mag())
                            trackMomResEta_e.Fill(momRes, trueEta)
                        elif( pdg == 13):
                            trackMomentumRes_mu.Fill(momRes)
                            trackMomResP_mu.Fill(momRes, trueMom.Mag())
                            trackMomResEta_mu.Fill(momRes, trueEta)
                        elif( pdg == 211):
                            trackMomentumRes_pi.Fill(momRes)
                            trackMomResP_pi.Fill(momRes, trueMom.Mag())
                            trackMomResEta_pi.Fill(momRes, trueEta)
                        elif( pdg == 321):
                            trackMomentumRes_K.Fill(momRes)
                            trackMomResP_K.Fill(momRes, trueMom.Mag())
                            trackMomResEta_K.Fill(momRes, trueEta)
                        elif( pdg == 2212):
                            trackMomentumRes_p.Fill(momRes)
                            trackMomResP_p.Fill(momRes, trueMom.Mag())
                            trackMomResEta_p.Fill(momRes, trueEta)
                            
                        matchedPartTrackDeltaEta.Fill(deltaEta)
                        matchedPartTrackDeltaPhi.Fill(deltaPhi)
                        matchedPartTrackDeltaR.Fill(deltaR)
                        matchedPartTrackDeltaMom.Fill(deltaMom)
                        
# Write output histograms to file below
trackMomentumRes.Write()
trackMomResP.Write()
trackMomResEta.Write()
trackMomentumRes_e.Write()
trackMomResP_e.Write()
trackMomResEta_e.Write()
trackMomentumRes_mu.Write()
trackMomResP_mu.Write()
trackMomResEta_mu.Write()
trackMomentumRes_pi.Write()
trackMomResP_pi.Write()
trackMomResEta_pi.Write()
trackMomentumRes_K.Write()
trackMomResP_K.Write()
trackMomResEta_K.Write()
trackMomentumRes_p.Write()
trackMomResP_p.Write()
trackMomResEta_p.Write()

matchedPartTrackDeltaEta.Write()
matchedPartTrackDeltaPhi.Write()
matchedPartTrackDeltaR.Write()
matchedPartTrackDeltaMom.Write()

# Close files
ofile.Close()
```
Insert your input file path and execute as the example code above.


## RDataFrames Example

Note that only the initial stage of the efficiency example is presented here in RDF format. This example was kindly created by [Simon](https://github.com/simonge/EIC_Analysis/blob/main/Analysis-Tutorial/EfficiencyAnalysisRDF.C).

### EfficiencyAnalysisRDF.C

Create a file called `EfficiencyAnalysisRDF.C` and paste the code below in. Remember to change the file path. 

Execute this script via - `root -l -q EfficiencyAnalysisRDF.C++`. Do this within eic-shell or somewhere else with the correct EDM4hep/EDM4eic libraries installed.

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
> - I'm not as familiar with RDF's and I have not created a "complete" example using this approach or a resolution analysis version. 
> - Working on a full version of this still as of 30/04/24 - Check back in the future for updates.
{: .callout}

{% include links.md %}
