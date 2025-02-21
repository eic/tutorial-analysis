---
title: "Tree Pruning Script"
---

Included below is a short script which can be utilised to prune an input tree. This can be utilised to trim down our large EICrecon output file to a smaller, simpler subset of trees. You can modify this as needed by following the comments in the script.

## TreePrune.C

Copy paste this to a new file called "TreePrune.C". Execute it via -

```console
root -l TreePrune.C
```

You will be prompted for an input file name. Alternatively run the script with the input already specified -

```console
root -l 'TreePrune.C("InputFilePath")'
```

```c+
+// Stephen JD Kay, University of York, 21/02/25
// A short script to read in a file and prune it to only retain a smaller subset of branches.
// This could be utilised to trim down a full EICrecon file for a new user to look at. This avoids the potentially overwhelming number of branches normally stored
#include <string>

void TreePrune(TString infile=""){

  if(infile == ""){
    cout << "Enter a filename to analyse: ";
    cin >> infile;
  }

  if(gSystem->AccessPathName(infile) == kTRUE){
    cerr << "!!!!! ERROR !!!!!" << endl << infile << " not found" << endl << "!!!!! ERROR !!!!!" << endl;
    exit(0);
  }
  else{
    cout << "Opening " << infile << endl;
  }

  TString ofile_name, ofile_tmp, BeamE;
  TObjArray *tmp_Name_arr;

  if(infile.Contains(".root") == false){
    cout << "!!!!!!!!!!! - ERROR - !!!!!!!!!!!!" << endl;
    cout << "Input files should be a root file!" << endl;
    cout << "!!!!!!!!!!! - ERROR - !!!!!!!!!!!!" << endl;
    exit(1);
  }
  
  if(infile.Contains("/")){
    tmp_Name_arr = infile.Tokenize("/");
    ofile_tmp = (((TObjString *)(tmp_Name_arr->At(tmp_Name_arr->GetLast())))->String()).ReplaceAll(".root","");
  }
  else{
    ofile_tmp = infile;
    ofile_tmp.ReplaceAll(".root","");
  }
  // Set output file name
  ofile_name = Form("%s_Pruned.root", ofile_tmp.Data());  

  // Open our full, unpruned file
  TFile *full_file = TFile::Open(infile);
  TTree* full_tree;
  // Get the events tree
  full_file->GetObject("events", full_tree);
  // Deactivate all branches
  full_tree->SetBranchStatus("*", 0);
  // Activate only the branches we want to keep, add more as needed via the line below
  // full_tree->SetBranchStatus("",1)
  full_tree->SetBranchStatus("MCParticles*",1);
  full_tree->SetBranchStatus("ReconstructedParticles**",1);
  full_tree->SetBranchStatus("ReconstructedChargedParticles*",1);
  full_tree->SetBranchStatus("ReconstructedParticleAssociations*",1);
  full_tree->SetBranchStatus("ReconstructedChargedParticleAssociations*",1);
  full_tree->SetBranchStatus("EcalEndcapNClusters*",1);
  full_tree->SetBranchStatus("EcalEndcapPClusters*",1);
  
  // Set and open output file for the histograms
  TFile *ofile = TFile::Open(ofile_name,"RECREATE");
  // Clone the tree
  TTree* pruned_tree = full_tree->CloneTree(0);
  // Copy the branches
  pruned_tree->CopyEntries(full_tree);

  // Write the new tree to file
  ofile->Write();
  ofile->Close(); // Close output file
  full_file->Close(); // Close input file

  delete ofile;
  delete full_file;
  
}
```

