{
  gSystem->AddIncludePath("-I/Users/souvik/CMSSW_4_2_3_FWLITE/src");
	gSystem->Load("libPhysicsToolsKinFitter.so");
	gROOT->ProcessLine(".L ../../../HbbHbb_Partitioning.cc++O");
//	HbbHbb_Partitioning("DiJetPt_RadionToHH_4b_M-270_TuneZ2star_8TeV-Madgraph_pythia6", 125., 1);
	HbbHbb_Partitioning("DiJetPt_glugluToX300ToHHTobbbb_8TeV_width1MeV_GEN_FASTSIM_HLT", 125., 1);
//	HbbHbb_Partitioning("DiJetPt_glugluToX400ToHHTobbbb_8TeV_width1MeV_GEN_FASTSIM_HLT", 125., 1);
	gROOT->ProcessLine(".q");	
}
