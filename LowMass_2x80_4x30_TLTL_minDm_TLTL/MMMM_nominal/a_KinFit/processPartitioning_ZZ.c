{
  gSystem->AddIncludePath("-I/Users/souvik/CMSSW_4_2_3_FWLITE/src");
	gSystem->Load("libPhysicsToolsKinFitter.so");
	gROOT->ProcessLine(".L ../../../HbbHbb_Partitioning.cc++O");
	HbbHbb_Partitioning("DiJetPt_ZZ_TuneZ2star_8TeV_pythia6_tauola", 91., 1);
	gROOT->ProcessLine(".q");	
}
