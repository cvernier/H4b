{
  gSystem->AddIncludePath("-I/Users/souvik/CMSSW_4_2_3_FWLITE/src");
	gSystem->Load("libPhysicsToolsKinFitter.so");
	gROOT->ProcessLine(".L ../../HbbHbb_Partitioning.cc++O");
	HbbHbb_Partitioning("BJetPlusX_Run2012BCD_Skim", 142.5, 1);
	gROOT->ProcessLine(".q");	
}
