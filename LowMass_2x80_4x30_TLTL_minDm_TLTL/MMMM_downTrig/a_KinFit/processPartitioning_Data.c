{
  gSystem->AddIncludePath("-I/Users/souvik/CMSSW_4_2_3_FWLITE/src");
	gSystem->Load("libPhysicsToolsKinFitter.so");
	
	.x ../../../HbbHbb_Partitioning.cc++("BJetPlusX_Run2012B-13Jul2012-v1_EDMNtuple_V42_ProcV1_Skim", 125., 1);
	.x ../../../HbbHbb_Partitioning.cc++("BJetPlusX_Run2012C-24Aug2012-v2_EDMNtuple_V42_ProcV1_Skim", 125., 1);
	.x ../../../HbbHbb_Partitioning.cc++("BJetPlusX_Run2012C-PromptReco-v2_EDMNtuple_V42_ProcV1_Skim", 125., 1);
	.x ../../../HbbHbb_Partitioning.cc++("BJetPlusX_Run2012D-PromptReco-v1_EDMNtuple_V42_ProcV1_Skim", 125., 1);
	
}
