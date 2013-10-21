{
  gSystem->AddIncludePath("-I/Users/souvik/CMSSW_4_2_3_FWLITE/src");
	gSystem->Load("libPhysicsToolsKinFitter.so");
	.x ../../../HbbHbb_Partitioning.cc++("glugluToX300ToHHTobbbb_8TeV_width1MeV_FASTSIM_PU", 125., 1);
	.x ../../../HbbHbb_Partitioning.cc++("glugluToX400ToHHTobbbb_8TeV_width1MeV_FASTSIM_PU", 125., 1);
	.x ../../../HbbHbb_Partitioning.cc++("glugluToX500ToHHTobbbb_8TeV_width1MeV_FASTSIM_PU", 125., 1);
	.x ../../../HbbHbb_Partitioning.cc++("glugluToX600ToHHTobbbb_8TeV_width1MeV_FASTSIM_PU", 125., 1);
	.x ../../../HbbHbb_Partitioning.cc++("glugluToX700ToHHTobbbb_8TeV_width1MeV_FASTSIM_PU", 125., 1);
	.x ../../../HbbHbb_Partitioning.cc++("glugluToX800ToHHTobbbb_8TeV_width1MeV_FASTSIM_PU", 125., 1);
	
}
