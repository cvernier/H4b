#root -l -b -q '../HbbHbbPileupWeight.c++("glugluToX600ToHHTobbbb_8TeV_FASTSIM_PU_Signal_2", "8TeVData2012BCD_V5_Data_2", "PUWeight.root")'
#echo "Pileup Weights for Signal Monte Carlo created."
#root -l -b -q 'HbbHbb_KinematicSelection.cc+("/gpfs/ddn/cms/user/cvernier/H4b_step2/CMSSW_5_3_3_patch2/src/VHbbAnalysis/VHbbDataFormats/bin/", "DiJetPt_BJetPlusX_Run2012B-13Jul2012-v1_EDMNtuple_V42_ProcV1", "")'
#root -l -b -q 'HbbHbb_KinematicSelection.cc+("/gpfs/ddn/cms/user/cvernier/H4b_step2/CMSSW_5_3_3_patch2/src/VHbbAnalysis/VHbbDataFormats/bin/", "DiJetPt_BJetPlusX_Run2012C-24Aug2012-v2_EDMNtuple_V42_ProcV1", "")'
#root -l -b -q 'HbbHbb_KinematicSelection.cc+("/gpfs/ddn/cms/user/cvernier/H4b_step2/CMSSW_5_3_3_patch2/src/VHbbAnalysis/VHbbDataFormats/bin/", "DiJetPt_BJetPlusX_Run2012C-PromptReco-v2_EDMNtuple_V42_ProcV1", "")'
#root -l -b -q 'HbbHbb_KinematicSelection.cc+("/gpfs/ddn/cms/user/cvernier/H4b_step2/CMSSW_5_3_3_patch2/src/VHbbAnalysis/VHbbDataFormats/bin/", "DiJetPt_BJetPlusX_Run2012D-PromptReco-v1_EDMNtuple_V42_ProcV1", "")'

root -l -b -q '../../HbbHbb_bTagging.cc++("Data", "DiJetPt_BJetPlusX_Run2012B-13Jul2012-v1_EDMNtuple_V42_ProcV1", "MMMM","nominal")' 
echo "2012B b-tagging done"
root -l -b -q '../../HbbHbb_bTagging.cc++("Data", "DiJetPt_BJetPlusX_Run2012C-24Aug2012-v2_EDMNtuple_V42_ProcV1", "MMMM","nominal")' 
echo "2012C b-tagging done"
root -l -b -q '../../HbbHbb_bTagging.cc++("Data", "DiJetPt_BJetPlusX_Run2012C-PromptReco-v2_EDMNtuple_V42_ProcV1", "MMMM","nominal")' 
echo "2012C part 2 b-tagging done"
root -l -b -q '../../HbbHbb_bTagging.cc++("Data", "DiJetPt_BJetPlusX_Run2012D-PromptReco-v1_EDMNtuple_V42_ProcV1", "MMMM","nominal")' 
echo "2012D b-tagging done"
