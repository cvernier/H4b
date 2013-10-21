source processKinSel_Signal.sh 
echo "nominal"
cd /gpfs/ddn/cms/user/cvernier/H4b/CMSSW_5_3_3_patch2/src/UserCode/SouvikDas/HbbHbb/Analysis/PDF/Archive/LowMass_2x80_4x30_TLTL_minDm_TLTL/MMMM_nominal
source processBTagging_Signal.sh
cd a_KinFit
root -l -b processPartitioning_Signal.c
echo "trigUp"
cd /gpfs/ddn/cms/user/cvernier/H4b/CMSSW_5_3_3_patch2/src/UserCode/SouvikDas/HbbHbb/Analysis/PDF/Archive/LowMass_2x80_4x30_TLTL_minDm_TLTL/MMMM_upTrig
source processBTagging_Signal.sh
cd a_KinFit
root -l -b processPartitioning_Signal.c
echo "trigDown"
cd /gpfs/ddn/cms/user/cvernier/H4b/CMSSW_5_3_3_patch2/src/UserCode/SouvikDas/HbbHbb/Analysis/PDF/Archive/LowMass_2x80_4x30_TLTL_minDm_TLTL/MMMM_downTrig
source processBTagging_Signal.sh
cd a_KinFit
root -l -b processPartitioning_Signal.c
echo "JECm1"
cd /gpfs/ddn/cms/user/cvernier/H4b/CMSSW_5_3_3_patch2/src/UserCode/SouvikDas/HbbHbb/Analysis/PDF/Archive/LowMass_2x80_4x30_TLTL_minDm_TLTL/KinSel_JECm1
source processKinSel_Signal.sh
cd MMMM_nominal
source processBTagging_Signal.sh
cd a_KinFit
root -l -b processPartitioning_Signal.c
echo "JECp1"
cd /gpfs/ddn/cms/user/cvernier/H4b/CMSSW_5_3_3_patch2/src/UserCode/SouvikDas/HbbHbb/Analysis/PDF/Archive/LowMass_2x80_4x30_TLTL_minDm_TLTL/KinSel_JECp1
source processKinSel_Signal.sh
cd MMMM_nominal
source processBTagging_Signal.sh
cd a_KinFit
root -l -b processPartitioning_Signal.c
echo "JERm1"
cd /gpfs/ddn/cms/user/cvernier/H4b/CMSSW_5_3_3_patch2/src/UserCode/SouvikDas/HbbHbb/Analysis/PDF/Archive/LowMass_2x80_4x30_TLTL_minDm_TLTL/KinSel_JERm1
source processKinSel_Signal.sh
cd MMMM_nominal
source processBTagging_Signal.sh
cd a_KinFit
root -l -b processPartitioning_Signal.c
echo "JERp1"
cd /gpfs/ddn/cms/user/cvernier/H4b/CMSSW_5_3_3_patch2/src/UserCode/SouvikDas/HbbHbb/Analysis/PDF/Archive/LowMass_2x80_4x30_TLTL_minDm_TLTL/KinSel_JERp1
source processKinSel_Signal.sh
cd MMMM_nominal
source processBTagging_Signal.sh
cd a_KinFit
root -l -b processPartitioning_Signal.c
echo "upBC"
cd /gpfs/ddn/cms/user/cvernier/H4b/CMSSW_5_3_3_patch2/src/UserCode/SouvikDas/HbbHbb/Analysis/PDF/Archive/LowMass_2x80_4x30_TLTL_minDm_TLTL/KinSel_upBC
source processKinSel_Signal.sh
cd MMMM_upBC
source processBTagging_Signal.sh
cd a_KinFit
root -l -b processPartitioning_Signal.c
echo "downBC"
cd /gpfs/ddn/cms/user/cvernier/H4b/CMSSW_5_3_3_patch2/src/UserCode/SouvikDas/HbbHbb/Analysis/PDF/Archive/LowMass_2x80_4x30_TLTL_minDm_TLTL/KinSel_downBC
source processKinSel_Signal.sh
cd MMMM_downBC
source processBTagging_Signal.sh
cd a_KinFit
root -l -b processPartitioning_Signal.c
echo "upL"
cd /gpfs/ddn/cms/user/cvernier/H4b/CMSSW_5_3_3_patch2/src/UserCode/SouvikDas/HbbHbb/Analysis/PDF/Archive/LowMass_2x80_4x30_TLTL_minDm_TLTL/KinSel_upL
source processKinSel_Signal.sh
cd MMMM_upL
source processBTagging_Signal.sh
cd a_KinFit
root -l -b processPartitioning_Signal.c
echo "downL"
cd /gpfs/ddn/cms/user/cvernier/H4b/CMSSW_5_3_3_patch2/src/UserCode/SouvikDas/HbbHbb/Analysis/PDF/Archive/LowMass_2x80_4x30_TLTL_minDm_TLTL/KinSel_downL
source processKinSel_Signal.sh
cd MMMM_downL
source processBTagging_Signal.sh
cd a_KinFit
root -l -b processPartitioning_Signal.c


