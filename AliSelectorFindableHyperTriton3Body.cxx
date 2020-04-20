#define AliSelectorFindableHyperTriton3Body_cxx

#include "Hyp3FindConfig.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliExternalTrackParam.h"
#include "AliPID.h"
#include "AliVTrack.h"
#include <TCanvas.h>
#include <TDirectory.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TObjArray.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVector3.h>
#include <cmath>

#include <KFVertex.h>
#include <KFParticle.h>

#include "Math/Vector3D.h"

#include <stdio.h>

#include "DCAFitterN.h"

#include "AliSelectorFindableHyperTriton3Body.h"


const float kHypMass = 2.99131;
/// usefull functions

template <typename T> double Sq(T a) { return a * a; }

template <typename T> double Point2PointDistance(T *p0, T *p1) {
  double d2 = 0.;
  for (int iDim = 0; iDim < 3; ++iDim) {
    d2 += Sq(p0[iDim] - p1[iDim]);
  }
  return std::sqrt(d2);
}

template <typename T> double Norm(T x, T y) { return std::sqrt(Sq(x) + Sq(y)); }

template <typename T> double Norm(T x, T y, T z) { return std::sqrt(Sq(x) + Sq(y) + Sq(z)); }

AliSelectorFindableHyperTriton3Body::AliSelectorFindableHyperTriton3Body(TString outputName, TString outputPath, int vertexer,TTree *):
fOutputFileName{outputName},
fOutputFilePath{outputPath},
fAlg{vertexer},
fHypertritonVertexer(){
  fESDtrackCuts = AliESDtrackCuts::GetStandardV0DaughterCuts();
  fESDtrackCuts->SetMinNClustersTPC(0);
  fESDtrackCuts->SetEtaRange(-0.9,0.9);
  //test to compare with the other selector

}

//______________________________________________________________________________
void AliSelectorFindableHyperTriton3Body::Begin(TTree * /*tree*/){
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
}

//______________________________________________________________________________
void AliSelectorFindableHyperTriton3Body::SlaveBegin(TTree * /*tree*/){
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

  const char lAM[3]{"AM"};
  const char *lSpecies[3]{"d", "p", "pi"};
  const char *lCuts[3]{"tree", "selection", "vertexer"};
  const char *lType[3]{"true", "fake", "clones"};
  const char lCoords[3]{'x','y','z'};
  const char *lProj[2]={"pT","ct"};
  const char *lVarName[3]={"NsigmaTPC", "NclusterTPC","NclusterITS"};

  /// histograms for efficiencies
  for(int iMatter=0; iMatter<2; iMatter++){
    fHistCtRec[iMatter] = new TH1D(Form("fHistCtRec%c",lAM[iMatter]),";;",10,0,100);
    fHistX[iMatter] = new TH1D(Form("fHistX%c",lAM[iMatter]),";;",100,-50,50);
    fHistY[iMatter] = new TH1D(Form("fHistY%c",lAM[iMatter]),";;",100,-50,50);
    fHistZ[iMatter] = new TH1D(Form("fHistZ%c",lAM[iMatter]),";;",100,-50,50);

    for(int iVar=0; iVar<kNVar; iVar++){
      fHistFakeVsCuts[iVar][iMatter] = new TH3D(Form("fHistFake_%s_%c",lVarName[iVar],lAM[iMatter]),";#it{p}_{T} (GeV/#it{c});#it{ct} (cm);",kNBinPt,0.,10.,kNBinCt,0.,100.,kNCut,0,kNCut);
      fHistClonesVsCuts[iVar][iMatter] = new TH3D(Form("fHistClones_%s_%c",lVarName[iVar],lAM[iMatter]),";#it{p}_{T} (GeV/#it{c});#it{ct} (cm);",kNBinPt,0.,10.,kNBinCt,0.,100.,kNCut,0,kNCut);
     
      fHistSingleRecVsCuts[iVar][iMatter] = new TH3D(Form("fHistRec_%s_%c",lVarName[iVar],lAM[iMatter]),";#it{p}_{T} (GeV/#it{c});#it{ct} (cm);",kNBinPt,0.,10.,kNBinCt,0.,100.,kNCut,0,kNCut);
      fHistGenVsCuts[iVar][iMatter] = new TH3D(Form("fHistGen_%s_%c",lVarName[iVar],lAM[iMatter]),";#it{p}_{T} (GeV/#it{c});#it{ct} (cm);",kNBinPt,0.,10.,kNBinCt,0.,100.,kNCut,0,kNCut);
      fHistResolutionVsCuts[iVar][iMatter] = new TH3D(Form("fHistRes_%s_%c",lVarName[iVar],lAM[iMatter]),";#Delta#it{p}_{T} (GeV/#it{c});#Delta#it{ct} (cm);",100,-2.,2.,100,-15.,15.,kNCut,0,kNCut);
      
      fHistMassRes[iVar][0][iMatter] = new TH3D(Form("fHistMassRes_%s_%s_%c",lVarName[iVar],lProj[0], lAM[iMatter]), ";#Delta#it{m} (dp#pi) (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});", 100, 2.95-kHypMass, 3.05-kHypMass, 10, 0, 10, kNCut,0,kNCut);
      GetOutputList()->Add(fHistMassRes[iVar][0][iMatter]);
      fHistMassRes[iVar][1][iMatter] = new TH3D(Form("fHistMassRes_%s_%s_%c",lVarName[iVar],lProj[1], lAM[iMatter]), ";#Delta#it{m} (dp#pi) (GeV/#it{c}^{2}); #it{c}t [cm];", 100, 2.95-kHypMass, 3.05-kHypMass, 10, 0, 100, kNCut,0,kNCut);
      GetOutputList()->Add(fHistMassRes[iVar][1][iMatter]);

      GetOutputList()->Add(fHistSingleRecVsCuts[iVar][iMatter]);
      GetOutputList()->Add(fHistGenVsCuts[iVar][iMatter]);    
      GetOutputList()->Add(fHistResolutionVsCuts[iVar][iMatter]); 

      GetOutputList()->Add(fHistFakeVsCuts[iVar][iMatter]);
      GetOutputList()->Add(fHistClonesVsCuts[iVar][iMatter]);
    }
      GetOutputList()->Add(fHistCtRec[iMatter]);
      GetOutputList()->Add(fHistX[iMatter]);
      GetOutputList()->Add(fHistY[iMatter]);
      GetOutputList()->Add(fHistZ[iMatter]);
  }

  /// Histograms for selection

  for(int iCuts = 0; iCuts < 3; iCuts++){
    for(int iMatter = 0; iMatter < 2; iMatter++){
      fHistInvMassPt[iMatter][iCuts] = new TH2D(Form("fHistInvMassPt_%c_%s", lAM[iMatter], lType[iCuts]), ";#it{m} (dp#pi) (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c})", 200, 2.95, 3.35, 20, 0, 10);
      GetOutputList()->Add(fHistInvMassPt[iMatter][iCuts]);
      fHistInvMassPtSel[iMatter][iCuts] = new TH2D(Form("fHistInvMassPt_%c_%s", lAM[iMatter], lCuts[iCuts]), ";#it{m} (dp#pi) (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c})", 200, 2.95, 3.35, 20, 0, 10);
      GetOutputList()->Add(fHistInvMassPtSel[iMatter][iCuts]);
    }
    for(int iSpecies = 0; iSpecies < 3; iSpecies++){
      fHistDaughterPt[iSpecies][iCuts] = new TH1D(Form("fHistPt_%s_%s", lSpecies[iSpecies], lCuts[iCuts]), ";#it{p}_{T} (GeV/#it{c});Counts", 100, 0, 10);
      fHistDaughterTPCchi2[iSpecies][iCuts] = new TH1D(Form("fHistTPCchi2_%s_%s", lSpecies[iSpecies], lCuts[iCuts]), ";#chi^{2};Counts", 100, 0, 10);
      fHistDaughterITSchi2[iSpecies][iCuts] = new TH1D(Form("fHistITSchi2_%s_%s", lSpecies[iSpecies], lCuts[iCuts]), ";#chi^{2};Counts", 100, 0, 10);
      fHistNclsITS[iSpecies][iCuts] = new TH1D(Form("fHistNclsITS_%s_%s", lSpecies[iSpecies], lCuts[iCuts]), ";N_{cluster ITS};Counts", 7, -0.5, 6.5);
      fHistNclsTPC[iSpecies][iCuts] = new TH1D(Form("fHistNclsTPC_%s_%s", lSpecies[iSpecies], lCuts[iCuts]), ";N_{cluster TPC};Counts", 101, -0.5, 200.5);
      fHistNSigmaTPC[iSpecies][iCuts] = new TH1D(Form("fHistNSigmaTPC_%s_%s", lSpecies[iSpecies], lCuts[iCuts]), "n#sigma_{TPC}", 10, 0, 10);
      fHistNSigmaTOF[iSpecies][iCuts] = new TH1D(Form("fHistNSigmaTOF_%s_%s", lSpecies[iSpecies], lCuts[iCuts]), "n#sigma_{TOF}", 10, 0, 10);
      GetOutputList()->Add(fHistDaughterPt[iSpecies][iCuts]);
      GetOutputList()->Add(fHistDaughterTPCchi2[iSpecies][iCuts]);
      GetOutputList()->Add(fHistDaughterITSchi2[iSpecies][iCuts]);
      GetOutputList()->Add(fHistNclsITS[iSpecies][iCuts]);
      GetOutputList()->Add(fHistNclsTPC[iSpecies][iCuts]);
      GetOutputList()->Add(fHistNSigmaTPC[iSpecies][iCuts]);
      GetOutputList()->Add(fHistNSigmaTOF[iSpecies][iCuts]);
    } 
  }
  /// Histograms after vertexer

  fHistVertexChi2 = new TH1D("fHistVertexChi2", "", 100, 0, 200);
  fHistCosPAngle  = new TH1D("fCosPointingAngle", ";#it{cos#theta_{pointing}};Counts", 5000, 0.5, 1.);
  GetOutputList()->Add(fHistVertexChi2);
  GetOutputList()->Add(fHistCosPAngle);

  for (int iCoord = 0; iCoord < 3; iCoord++) {
    fHistResDecayVtx[iCoord] =
        new TH1D(Form("fHistResDecayVtx%c", lCoords[iCoord]), Form(";#Delta%c (mm)", lCoords[iCoord]), 500, -10, 10);
    GetOutputList()->Add(fHistResDecayVtx[iCoord]);
  }

  /// DCA to primary vertex
  for (int iSpecies = 0; iSpecies < 3; iSpecies++) {
    fHistDCA2pvXY[iSpecies] = new TH1D(Form("fDCA2PrimaryvtxXY_%s", lSpecies[iSpecies]), ";DCA_{xy} (mm);Counts", 600, 0, 30);
    fHistDCA2pvZ[iSpecies]  = new TH1D(Form("fDCA2PrimaryvtxZ_%s", lSpecies[iSpecies]), ";DCA_{z} (mm);Counts", 600, 0, 30);
    fHistDCA2pv[iSpecies]   = new TH1D(Form("fDCA2Primaryvtx_%s", lSpecies[iSpecies]), ";DCA (mm);Counts", 600, 0, 30);
    fHistDCA2dvXY[iSpecies] = new TH1D(Form("fDCA2DecayvtxXY_%s", lSpecies[iSpecies]), ";DCA_{xy} (mm);Counts", 600, 0, 30);
    fHistDCA2dvZ[iSpecies]  = new TH1D(Form("fDCA2DecayvtxZ_%s", lSpecies[iSpecies]), ";DCA_{z} (mm);Counts", 600, 0, 30);
    fHistDCA2dv[iSpecies]   = new TH1D(Form("fDCA2Decayvtx_%s", lSpecies[iSpecies]), ";DCA (mm);Counts", 600, 0, 30);
    fHistTrackDistance[iSpecies] = new TH1D(Form("lTrackDistance_%s-%s",lSpecies[iSpecies], lSpecies[(iSpecies + 1) % 3]), ";distance (mm);counts", 200, 0, 200);
    GetOutputList()->Add(fHistDCA2pvXY[iSpecies]);
    GetOutputList()->Add(fHistDCA2pvZ[iSpecies]);
    GetOutputList()->Add(fHistDCA2pv[iSpecies]);
    GetOutputList()->Add(fHistDCA2dvXY[iSpecies]);
    GetOutputList()->Add(fHistDCA2dvZ[iSpecies]);
    GetOutputList()->Add(fHistDCA2dv[iSpecies]);
    GetOutputList()->Add(fHistTrackDistance[iSpecies]);
  }
}

//______________________________________________________________________________
Bool_t AliSelectorFindableHyperTriton3Body::Process(Long64_t entry){
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // When processing keyed objects with PROOF, the object is already loaded
  // and is available via the fObject pointer.
  //
  // This function should contain the \"body\" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.

  fReader.SetEntry(entry);

  // Support variables

  double lPrimaryVtxCov[6] = {0.};
  double lRecPrimaryVtx[3] = {0.};
  double lRecDecayVtx[3] = {0.};
  double lRecDecayLenght[3] = {0.};

  TLorentzVector lLVhyp = {0., 0., 0., 0.};
  TLorentzVector lLVdaughter[3];

  const float lMasses[3]{AliPID::ParticleMass(AliPID::kDeuteron), AliPID::ParticleMass(AliPID::kProton), AliPID::ParticleMass(AliPID::kPion)};

  //------------------------------------------------------------
  // Get main observables from tree
  //------------------------------------------------------------
  float lMagField = *fTreeHyp3BodyVarMagneticField;
  fPrimaryVertex->GetXYZ(lRecPrimaryVtx);
  double lTruePrimaryVtx[3]{*fTreeHyp3BodyVarPVtx[0],*fTreeHyp3BodyVarPVtx[1],*fTreeHyp3BodyVarPVtx[2]};

  AliESDtrack *lTrack[3]{&*fTreeHyp3BodyVarTracks[0],&*fTreeHyp3BodyVarTracks[1],&*fTreeHyp3BodyVarTracks[2]};
  int lDaughterNclsITS[3]{lTrack[0]->GetITSNcls(), lTrack[1]->GetITSNcls(), lTrack[2]->GetITSNcls()};
  int lDaughterNclsTPC[3]{lTrack[0]->GetTPCNcls(), lTrack[1]->GetTPCNcls(), lTrack[2]->GetTPCNcls()};
  double lDaughterChi2TPC[3]{lTrack[0]->GetTPCchi2(), lTrack[1]->GetTPCchi2(), lTrack[2]->GetTPCchi2()};
  double lDaughterChi2ITS[3]{lTrack[0]->GetITSchi2(), lTrack[1]->GetITSchi2(), lTrack[2]->GetITSchi2()};
  float lDaughterNsigmaTPC[3] = {*fTreeHyp3BodyVarNsigmaTPC[0],*fTreeHyp3BodyVarNsigmaTPC[1],*fTreeHyp3BodyVarNsigmaTPC[2]};
  float lDaughterNsigmaTOF[3] = {*fTreeHyp3BodyVarNsigmaTOF[0],*fTreeHyp3BodyVarNsigmaTOF[1],*fTreeHyp3BodyVarNsigmaTOF[2]};

  double lTrueDecayVtx[3]{*fTreeHyp3BodyVarDecayVtx[0],*fTreeHyp3BodyVarDecayVtx[1],*fTreeHyp3BodyVarDecayVtx[2]};
  fCurrentEventId = *fTreeHyp3BodyVarEventId;
  fCurrentMotherId = *fTreeHyp3BodyVarMotherId;
  fFakeCand = *fTreeHyp3BodyVarIsFakeCand;
  bool lCharge = lTrack[0]->GetSign() < 0;

  float lDaughterPt[3] = {0.};
  for (int iTrack = 0; iTrack < 3; iTrack++) {
    lLVdaughter[iTrack].SetXYZM(lTrack[iTrack]->Px(), lTrack[iTrack]->Py(), lTrack[iTrack]->Pz(), lMasses[iTrack]);
    lDaughterPt[iTrack] = lLVdaughter[iTrack].Pt();
    lLVhyp += lLVdaughter[iTrack];
  }

  double lHypMassRec = lLVhyp.M();
  double lHypPtRec   = lLVhyp.Pt();

  //------------------------------------------------------------
  //                  Generated hypertritons                   
  //------------------------------------------------------------

  double lHypPtGen = std::hypot(*fTreeHyp3BodyVarTrueP[0],*fTreeHyp3BodyVarTrueP[1]);
  double lHypPGen = std::hypot(lHypPtGen,*fTreeHyp3BodyVarTrueP[2]);  
  double lHypCtGen = Norm(lTrueDecayVtx[0]-lTruePrimaryVtx[0],lTrueDecayVtx[1]-lTruePrimaryVtx[1],lTrueDecayVtx[2]-lTruePrimaryVtx[2])*kHypMass/lHypPGen;

  if(fCurrentMotherId != fLastMotherId){
    fLastMotherId = fCurrentMotherId;

    for(int iVar=0; iVar<kNVar; iVar++)
      for(int iCut=0; iCut<kNCut; iCut++)
        fHistGenVsCuts[iVar][lCharge]->Fill(lHypPtGen,lHypCtGen,iCut);

    if(fCurrentEventId != fLastEventId) fLastEventId = fCurrentEventId;
    fNclones = 0;
  } else {
    if(fCurrentEventId != fLastEventId){
      fLastEventId = fCurrentEventId;
      
      for(int iVar=0; iVar<kNVar; iVar++)
        for(int iCut=0; iCut<kNCut; iCut++)
          fHistGenVsCuts[iVar][lCharge]->Fill(lHypPtGen,lHypCtGen,iCut);
      
      fNclones = 0;
    }
  }
  

  //------------------------------------------------------------
  //                    Before selection
  //------------------------------------------------------------

  fHistInvMassPtSel[lCharge][0]->Fill(lHypMassRec,lHypPtRec);
  for(int iTrack = 0; iTrack < 3; iTrack++){
    fHistDaughterPt[iTrack][0]->Fill(lDaughterPt[iTrack]);
    fHistNclsITS[iTrack][0]->Fill(lDaughterNclsITS[iTrack]);
    fHistNclsTPC[iTrack][0]->Fill(lDaughterNclsTPC[iTrack]);
    fHistNSigmaTPC[iTrack][0]->Fill(lDaughterNsigmaTPC[iTrack]);
    fHistNSigmaTOF[iTrack][0]->Fill(lDaughterNsigmaTOF[iTrack]);
    fHistDaughterTPCchi2[iTrack][0]->Fill(lDaughterChi2TPC[iTrack]);
    fHistDaughterITSchi2[iTrack][0]->Fill(lDaughterChi2ITS[iTrack]);
  }

  //------------------------------------------------------------
  //                      After selection
  //------------------------------------------------------------
  


  if(AcceptCandidate(0,0)){
    fHistInvMassPtSel[lCharge][1]->Fill(lHypMassRec,lHypPtRec);
    for(int iTrack = 0; iTrack < 3; iTrack++){
      fHistDaughterPt[iTrack][1]->Fill(lDaughterPt[iTrack]);
      fHistNclsITS[iTrack][1]->Fill(lDaughterNclsITS[iTrack]);
      fHistNclsTPC[iTrack][1]->Fill(lDaughterNclsTPC[iTrack]);
      fHistNSigmaTPC[iTrack][1]->Fill(lDaughterNsigmaTPC[iTrack]);
      fHistNSigmaTOF[iTrack][1]->Fill(lDaughterNsigmaTOF[iTrack]);
      fHistDaughterTPCchi2[iTrack][1]->Fill(lDaughterChi2TPC[iTrack]);
      fHistDaughterITSchi2[iTrack][1]->Fill(lDaughterChi2ITS[iTrack]);
    }
  }

  //------------------------------------------------------------
  // Secondary vertex reconstruction
  //------------------------------------------------------------

  // reconstruct the decay vertex with the dedicated vertexer
  RParticles lCollection;
  RHyperTritonO2 recHyp;
  if(fAlg==0){
    if(!KFVertexer(lTrack, lCollection, 0)) return true;
    fHistX[lCharge]->Fill(lCollection.hypertriton.X()-lTrueDecayVtx[0]);
    fHistY[lCharge]->Fill(lCollection.hypertriton.Y()-lTrueDecayVtx[1]);
    fHistZ[lCharge]->Fill(lCollection.hypertriton.Z()-lTrueDecayVtx[2]);
  }
  else if(fAlg==1){
    if(!O2Vertexer(lTrack, recHyp, lRecPrimaryVtx)) return true;
  }
  else if (fAlg==2){
    if(!fHypertritonVertexer.FindDecayVertex(lTrack[0], lTrack[1], lTrack[2], lMagField)) return true;
  }

  if(AcceptCandidate(0,0)){
    //lHypPt is the same for everyone?
    fHistInvMassPtSel[lCharge][2]->Fill(lHypMassRec,lHypPtRec);
    for(int iTrack = 0; iTrack < 3; iTrack++){
      fHistDaughterPt[iTrack][2]->Fill(lDaughterPt[iTrack]);
      fHistNclsITS[iTrack][2]->Fill(lDaughterNclsITS[iTrack]);
      fHistNclsTPC[iTrack][2]->Fill(lDaughterNclsTPC[iTrack]);
      fHistNSigmaTPC[iTrack][2]->Fill(lDaughterNsigmaTPC[iTrack]);
      fHistNSigmaTOF[iTrack][2]->Fill(lDaughterNsigmaTOF[iTrack]);
      fHistDaughterTPCchi2[iTrack][2]->Fill(lDaughterChi2TPC[iTrack]);
      fHistDaughterITSchi2[iTrack][2]->Fill(lDaughterChi2ITS[iTrack]);
    }
  }
  
  float lHypCtRec;
  AliESDVertex *lDecayVertex;
  TVector3 vRecDecayLenght;
  KFPVertex kfPVertex;
  

  if(fAlg==2){  
    lDecayVertex = static_cast<AliESDVertex *>(fHypertritonVertexer.GetCurrentVertex());   
    lDecayVertex->GetXYZ(lRecDecayVtx);
    for (int iCoord = 0; iCoord < 3; iCoord++) {
      fHistResDecayVtx[iCoord]->Fill((lRecDecayVtx[iCoord] - lTrueDecayVtx[iCoord]) * 10.);
      lRecDecayLenght[iCoord] = lRecDecayVtx[iCoord] - lRecPrimaryVtx[iCoord];
    }
    fHistX[lCharge]->Fill(lRecDecayVtx[0]-lTrueDecayVtx[0]);
    fHistY[lCharge]->Fill(lRecDecayVtx[1]-lTrueDecayVtx[1]);
    fHistZ[lCharge]->Fill(lRecDecayVtx[2]-lTrueDecayVtx[2]);
    
    vRecDecayLenght = TVector3(lRecDecayLenght[0], lRecDecayLenght[1], lRecDecayLenght[2]);
    lHypCtRec = vRecDecayLenght.Mag()*kHypMass/lLVhyp.P();
    fHistCtRec[lCharge]->Fill(lHypCtRec);
  }
  else if(fAlg==0){
    kfPVertex.SetXYZ(lRecPrimaryVtx[0],lRecPrimaryVtx[1],lRecPrimaryVtx[2]);
    KFParticle prodVertex{kfPVertex};
    lCollection.hypertriton.SetProductionVertex(prodVertex); 
    lHypCtRec = lCollection.hypertriton.S()*kHypMass;
    fHistCtRec[lCharge]->Fill(lHypCtRec);
  }
  else{
    fHistCtRec[lCharge]->Fill(recHyp.ct);
    lHypMassRec=recHyp.m;
    lHypCtRec=recHyp.ct;
    lHypPtRec=recHyp.pt;
  }
  /// Efficiency histograms
  int fOldClones = fNclones;

  for(int iVar=0; iVar<kNVar; iVar++)
    for(int iCut=0; iCut<kNCut; iCut++)
      if(AcceptCandidate(iCut,iVar)){
        fNclones=fOldClones;
        if(fFakeCand){
          fHistFakeVsCuts[iVar][lCharge]->Fill(lHypPtGen,lHypCtGen,iCut);
          fHistInvMassPt[lCharge][1]->Fill(lHypMassRec,lHypPtRec);
        }
        else{
          if(fNclones==0){
            fHistSingleRecVsCuts[iVar][lCharge]->Fill(lHypPtGen,lHypCtGen,iCut);
            fHistResolutionVsCuts[iVar][lCharge]->Fill(lHypPtGen-lHypPtRec,lHypCtGen-lHypCtRec,iCut);         
            fHistInvMassPt[lCharge][0]->Fill(lHypMassRec,lHypPtRec);
            fNclones++;
          }
          else{
            fHistClonesVsCuts[iVar][lCharge]->Fill(lHypPtGen,lHypCtGen,iCut);
            fHistInvMassPt[lCharge][2]->Fill(lHypMassRec,lHypPtRec);
            fHistMassRes[iVar][0][lCharge]->Fill(lHypMassRec-kHypMass,lHypPtRec,iCut);
            fHistMassRes[iVar][1][lCharge]->Fill(lHypMassRec-kHypMass,lHypCtRec,iCut);
          }
        }
      }

  //only for the standard cuts
  if(!AcceptCandidate(0,0)) return true;

  fHistInvMassPt[lCharge][2]->Fill(lHypMassRec,lHypPtRec);
  //this part is not implemented for the o2
  float cospa;
  if(fAlg==0){
    ROOT::Math::XYZVectorF decayVtx{lCollection.hypertriton.X(), lCollection.hypertriton.Y(), lCollection.hypertriton.Z()};
    ROOT::Math::XYZVectorF mom{lCollection.hypertriton.Px(), lCollection.hypertriton.Py(), lCollection.hypertriton.Pz()};
    cospa = mom.Dot(decayVtx) / std::sqrt(decayVtx.Mag2() * mom.Mag2());

    float dec_vertex[3] = {lCollection.hypertriton.X(), lCollection.hypertriton.Y(), lCollection.hypertriton.Z()};
    double dcaXY[3],dca[3],dcaPP[3];
    dcaXY[0] = lCollection.deuteron.GetDistanceFromVertexXY(dec_vertex);
    dca[0] = lCollection.deuteron.GetDistanceFromVertex(dec_vertex);
    dcaXY[1] = lCollection.proton.GetDistanceFromVertexXY(dec_vertex);
    dca[1] = lCollection.proton.GetDistanceFromVertex(dec_vertex);
    dcaXY[2] = lCollection.pion.GetDistanceFromVertexXY(dec_vertex);
    dca[2] = lCollection.pion.GetDistanceFromVertex(dec_vertex);

    dcaPP[0] = lCollection.deuteron.GetDistanceFromParticle(lCollection.proton);
    dcaPP[1] = lCollection.pion.GetDistanceFromParticle(lCollection.deuteron);
    dcaPP[2] = lCollection.pion.GetDistanceFromParticle(lCollection.proton);
    //still to check the units
    for(int iPar=0; iPar<3; iPar++){
      fHistDCA2dvXY[iPar]->Fill(dcaXY[iPar]);
      fHistDCA2dvZ[iPar]->Fill(TMath::Sqrt(dca[iPar]*dca[iPar]-dcaXY[iPar]*dcaXY[iPar]));
      fHistDCA2dv[iPar]->Fill(dca[iPar]);
      fHistTrackDistance[iPar]->Fill(dcaPP[iPar]);
    }
    fHistVertexChi2->Fill(lCollection.hypertriton.GetChi2()/lCollection.hypertriton.GetNDF());
  }
  
  //this part is implemented only for the std vertexer

  if(fAlg==2){
    float pointAngle = lLVhyp.Angle(vRecDecayLenght);
    cospa = std::cos(pointAngle);
    /// compute the DCA of the 3 tracks from the primary and decay vertex
    AliESDVertex lPV(lTruePrimaryVtx, lPrimaryVtxCov, 1., 1000);
    for (int iTrack = 0; iTrack < 3; iTrack++) {
      double dca2dv[2]    = {0.};
      double dca2pv[2]    = {0.};
      double dca2dvcov[3] = {0.};
      double dca2pvcov[3] = {0.};
      lTrack[iTrack]->PropagateToDCA(lDecayVertex, lMagField, 1000., dca2dv, dca2dvcov);
      lTrack[iTrack]->PropagateToDCA(&lPV, lMagField, 1000., dca2pv, dca2pvcov);

      float dcaXYdv = std::abs(dca2dv[0]) * 10.;    // in mm
      float dcaZdv  = std::abs(dca2dv[1]) * 10.;    // in mm
      float dcadv   = Norm(dcaXYdv, dcaZdv) * 10.;  // in mm

      fHistDCA2dvXY[iTrack]->Fill(dcaXYdv);
      fHistDCA2dvZ[iTrack]->Fill(dcaZdv);
      fHistDCA2dv[iTrack]->Fill(dcadv);

      float dcaXYpv = std::abs(dca2pv[0]) * 10.;    // in mm
      float dcaZpv  = std::abs(dca2pv[1]) * 10.;    // in mm
      float dcapv   = Norm(dcaXYpv, dcaZpv) * 10.;  // in mm

      fHistDCA2pvXY[iTrack]->Fill(dcaXYpv);
      fHistDCA2pvZ[iTrack]->Fill(dcaZpv);
      fHistDCA2pv[iTrack]->Fill(dcapv);
    }
    /// compute the track2track distance used in the vertexer
    float pPM[3][3];
    for (int iPerm = 0; iPerm < 3; iPerm++) {
      fHypertritonVertexer.Find2ProngClosestPoint(lTrack[iPerm], lTrack[(iPerm + 1) % 3], lMagField, pPM[iPerm]);
      fHistTrackDistance[iPerm]->Fill(Point2PointDistance(pPM[iPerm], pPM[(iPerm + 1) % 3]) * 10.);
    }
    double vertexChi2NDF = lDecayVertex->GetChi2perNDF();
    fHistVertexChi2->Fill(vertexChi2NDF);
  }
  
  if(fAlg!=1)
    fHistCosPAngle->Fill(cospa);
  return true;
}

//______________________________________________________________________________
void AliSelectorFindableHyperTriton3Body::SlaveTerminate(){
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
}

//______________________________________________________________________________
void AliSelectorFindableHyperTriton3Body::Terminate(){
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  TFile output(Form("%s/%s", fOutputFilePath.Data(), fOutputFileName.Data()), "RECREATE");
  GetOutputList()->Write();
  output.Close();
}

//______________________________________________________________________________
bool AliSelectorFindableHyperTriton3Body::AcceptCandidate(int iCut = 0, int iVar = 0){
  int CutsSet[3] = {0,0,0};
  CutsSet[iVar] = iCut;


  for(int iTrack = 0; iTrack < 3; iTrack++){
    /*
    bool hasTOF{HasTOF(&*fTreeHyp3BodyVarTracks[iTrack])};
    float dca[2];
    float  fRequireTOFpid[3] = {10.,10.,10.}; /// momentum after which the TOF matching is required
    int   fMinTPCpidClusters[3] = {70, 70, 70};
    float fMinTrackDCA[3] = {0., 0., 0.};
    float fTPCsigmas[3] = {3.5f, 3.5f, 3.5f};
    float fTOFsigmas[3] = {4.f, 4.f, 4.f};
    float fTrackPtRange[3][2] = {{0.f, 7.f},{0.f, 4.f},{0.f, 1.f}};
    constexpr AliPID::EParticleType kAliPID[3]{AliPID::kDeuteron,AliPID::kProton,AliPID::kPion};

    fTreeHyp3BodyVarTracks[iTrack]->GetImpactParameters(dca[0], dca[1]);
    double dcaNorm = std::hypot(dca[0], dca[1]);

    for (int iT{0}; iT < 3; ++iT) {
      nSigmasTPC[iT] = fPIDResponse->NumberOfSigmasTPC(&*fTreeHyp3BodyVarTracks[iTrack], kAliPID[iT]);
      nSigmasTOF[iT] = fPIDResponse->NumberOfSigmasTOF(&*fTreeHyp3BodyVarTracks[iTrack], kAliPID[iT]);
      bool requireTOFpid = &*fTreeHyp3BodyVarTracks[iTrack]->P() > fRequireTOFpid[iT];
      if (std::abs(nSigmasTPC[iT]) < fTPCsigmas[iT] && dcaNorm > fMinTrackDCA[iT] && track->Pt() < fTrackPtRange[iT][1] && 
          track->Pt() > fTrackPtRange[iT][0] && track->GetTPCsignalN() >= fMinTPCpidClusters[iT])
        candidate[iT] = (std::abs(nSigmasTOF[iT]) < fTOFsigmas[iT]) || (!hasTOF && !requireTOFpid);
    }
    */


    if(!fESDtrackCuts->AcceptTrack(&*fTreeHyp3BodyVarTracks[iTrack])) return false;
    //cut on NsigmaTPC
    if(*fTreeHyp3BodyVarNsigmaTPC[iTrack] > kCuts[0][CutsSet[0]][iTrack]) return false;
    //cut on NclusterTPC
    if((*fTreeHyp3BodyVarTracks)->GetTPCNcls() < kCuts[1][CutsSet[1]][iTrack]) return false;
    //cut on NclusterITS
    if((*fTreeHyp3BodyVarTracks)->GetITSNcls() < kCuts[2][CutsSet[2]][iTrack]) return false;
  }
  return true;
}



bool AliSelectorFindableHyperTriton3Body::KFVertexer(AliESDtrack* kTrack [], RParticles &kResult, int Chi2Set){

  float kDeuMass{1.87561};
  float kPMass{0.938272};
  float kPiMass{0.13957};
  float kMasses[3]{kDeuMass,kPMass,kPiMass};
  
  double posmom[6],cov[21];
  RHyperTriton3KFSel recHyp;
  KFParticle helper[3];

  for (int iT=0; iT < 3; iT++) {
    kTrack[iT]->GetXYZ(posmom);
    kTrack[iT]->GetPxPyPz(posmom+3);
    kTrack[iT]->GetCovarianceXYZPxPyPz(cov);
    helper[iT].Create(posmom,cov,kTrack[iT]->Charge(),kMasses[iT]);
    helper[iT].Chi2() = kTrack[iT]->GetTPCchi2();
    helper[iT].NDF() = kTrack[iT]->GetNumberOfTPCClusters() * 2;
  }
  
  KFParticle oneCandidate;  
  //oneCandidate.SetConstructMethod(2);
  oneCandidate.AddDaughter(helper[0]);

  if (kTrack[2] == kTrack[1] || kTrack[1]->Charge() * kTrack[0]->Charge() < 0)
    return false;
  KFParticle twoCandidate{oneCandidate};
  
  twoCandidate.AddDaughter(helper[1]);

  
  recHyp.chi2_deuprot = twoCandidate.GetChi2() / twoCandidate.GetNDF();
  if (recHyp.chi2_deuprot > kMaxKFchi2[0][0])
    return false;

  if (kTrack[1] == kTrack[2] || kTrack[0] == kTrack[2] || kTrack[2]->Charge() * kTrack[1]->Charge() > 0) 
    return false;

  KFParticle hyperTriton{twoCandidate};
  hyperTriton.AddDaughter(helper[2]);
  recHyp.chi2_3prongs = hyperTriton.GetChi2() / hyperTriton.GetNDF();
  if (recHyp.chi2_3prongs > kMaxKFchi2[0][1])
    return false;
  
  recHyp.chi2_topology = hyperTriton.GetChi2() / hyperTriton.GetNDF();

  if (recHyp.chi2_topology > kMaxKFchi2[0][2])
    return false;

  RParticles kCollection;
  kCollection.deuteron = helper[0];
  kCollection.proton = helper[1];
  kCollection.pion = helper[2];
  kCollection.hypertriton = hyperTriton;

  float dca_de_pr = helper[0].GetDistanceFromParticle(helper[1]);
  float dca_de_pi = helper[0].GetDistanceFromParticle(helper[2]);
  float dca_pr_pi = helper[1].GetDistanceFromParticle(helper[2]);
  if(dca_de_pr > 0.05 || dca_de_pi > 0.05 || dca_pr_pi > 0.05)
    return false;
  kResult = kCollection;
  return true;
}

struct HelperParticleO2 {
  o2::track::TrackParCov* track = nullptr;
  float        nSigmaTPC = -1.f;
  float        nSigmaTOF = -1.f;
};

bool AliSelectorFindableHyperTriton3Body::O2Vertexer(AliESDtrack* kTrack [],RHyperTritonO2 &recHyp, double pvPos[]){
  o2::vertexing::DCAFitter3 fVertexer;

  float kDeuMass{1.87561};
  float kPMass{0.938272};
  float kPiMass{0.13957};

  HelperParticleO2 helper[3];
  for (int iT=0; iT < 3; iT++)
    helper[iT].track = static_cast<o2::track::TrackParCov*>((AliExternalTrackParam*)kTrack[iT]);

  if(helper[0].track == helper[1].track || helper[0].track == helper[2].track || helper[1].track == helper[2].track)
    return false;

  int nVert{0};
  try {
    nVert = fVertexer.process(*helper[0].track, *helper[1].track, *helper[2].track);
  }
  catch (std::runtime_error& e) {}//?
  if (nVert) {
    auto vert = fVertexer.getPCACandidate();

    auto& deuTrack = fVertexer.getTrack(0);
    auto& prTrack = fVertexer.getTrack(1);
    auto& piTrack = fVertexer.getTrack(2);

    ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> ldeu{deuTrack.Pt(), deuTrack.Eta(), deuTrack.Phi(), kDeuMass};
    ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> lpro{prTrack.Pt(), prTrack.Eta(), prTrack.Phi(), kPMass};
    ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> lpi{piTrack.Pt(), piTrack.Eta(), piTrack.Phi(), kPiMass};
    ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> hypertriton{ldeu + lpro + lpi};

    const float mass = hypertriton.mass();
    ROOT::Math::XYZVectorF decayVtx{(float)(vert[0] - pvPos[0]), (float)(vert[1] - pvPos[1]), (float)(vert[2] - pvPos[2])};
    const float totalMom = hypertriton.P();
    const float len = std::sqrt(decayVtx.Mag2());

    recHyp.ct = len * kHypMass / totalMom; 
    recHyp.candidates = nVert;
    recHyp.pt = hypertriton.pt();
    recHyp.m = mass;
    return true;
  }
  else 
    return false;
}