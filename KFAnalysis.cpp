#include <TChain.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <vector>
#include "TH1D.h"
#include "TH2D.h"
#include <TCanvas.h>

void Rainbow_Plot(TH2D* [][5],TH2D* ,const char* [],float [],TFile*);
void KFAnalysis()
{
    TChain mcChain("Hyp3KF");
    mcChain.AddFile("HyperTritonTree3KF.root");

    TTreeReader fReader(&mcChain);
    TTreeReaderArray<RHyperTriton3KF> RHyperVec = {fReader, "RHyperTriton"};
    TTreeReaderArray<SHyperTriton3KF> SHyperVec = {fReader, "SHyperTriton"};
    //TTreeReaderValue<RCollision> RColl = {fReader, "RCollision"};


    const int n_cut_dca = 5;
    const int n_cut_chi2 = 5;
    const char *lCutDCA[3]{"DCA_pd", "DCA_ppi","DCA_dpi"};
    const char *lCutChi2[3]{"Chi2_deuprot","Chi2_3prongs","Chi2_topology"};
    
    TH1D* fHistChi[3];
    TH1D* fHistDCA[3];
    
    for(int iCut=0; iCut<3; iCut++)
    {
      fHistChi[iCut] = new TH1D(Form("fHist%s",lCutChi2[iCut]),";#Chi^2",200,0,10);
      fHistDCA[iCut] = new TH1D(Form("fHist%s",lCutDCA[iCut]),";DCA",200,0,2.2);
    }

    TH2D* fHistGen = new TH2D("fHistGen",";p_{T} [GeV/c];ct [cm]",10,0,10,10,0,100);
    TH2D* fHistEffDCA[3][n_cut_dca];
    TH2D* fHistEffChi2[3][n_cut_chi2];

    float cut_dca[]={0.,0.1,0.2,0.3,0.4};
    float cut_chi2[]={0.,0.1,0.2,0.3,0.4};
    //inizializzo gli istogrammi
    for(int iCut=0; iCut<n_cut_dca; iCut++)
    {
      for(int iPar=0; iPar<3;iPar++)
        fHistEffDCA[iPar][iCut] = new TH2D(Form("fHistEff_%.1f_%s",cut_dca[iCut],lCutDCA[iPar]),";p_{T} [GeV/c];ct [cm]",10,0,10,10,0,100);
    }

    for(int iCut=0; iCut<n_cut_chi2; iCut++)
    {
      for(int iPar=0; iPar<3;iPar++)
        fHistEffChi2[iPar][iCut] = new TH2D(Form("fHistEff_%.1f_%s",cut_chi2[iCut],lCutChi2[iPar]),";p_{T} [GeV/c];ct [cm]",10,0,10,10,0,100);
    }
    //leggo i file e riempio gli istogrammi
    while (fReader.Next())
    {
      for (auto &SHyper : SHyperVec)
      {
        fHistGen->Fill(TMath::Sqrt(SHyper.px*SHyper.px+SHyper.px*SHyper.px),SHyper.l);
      }

      for (auto &RHyper : RHyperVec)
      {

        double dca[3]={RHyper.dca_de_pr,RHyper.dca_de_pi,RHyper.dca_pr_pi};
        double chi2[3]={RHyper.chi2_deuprot,RHyper.chi2_3prongs,RHyper.chi2_topology}; 

        for(int iCut=0; iCut<3; iCut++)
        {
          fHistChi[iCut]->Fill(chi2[iCut]);
          fHistDCA[iCut]->Fill(dca[iCut]);
        }

        for(int cut=0; cut<n_cut_dca; cut++)
        {
          for(int iPar=0; iPar<3;iPar++)
          {
            if(dca[iPar]<cut_dca[cut])
              fHistEffDCA[iPar][cut]->Fill(TMath::Sqrt(RHyper.px*RHyper.px+RHyper.px*RHyper.px),RHyper.l);
          }
        }
        
        for(int cut=0; cut<n_cut_chi2; cut++)
        {
          for(int iPar=0; iPar<3;iPar++)
          {
            if(chi2[iPar]<cut_dca[cut])
              fHistEffChi2[iPar][cut]->Fill(TMath::Sqrt(RHyper.px*RHyper.px+RHyper.px*RHyper.px),RHyper.l);
          }
        }
      }

    }


  TFile* rfile = new TFile("KFResults.root","RECREATE");
  rfile->Write();
  Rainbow_Plot(fHistEffChi2,fHistGen,lCutChi2,cut_chi2,rfile);
  Rainbow_Plot(fHistEffDCA,fHistGen,lCutDCA,cut_dca,rfile);
    
}

void Rainbow_Plot(TH2D* fHistRec[][5],TH2D* fHistGen,const char* lVar[],float cut[],TFile* rfile)
{

  rfile->ReOpen("UPDATE");
  //file where results will be saved
  
  //proiezioni e divisioni
  TH1D* fHistEff;
  TH1D* fHistGenProj;
  TCanvas RainbowPt("","");
  TCanvas RainbowCt("","");

  //faccio le proiezioni e calcolo le efficienze
  for(int iPar=0; iPar<3;iPar++)
  {
    for(int iCut=0; iCut<5; iCut++)
    {
      //ct
      fHistEff = (TH1D*) fHistRec[iPar][iCut]->ProjectionY(Form("fHistEff_Ct_%s_%.1f",lVar[iPar],cut[iCut]),1,fHistRec[iPar][iCut]->GetNbinsY());  
      fHistGenProj = (TH1D*) fHistGen->ProjectionY("fHistGenProj",1,fHistGen->GetNbinsY());  

      fHistEff->SetMarkerStyle(8);
      fHistEff->Divide(fHistGenProj);
      fHistEff->GetYaxis()->SetTitle("efficiency");
      rfile->cd();
      fHistEff->Write();
      RainbowCt.cd();
      if(iCut==0)
        fHistEff->Draw("PMC PLC");
      else
        fHistEff->Draw("SAME PMC PLC");

      //pt
      fHistEff = (TH1D*) fHistRec[iPar][iCut]->ProjectionX(Form("fHistEff_Pt_%s_%.1f",lVar[iPar],cut[iCut]),1,fHistRec[iPar][iCut]->GetNbinsX());  
      fHistGenProj = (TH1D*) fHistGen->ProjectionX("fHistGenProj",1,fHistGen->GetNbinsX());  

      fHistEff->SetMarkerStyle(8);
      fHistEff->Divide(fHistGenProj);
      fHistEff->GetYaxis()->SetTitle("efficiency");
      rfile->cd();
      fHistEff->Write();
      RainbowPt.cd();
      if(iCut==0)
        fHistEff->Draw("PMC PLC");
      else
        fHistEff->Draw("SAME PMC PLC");

    }
    //dca/chi2_taglio
    RainbowCt.SetName(Form("RainbowPlotCt_%s",lVar[iPar])); 
    RainbowCt.BuildLegend();

    RainbowPt.SetName(Form("RainbowPlotPt_%s",lVar[iPar])); 
    RainbowPt.BuildLegend();

    rfile->cd();

    RainbowPt.Write();
    RainbowCt.Write();
    
  }

}