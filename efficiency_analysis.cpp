#include <TCanvas.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <string.h>
#include "Hyp3FindConfig.h"
#include "utils.h"
void ordered_indexes(int ,int []);

void selector_efficiencies(char * input_name = "selector_results.root",char * output_name = "efficiency.root",char* pdf_file = "efficiency_plots.pdf",char * folder_name = "plot_folder")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  gSystem->Exec(Form("mkdir %s",folder_name));

  const char lAM[3]{"AM"};
  const char *lPrintMode[3]={"[","","]"};
  const char *lProj[2]{"pT", "ct"};
  const char *lProjSet[2]{"zx", "zy"};
  const char *lVarName[3]={ "NsigmaTPC","NclusterTPC","NclusterITS"};
  const char *lMeasure[3]={"Efficiency","FakeRate","CloneRate"};
  const char* lProjMeas[3] = {"Rec","Fake","Clones"};

  TFile input_file(input_name);
  TFile output_file(output_name,"RECREATE");
  const int ndir=5;
  TDirectory* subdir;
  for(int iVar=0; iVar<3; iVar++){
    for(int iMeas=0; iMeas<4; iMeas++)
      for(int iProj=0; iProj<2; iProj++)
        subdir = output_file.mkdir(Form("%s/%s/%s",lVarName[iVar],lMeasure[iMeas],lProj[iProj]));
  }

  TH3D* fHistVsCuts = nullptr;
  TH2D* fHistProj = nullptr;
  TH1D* fHist[3][number_of_variations] = {{nullptr}};
  TH2D* fHistGenTot = nullptr;
  TH1D* fHistGen[2] = {nullptr};//for pt and ct
  TCanvas cv;
  cv.Print(Form("%s[",pdf_file));
  //array for the indexes of the cuts, with the cuts ordered from the tighter to the looser
  int ordered_cuts[number_of_variations];
  //matter or antimatter
  for(int iMat=0; iMat<2; iMat++){
    // histograms of the genereted hypetritons
    fHistGenTot = (TH2D*) input_file.Get(Form("fHistGen_%c",lAM[iMat]));
    fHistGen[0] = (TH1D*) fHistGenTot->ProjectionX("fHistGenPt",1,fHistGenTot->GetNbinsY());
    fHistGen[1] = (TH1D*) fHistGenTot->ProjectionY("fHistGenCt",1,fHistGenTot->GetNbinsX()); 

    //projection on pt or ct
    for(int iProj=0; iProj<2; iProj++){
      //variable of the cut
      for(int iVar=0; iVar<number_of_cuts; iVar++){
        //meaurements (efficiency,fake-rate,clone-rate)
        for(int iMeas=0; iMeas<3; iMeas++){
          //get the TH3(pt,ct,cut)
          fHistVsCuts = (TH3D*) input_file.Get(Form("fHist%s_%s_%c",lProjMeas[iMeas],lVarName[iVar],lAM[iMat]));
          fHistProj = (TH2D*) fHistVsCuts->Project3D(lProjSet[iProj]);
          //ordering of the cuts
          ordered_indexes(iVar,ordered_cuts);
          //different cuts
          for(auto iCut :ordered_cuts){
            fHist[iMeas][iCut] = (TH1D*)fHistProj->ProjectionX(Form("fHist%s_%s_%s_cut_%.1f_%c",lProjMeas[iMeas],lVarName[iVar],lProj[iProj],kCuts[iVar][iCut][0],lAM[iMat]),iCut+1,iCut+1);
            if(lProjMeas[iMeas]=="Rec"){
              fHist[iMeas][iCut]->Divide(fHistGen[iProj]);
              SetEfficiencyErrors(fHist[iMeas][iCut],fHistGen[iProj]);
              fHist[iMeas][iCut]->GetYaxis()->SetTitle("efficiency");
            }else if(lProjMeas[iMeas]=="Fake" || lProjMeas[iMeas]=="Clones"){
              fHist[iMeas][iCut]->Divide(fHistGen[iProj]);
              SetRateErrors(fHist[iMeas][iCut],fHistGen[iProj]); 
              fHist[iMeas][iCut]->GetYaxis()->SetTitle(Form("#%s/#Gen",lProjMeas[iMeas]));
            }
            else{
              fHist[iMeas][iCut]->Scale(1./fHistGen[iProj]->GetEntries());                
              fHist[iMeas][iCut]->GetYaxis()->SetTitle("counts/N_{ev}");
            }
            //function to set some graphic features
            set_style(fHist[iMeas][iCut],iCut);
            fHist[iMeas][iCut]->SetTitle(Form("cut on %s : %.1f",lVarName[iVar],kCuts[iVar][iCut][0]));
            output_file.cd(Form("%s/%s/%s",lVarName[iVar],lMeasure[iMeas],lProj[iProj]));
            fHist[iMeas][iCut]->Write();
          }
          //set the range of the histograms
          set_range(fHist[iMeas]);
          //plot of the histograms in a single canvas
          rainbow_plot(fHist[iMeas],cv,Form("RainbowPlot_%s_%s_%s_%c",lVarName[iVar],lMeasure[iMeas],lProj[iProj],lAM[iMat]),true);          
          output_file.cd();
          cv.Write();
          std::cout<<Form("RainbowPlot_%s_%s_%s_%c",lVarName[iVar],lMeasure[iMeas],lProj[iProj],lAM[iMat])<<" ";
          
          cv.Print(Form("%s",pdf_file)); 
          cv.Clear();
        }
        //plot of the three measurements in a single canvas
        multirainbow_plot(fHist,cv,Form("MultiRainbowPlot_%s_%s_%c",lVarName[iVar],lProj[iProj],lAM[iMat]),true);
        cv.Print(Form("%s",pdf_file));
        cv.Clear();
      }
    }
  }
  cv.Print(Form("%s]",pdf_file));
}
/*
void compare_efficiency(char * std_name = "efficiency.root",char * kf_name = "efficiencyKF.root",char * output_name = "comparison.root")
{
  gSystem->Exec(Form("mkdir comparison_images"));
  gStyle->SetOptStat(0);
  TFile std_file(std_name);
  TFile kf_file(kf_name);
  //TODO:  hardcoded names
  TFile stdmass_file("selector_resultsStd.root");
  TFile kfmass_file("selector_resultsKF.root");

  const char lAM[3]{"AM"};
  const char *lProj[2]{"pT", "ct"};
  const char *lProjSet[2]{"zx", "zy"};
  const char *lVarName[3]={"NsigmaTPC", "NclusterTPC","NclusterITS"};
  const char *lMeasure[4]={"Efficiency","Resolution","FakeRate","CloneRate"};
  //Gen must be the first of the list
  const char* lProjMeas[5] = {"Gen","Rec","Res","Fake","Clones"};
  const char *lType[3]{"true", "fake", "clones"};

  TFile output_file(output_name,"RECREATE");
  output_file.cd();

  TCanvas cv("","");

  TH1D* fHistkf;
  TH1D* fHiststd;
  TH2D* fHistInvMass;
  //matter or antimatter
  for(int iMatter=0; iMatter<2; iMatter++){
    output_file.cd();
    for(int iType=0; iType<3; iType++){
      fHistInvMass = (TH2D*) stdmass_file.Get(Form("fHistInvMassPt_%c_%s",lAM[iMatter],lType[iType]));
      fHistInvMass->SetName(Form("STDfHistInvMassPt_%s_%c",lType[iType],lAM[iMatter]));
      fHistInvMass->Write();
      fHistInvMass = (TH2D*) kfmass_file.Get(Form("fHistInvMassPt_%c_%s",lAM[iMatter],lType[iType]));
      fHistInvMass->SetName(Form("KFfHistInvMassPt_%s_%c",lType[iType],lAM[iMatter]));
      fHistInvMass->Write();
    }
    TLine zeroline(0,0,0,1);
    //projection on pt or ct
    for(int iProj=0; iProj<2; iProj++){
      //measurements
      for(int iHist=1; iHist<5; iHist++){
        //if(iHist==2) continue;
        fHiststd = (TH1D*) std_file.Get(Form("fHist%s_%s_%s_cut_%.1f_%c",lProjMeas[iHist],lVarName[0],lProj[iProj],kCuts[0][0][0],lAM[iMatter]));
        fHistkf = (TH1D*) kf_file.Get(Form("fHist%s_%s_%s_cut_%.1f_%c",lProjMeas[iHist],lVarName[0],lProj[iProj],kCuts[0][0][0],lAM[iMatter]));
        if(lMeasure[iHist-1]=="Efficiency"){
          fHiststd->GetYaxis()->SetRange(0.,1.);
          fHistkf->GetYaxis()->SetRange(0.,1.);
        }
        fHiststd->SetTitle("std vertexer");
        fHistkf->SetTitle("kf vertexer");
        cv.cd();
        fHistkf->Draw("PMC PLC");
        fHiststd->Draw("SAME PMC PLC");
        cv.SetName(Form("Comparison_%s_%s_%c",lMeasure[iHist-1],lProj[iProj],lAM[iMatter])); 
        cv.BuildLegend();
        if(lMeasure[iHist-1]=="Resolution"){
          float max1=fHiststd->GetBinContent(fHiststd->GetMaximumBin());
          float max2=fHistkf->GetBinContent(fHistkf->GetMaximumBin());
          zeroline.SetY2(TMath::Max(max1+0.1,max2+0.1));
          zeroline.SetLineColor(kRed);
          zeroline.SetLineStyle(2);
          zeroline.Draw("SAME");
        }        
        
        cv.SaveAs(Form("comparison_images/comparison_%s_%s_%c.pdf",lMeasure[iHist-1],lProj[iProj],lAM[iMatter])); 

        output_file.cd();
        cv.Write();
      }  
    }
  }
}
*/
void ordered_indexes(int Var,int kCutsOrdered[]){
  float tmp;
  float OriginalOrder[number_of_variations];
  float Ordered[number_of_variations];

  for(int i=0;i<number_of_variations;i++){
    Ordered[i]=kCuts[Var][i][0];
    OriginalOrder[i]=kCuts[Var][i][0];
  }

  for (int j=0;j<number_of_variations;j++){
    for (int i=number_of_variations-2;i>=j;i--){
      if (Ordered[i]>Ordered[i+1])
      {
        tmp = Ordered[i];
        Ordered[i] = Ordered[i+1];
        Ordered[i+1] = tmp;
      }
    }
  }

  for(int i=0;i<number_of_variations;i++){
    for(int j=0;j<number_of_variations;j++){
      if(Ordered[i]==OriginalOrder[j])
        kCutsOrdered[i]=j;
    }
  }
}
/*
//function to compare efficiencies
void compare_system(char * findable_name = "efficiencyKF.root",char * task_name = "KFResults.root",char * output_name = "compare_method.root")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TFile findable_file(findable_name);
  TFile task_file(task_name);
  const char lAM[3]{"AM"};
  const char* lName[2]={"Antimatter","Matter"};
  TH1D* histos[2][2];
  for(int iMat=0;iMat<2;iMat++){
    histos[0][0] = (TH1D*) task_file.Get(Form("efficiency_chi2/fHistEff_Pt_Chi2_deuprot_50.0_%c",lAM[iMat]));
    histos[0][1] = (TH1D*) findable_file.Get(Form("NsigmaTPC/Efficiency/pT/fHistRec_NsigmaTPC_pT_cut_3.0_%c",lAM[iMat]));
    histos[1][0] = (TH1D*) task_file.Get(Form("efficiency_chi2/fHistEff_Ct_Chi2_deuprot_50.0_%c",lAM[iMat]));
    histos[1][1] = (TH1D*) findable_file.Get(Form("NsigmaTPC/Efficiency/ct/fHistRec_NsigmaTPC_ct_cut_3.0_%c",lAM[iMat]));


    histos[0][0]->SetTitle("Task 3KF");
    histos[0][1]->SetTitle("Findable");
    multirainbow_plot(histos,output_name,lName[iMat],true,iMat==0);
  }
}
*/