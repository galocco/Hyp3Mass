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
  gStyle->SetOptTitle(1);
  
  gSystem->Exec(Form("mkdir %s",folder_name));

  const char  lAM[2] = {'A','M'};
  const char* lProj[2] = {"pT","ct"};
  const char* lProjSet[2] = {"zx","zy"};
  const char* lVarName[3] = {"NsigmaTPC","NclusterTPC","NclusterITS"};
  const char* lMeasure[3] = {"Efficiency","FakeRate","CloneRate"};
  const char* lProjMeas[3] = {"Rec","Fake","Clones"};

  TFile input_file(input_name);
  TFile output_file(Form("%s/%s",folder_name,output_name),"RECREATE");
  const int ndir=5;
  TDirectory* subdir;
  for(int iVar=0; iVar<3; iVar++){
    for(int iMeas=0; iMeas<3; iMeas++){
      for(int iProj=0; iProj<2; iProj++){
        subdir = output_file.mkdir(Form("%s/%s/%s",lVarName[iVar],lMeasure[iMeas],lProj[iProj]));
      }
    }
  }

  TH3D* fHistRecVsCuts = nullptr;
  TH2D* fHistRecProj = nullptr;
  TH1D* fHistRec[3][kNvariations] = {{nullptr}};
  TH2D* fHistGenTot = nullptr;
  TH1D* fHistGen[2] = {nullptr}; //for pt and ct
  TCanvas cv;
  cv.Print(Form("%s/%s[",folder_name,pdf_file));
  //array for the indexes of the cuts, with the cuts ordered from the tighter to the looser
  int ordered_cuts[kNvariations];
  //matter or antimatter
  for(int iMat=0; iMat<2; iMat++){
    // histograms of the genereted hypetritons
    fHistGenTot = (TH2D*) input_file.Get(Form("fHistGen_%c",lAM[iMat]));
    fHistGen[0] = (TH1D*) fHistGenTot->ProjectionX("fHistGenPt",1,fHistGenTot->GetNbinsY());
    fHistGen[1] = (TH1D*) fHistGenTot->ProjectionY("fHistGenCt",1,fHistGenTot->GetNbinsX()); 

    //projection on pt or ct
    for(int iProj=0; iProj<2; iProj++){
      //variable of the cut
      for(int iVar=0; iVar<kNcuts; iVar++){
        //meaurements (efficiency,fake-rate,clone-rate)
        for(int iMeas=0; iMeas<3; iMeas++){
          //get the TH3(pt,ct,cut)
          fHistRecVsCuts = (TH3D*) input_file.Get(Form("fHist%s_%s_%c",lProjMeas[iMeas],lVarName[iVar],lAM[iMat]));
          fHistRecProj = (TH2D*) fHistRecVsCuts->Project3D(lProjSet[iProj]);
          //ordering of the cuts
          ordered_indexes(iVar,ordered_cuts);
          //different cuts
          for(auto iCut : ordered_cuts){
            fHistRec[iMeas][iCut] = (TH1D*)fHistRecProj->ProjectionX(Form("fHist%s_%s_%s_cut_%.1f_%c",lProjMeas[iMeas],lVarName[iVar],lProj[iProj],kCuts[iVar][iCut][0],lAM[iMat]),iCut+1,iCut+1);
            if(lProjMeas[iMeas]=="Rec"){
              fHistRec[iMeas][iCut]->Divide(fHistGen[iProj]);
              SetEfficiencyErrors(fHistRec[iMeas][iCut],fHistGen[iProj]); //TODO: controlllare che faccia cosa giusta
              fHistRec[iMeas][iCut]->GetYaxis()->SetTitle("efficiency");
            } else { //TODO: cambiare questa c
              fHistRec[iMeas][iCut]->Divide(fHistGen[iProj]);
              SetRateErrors(fHistRec[iMeas][iCut],fHistGen[iProj]); 
              fHistRec[iMeas][iCut]->GetYaxis()->SetTitle(Form("#%s/#Gen",lProjMeas[iMeas]));
            }
            //function to set some graphic features
            set_style(fHistRec[iMeas][iCut],iCut);
            fHistRec[iMeas][iCut]->SetTitle(Form("cut on %s : %.1f",lVarName[iVar],kCuts[iVar][iCut][0]));
            output_file.cd(Form("%s/%s/%s",lVarName[iVar],lMeasure[iMeas],lProj[iProj]));
            fHistRec[iMeas][iCut]->Write();
          }
          //set the range of the histograms
          set_range(fHistRec[iMeas]);
          //plot of the histograms in a single canvas
          rainbow_plot(fHist[iMeas],cv,Form("RainbowPlot_%s_%s_%s_%c",lVarName[iVar],lMeasure[iMeas],lProj[iProj],lAM[iMat]),true);          
          output_file.cd();
          cv.Write();          
          cv.Print(Form("%s/%s",folder_name,pdf_file)); 
          cv.Clear();
        }
        //plot of the three measurements in a single canvas
        multirainbow_plot(fHist,cv,Form("MultiRainbowPlot_%s_%s_%c",lVarName[iVar],lProj[iProj],lAM[iMat]),true);
        cv.Print(Form("%s/%s",folder_name,pdf_file));
        cv.Clear();
      }
    }
  }
  cv.Print(Form("%%s/s]",folder_name,pdf_file));
}
///
void compare_efficiencies(char* Std_name="efficiencyStd.root", char* KF_name="efficiencyKF.root", char* O2_name="efficiency02.root", char* output_name="efficiency_comparison.root", char* pdf_file="efficiency_comparison.pdf", char* folder_name="nome")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  
  gSystem->Exec(Form("mkdir %s",folder_name));

  const char lAM[3]{"AM"};
  const char *lProj[2]{"pT", "ct"};
  const char *lMeasure[3]={"Efficiency","FakeRate","CloneRate"};
  const char* lProjMeas[3] = {"Rec","Fake","Clones"};
  const char* lTitle[3] = {"Standard vertexer","0^{2} vertexer","Kalman Filter"};

  TFile* input_file[3];
  input_file[0] = new TFile(Std_name);
  input_file[1] = new TFile(O2_name);
  input_file[2] = new TFile(KF_name);
  TFile output_file(output_name,"RECREATE");
  TH1D* fHist[3][3];

  TCanvas cv;
  cv.Print(Form("%s[",pdf_file));
  for(int iMat=0; iMat<2; iMat++){
    for(int iProj=0; iProj<2; iProj++){
      for(int iMeas=0; iMeas<3; iMeas++){
        for(int iFile=0; iFile<3; iFile++){
          fHist[iMeas][iFile] = (TH1D*) input_file[iFile]->Get(Form("NsigmaTPC/%s/%s/fHist%s_NsigmaTPC_%s_cut_%.1f_%c",lMeasure[iMeas],lProj[iProj],lProjMeas[iMeas],lProj[iProj],kCuts[0][0][0],lAM[iMat]));
          set_style(fHist[iMeas][iFile],iFile);
          fHist[iMeas][iFile]->SetTitle(lTitle[iFile]);
          if(iMeas==0){
            fHist[iMeas][iFile]->GetYaxis()->SetRange(0,1);
          }
        }
        //set_range(fHist[iMeas]);
        rainbow_plot(fHist[iMeas],cv,Form("RainbowPlot_%s_%s_%c",lMeasure[iMeas],lProj[iProj],lAM[iMat]),true);
        output_file.cd();
        cv.Write();          
        cv.Print(Form("%s",pdf_file)); 
        cv.Clear();
      }
      multirainbow_plot(fHist,cv,Form("MultiRainbowPlot_%s_%c",lProj[iProj],lAM[iMat]),true);
      cv.Write();          
      cv.Print(Form("%s",pdf_file)); 
      cv.Clear();    
    }
  }
  cv.Print(Form("%s]",pdf_file));
}

void compare_resolutions(char* Std_name="resolutionStd.root", char* KF_name="resolutionKF.root", char* O2_name="resolution02.root", char* output_name="resolution_comparison.root", char* pdf_file="resolution_comparison.pdf", char* folder_name="nome")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  
  gSystem->Exec(Form("mkdir %s",folder_name));

  const char lAM[3]{"AM"};
  const char *lProj[2]{"pT", "ct"};
  const char *lMeasure[3]={"Efficiency","FakeRate","CloneRate"};
  const char* lProjMeas[3] = {"Rec","Fake","Clones"};
  const char* lTitle[3] = {"Standard vertexer","0^{2} vertexer","Kalman Filter"};

  TFile* input_file[3];
  input_file[0] = new TFile(Std_name);
  input_file[1] = new TFile(O2_name);
  input_file[2] = new TFile(KF_name);
  TFile output_file(output_name,"RECREATE");
  TH1D* fHist[3][3];

  TCanvas cv;
  cv.Print(Form("%s[",pdf_file));
  for(int iMat=0; iMat<2; iMat++){
    for(int iProj=0; iProj<2; iProj++){
      for(int iMeas=0; iMeas<3; iMeas++){
        for(int iFile=0; iFile<3; iFile++){
          fHist[iMeas][iFile] = (TH1D*) input_file[iFile]->Get(Form("NsigmaTPC/%s/%s/fHist%s_NsigmaTPC_%s_cut_%.1f_%c",lMeasure[iMeas],lProj[iProj],lProjMeas[iMeas],lProj[iProj],kCuts[0][0][0],lAM[iMat]));
          set_style(fHist[iMeas][iFile],iFile);
          fHist[iMeas][iFile]->SetTitle(lTitle[iFile]);
          if(iMeas==0)
            fHist[iMeas][iFile]->GetYaxis()->SetRange(0,1);
        }
        //set_range(fHist[iMeas]);
        rainbow_plot(fHist[iMeas],cv,Form("RainbowPlot_%s_%s_%c",lMeasure[iMeas],lProj[iProj],lAM[iMat]),true);
        output_file.cd();
        cv.Write();          
        cv.Print(Form("%s",pdf_file)); 
        cv.Clear();
      }
      multirainbow_plot(fHist,cv,Form("MultiRainbowPlot_%s_%c",lProj[iProj],lAM[iMat]),true);
      cv.Write();          
      cv.Print(Form("%s",pdf_file)); 
      cv.Clear();    
    }
  }
  cv.Print(Form("%s]",pdf_file));
}



void ordered_indexes(int Var,int kCutsOrdered[]){
  float tmp;
  float OriginalOrder[kNvariations];
  float Ordered[kNvariations];

  for(int i=0;i<kNvariations;i++){
    Ordered[i]=kCuts[Var][i][0];
    OriginalOrder[i]=kCuts[Var][i][0];
  }

  for (int j=0;j<kNvariations;j++){
    for (int i=kNvariations-2;i>=j;i--){
      if (Ordered[i]>Ordered[i+1])
      {
        tmp = Ordered[i];
        Ordered[i] = Ordered[i+1];
        Ordered[i+1] = tmp;
      }
    }
  }

  for(int i=0;i<kNvariations;i++){
    for(int j=0;j<kNvariations;j++){
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
