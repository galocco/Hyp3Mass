#include <TCanvas.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <string.h>
#include "Hyp3FindConfig.h"
#include "utils.h"
#include "TString.h"
void ordered_indexes(int ,int []);

void selector_efficiencies(TString input_name = "selector_results.root",TString output_name = "efficiency.root",TString pdf_file = "efficiency_plots.pdf",TString folder_name = "plot_folder")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  gSystem->Exec(Form("mkdir %s",folder_name.Data()));

  const char  lAM[2] = {'A','M'};
  const char* lProj[2] = {"pT","ct"};
  const char* lProjSet[2] = {"zx","zy"};
  const char* lVarName[3] = {"NsigmaTPC","NclusterTPC","NclusterITS"};
  const char* lMeasure[3] = {"Efficiency","FakeRate","CloneRate"};
  const char* lProjMeas[3] = {"Rec","Fake","Clones"};

  TFile input_file(input_name.Data());
  TFile gen_file("../AnalysisResults.root");
  TList* list_gen = (TList*) gen_file.Get("Hyp3FindTask_summary");
  TH3D* fHistGen3D[2];
  fHistGen3D[0] = (TH3D*) list_gen->FindObject("fHistGeneratedPtVsCtVsCentralityHypTrit3");
  fHistGen3D[1] = (TH3D*) list_gen->FindObject("fHistGeneratedPtVsCtVsCentralityAntiHypTrit3");
  //get the blastwave and normalize it
  TFile bw_file(bw_name.Data());
  TF1* blast_wave = (TF1*) bw_file.Get(blast_wave_name);
  float max_bw = blast_wave->GetMaximum(0,10);

  TFile output_file(Form("%s/%s",folder_name.Data(),output_name.Data()),"RECREATE");
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
  TCanvas cv("","",800,450);
  cv.Print(Form("%s/%s[",folder_name.Data(),pdf_file.Data()));
  //array for the indexes of the cuts, with the cuts ordered from the tighter to the looser
  int ordered_cuts[kNvariations];
  //matter or antimatter
  for(int iMat=0; iMat<2; iMat++){
    // histograms of the genereted hypetritons
    fHistGen[0] = (TH1D*)fHistGen3D[iMat]->ProjectionX();//pt
    fHistGen[1] = (TH1D*)fHistGen3D[iMat]->ProjectionY();//ct
    //the gen hyps are flat in pT -> they must be scale with the true distribution
    fHistGen[0]->Multiply(blast_wave,1./max_bw);//pt
    fHistGen[1]->Scale(blast_wave->Integral(0,10)/max_bw);//ct

    fHistGen[0]->Draw();
    cv.Print(Form("%s/%s",folder_name.Data(),pdf_file.Data())); 
    fHistGen[1]->Draw();
    cv.Print(Form("%s/%s",folder_name.Data(),pdf_file.Data())); 
    TString rec_string = "Rec";
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
            if(lProjMeas[iMeas]==rec_string.Data()){
              fHistRec[iMeas][iCut]->Divide(fHistGen[iProj]);
              SetEfficiencyErrors(fHistRec[iMeas][iCut],fHistGen[iProj]);
              fHistRec[iMeas][iCut]->GetYaxis()->SetTitle("efficiency");
            } else {
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
          rainbow_plot(fHistRec[iMeas],cv,Form("RainbowPlot_%s_%s_%s_%c",lVarName[iVar],lMeasure[iMeas],lProj[iProj],lAM[iMat]),true);          
          output_file.cd();
          cv.Write();          
          cv.Print(Form("%s/%s",folder_name.Data(),pdf_file.Data())); 
          cv.Clear();
        }
        //plot of the three measurements in a single canvas
        multirainbow_plot(fHistRec,cv,Form("MultiRainbowPlot_%s_%s_%c",lVarName[iVar],lProj[iProj],lAM[iMat]),true);
        cv.Print(Form("%s/%s",folder_name.Data(),pdf_file.Data()));
        cv.Clear();
      }
    }
  }
  cv.Print(Form("%s/%s]",folder_name.Data(),pdf_file.Data()));
}


///
void compare_efficiencies(TString Std_name="plot_folderStd/efficiencyStd.root", TString KF_name="plot_folderKF/efficiencyKF.root", TString O2_name="plot_folderO2/efficiencyO2.root", TString output_name="efficiency_comparison.root", TString pdf_file="efficiency_comparison.pdf", TString folder_name="eff_comparsion")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  gSystem->Exec(Form("mkdir %s",folder_name.Data()));

  const char  lAM[2] = {'A','M'};
  const char* lProj[2] = {"pT","ct"};
  const char* lMeasure[3] = {"Efficiency","FakeRate","CloneRate"};
  const char* lProjMeas[3] = {"Rec","Fake","Clones"};
  const char* lTitle[3] = {"Standard vertexer","O^{2} vertexer","Kalman Filter"};

  TFile* input_file[3];
  input_file[0] = new TFile(Std_name.Data());
  input_file[1] = new TFile(O2_name.Data());
  input_file[2] = new TFile(KF_name.Data());
  TFile output_file(output_name.Data(),"RECREATE");
  TH1D* fHist[3][3];//[measurement(eff,fake,clones)][vertexer]

  TCanvas cv("","",800,450);
  cv.Print(Form("%s[",pdf_file.Data()));
  //antimatter or matter
  for(int iMat=0; iMat<2; iMat++){
    //pt or ct
    for(int iProj=0; iProj<2; iProj++){
      //compare efficiency, fake rate and clones rate
      for(int iMeas=0; iMeas<3; iMeas++){
        //loop on the three vertexer results
        for(int iFile=0; iFile<3; iFile++){
          fHist[iMeas][iFile] = (TH1D*) input_file[iFile]->Get(Form("NsigmaTPC/%s/%s/fHist%s_NsigmaTPC_%s_cut_%.1f_%c",lMeasure[iMeas],lProj[iProj],lProjMeas[iMeas],lProj[iProj],kCuts[0][0][0],lAM[iMat]));
          set_style(fHist[iMeas][iFile],iFile);
          fHist[iMeas][iFile]->SetTitle(lTitle[iFile]);
        }
        //set y axis range
        set_range(fHist[iMeas]);
        //plots the all the iMeas-measurements in the same canvas
        rainbow_plot(fHist[iMeas],cv,Form("RainbowPlot_%s_%s_%c",lMeasure[iMeas],lProj[iProj],lAM[iMat]),true);
        
        output_file.cd();
        cv.Write();          
        cv.Print(Form("%s",pdf_file.Data())); 
        cv.Clear();
      }
      //plots the all the measurements in the same canvas
      multirainbow_plot(fHist,cv,Form("MultiRainbowPlot_%s_%c",lProj[iProj],lAM[iMat]),true);
      cv.Write();          
      cv.Print(Form("%s",pdf_file.Data())); 
      cv.Clear();    
    }
  }
  cv.Print(Form("%s]",pdf_file.Data()));
}

void compare_resolutions(TString Std_name="MassResolutionStd.root", TString KF_name="MassResolutionKF.root", TString O2_name="MassResolutionO2.root", TString output_name="resolution_comparison.root", TString pdf_file="resolution_comparison.pdf", TString folder_name="resolutionAll")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadRightMargin(0.1); 
  gStyle->SetPadLeftMargin(0.15);
  
  gSystem->Exec(Form("mkdir %s",folder_name.Data()));

  const char  lAM[2] = {'A','M'};
  const char* lTitle[3] = {"Standard vertexer","O^{2} vertexer","Kalman Filter #chi^{2}/NDF<50"};
  //const char* lTitle[3] = {"No cuts","#chi^{2} <0.5","#chi^{2}<5"};
  
  const char* lMeasure[2] = {"Average","StdDev"};
  const char* lRes[8] = {"Mass_pT","Mass_ct","Pt_pT","P_p","Ct_ct","X_X","Y_Y","Z_Z"};
  const char* lDir[2] = {"mean_value","stddev"};

  TFile* input_file[3];
  input_file[0] = new TFile(Std_name.Data());
  input_file[1] = new TFile(O2_name.Data());
  input_file[2] = new TFile(KF_name.Data());
  TFile output_file(Form("%s/%s",folder_name.Data(),output_name.Data()),"RECREATE");
  TH1D* fHist[2][3];//[mean/stddev][vertexer]

  TCanvas cv("","",800,450);
  cv.Print(Form("%s/%s[",folder_name.Data(),pdf_file.Data()));
  //antimatter and matter
  for(int iMat=0; iMat<2; iMat++){
    //all the resolutions 
    for(int iRes=0; iRes<8; iRes++){
      //loop on the measurements: mean and stddev
      for(int iMeas=0; iMeas<2; iMeas++){          
        //Kalman filter, O2 ,standard vertexer
        for(int iFile=0; iFile<3; iFile++){
          std::cout<<Form("%s/f%s_%s_%c",lDir[0],lMeasure[0],lRes[1],lAM[0])<<"\n";
          fHist[iMeas][iFile] = (TH1D*) input_file[iFile]->Get(Form("%s/f%s_%s_%c",lDir[iMeas],lMeasure[iMeas],lRes[iRes],lAM[iMat]));
          set_style(fHist[iMeas][iFile],iFile);
          fHist[iMeas][iFile]->SetTitle(lTitle[iFile]);
        }
        set_range(fHist[iMeas],5);
        rainbow_plot(fHist[iMeas],cv,Form("RainbowPlot_%s_%s_%c",lMeasure[iMeas],lRes[iRes],lAM[iMat]),true);
        output_file.cd();
        cv.Write();          
        cv.Print(Form("%s/%s",folder_name.Data(),pdf_file.Data())); 
        cv.Clear();
      }
      multirainbow_plot(fHist,cv,Form("MultiRainbowPlot_%s_%c",lRes[iRes],lAM[iMat]),true);
      cv.Write();          
      cv.Print(Form("%s/%s",folder_name.Data(),pdf_file.Data())); 
      cv.Clear();    
    }
  }
  cv.Print(Form("%s/%s]",folder_name.Data(),pdf_file.Data()));
}
void compare_chi2(TString KF_name1="MassResolutionKFChi2_1.root",TString KF_name2="MassResolutionKFChi2_2.root",TString KF_name3="MassResolutionKFChi2_5.root",TString KF_name4="MassResolutionKFChi2_50.root",TString output_name="chi2_comparison.root", TString pdf_file="chi2_comparison.pdf", TString folder_name="chi2")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadRightMargin(0.1); 
  gStyle->SetPadLeftMargin(0.15);
  gSystem->Exec(Form("mkdir %s",folder_name.Data()));

  const char  lAM[2] = {'A','M'};
  const char* lTitle[4] = {"#chi^{2}/NDF<1","#chi^{2}/NDF<2","#chi^{2}/NDF<5","#chi^{2}/NDF<50"};
  //const char* lTitle[3] = {"No cuts","#chi^{2} <0.5","#chi^{2}<5"};
  
  const char* lMeasure[2] = {"Average","StdDev"};
  const char* lRes[8] = {"Mass_pT","Mass_ct","Pt_pT","P_p","Ct_ct","X_X","Y_Y","Z_Z"};
  const char* lDir[2] = {"mean_value","stddev"};

  TFile* input_file[4];
  input_file[0] = new TFile(KF_name1.Data());
  input_file[1] = new TFile(KF_name2.Data());
  input_file[2] = new TFile(KF_name3.Data());
  input_file[3] = new TFile(KF_name4.Data());
  TFile output_file(output_name,"RECREATE");
  TH1D* fHist[2][4];//[mean/stddev][cuts]

  TLegend leg(0.2,0.2,0.4,0.4);
  leg.SetHeader("Kalman filter","C"); // option "C" allows to center the header
  for(int iFile=0; iFile<4; iFile++){
    TH1D* help = new TH1D("help","",10,0,10);
    set_style(help,iFile);
    leg.AddEntry(help,lTitle[iFile],"lp");  
  } 

  TCanvas cv("","",800,450);
  cv.Print(Form("%s[",pdf_file.Data()));
  //antimatter and matter
  for(int iMat=0; iMat<2; iMat++){
    //all the resolutions 
    for(int iRes=0; iRes<8; iRes++){
      //loop on the measurements: mean and stddev
      for(int iMeas=0; iMeas<2; iMeas++){          
        //Kalman filter, O2 ,standard vertexer
        for(int iFile=0; iFile<4; iFile++){
          std::cout<<Form("%s/f%s_%s_%c",lDir[0],lMeasure[0],lRes[1],lAM[0])<<"\n";
          fHist[iMeas][iFile] = (TH1D*) input_file[iFile]->Get(Form("%s/f%s_%s_%c",lDir[iMeas],lMeasure[iMeas],lRes[iRes],lAM[iMat]));
          set_style(fHist[iMeas][iFile],iFile);
          fHist[iMeas][iFile]->SetTitle(lTitle[iFile]);
        }
        set_range(fHist[iMeas],5);
        rainbow_plot(fHist[iMeas],cv,Form("RainbowPlot_%s_%s_%c",lMeasure[iMeas],lRes[iRes],lAM[iMat]),false);
        cv.cd();
        leg.Draw("SAME");
        output_file.cd();
        cv.Write();          
        cv.Print(Form("%s",pdf_file.Data())); 
        cv.Clear();
      }
      multirainbow_plot(fHist,cv,Form("MultiRainbowPlot_%s_%c",lRes[iRes],lAM[iMat]),false);
      cv.cd(1);
      leg.Draw("SAME");
      cv.Write();          
      cv.Print(Form("%s",pdf_file.Data())); 
      cv.Clear();    
    }
  }
  cv.Print(Form("%s]",pdf_file.Data()));
}

void compare_KF_O2(TString KF_name="MassResolutionKFChi2_1.root", TString O2_name="MassResolutionO2.root", TString output_name="vertexer_comparison1.root", TString pdf_file="vertexer_comparison1.pdf", TString folder_name="KF_O2_1")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadRightMargin(0.1); 
  gStyle->SetPadLeftMargin(0.15);
  
  gSystem->Exec(Form("mkdir %s",folder_name.Data()));

  const char  lAM[2] = {'A','M'};
  const char* lTitle[2] = {"O^{2} vertexer","Kalman Filter #chi^{2}/NDF<1"};
  //const char* lTitle[3] = {"No cuts","#chi^{2} <0.5","#chi^{2}<5"};
  
  const char* lMeasure[2] = {"Average","StdDev"};
  const char* lRes[8] = {"Mass_pT","Mass_ct","Pt_pT","P_p","Ct_ct","X_X","Y_Y","Z_Z"};
  const char* lDir[2] = {"mean_value","stddev"};

  TFile* input_file[2];
  input_file[0] = new TFile(O2_name.Data());
  input_file[1] = new TFile(KF_name.Data());
  TFile output_file(Form("%s/%s",folder_name.Data(),output_name.Data()),"RECREATE");
  TH1D* fHist[2][2];//[mean/stddev][vertexer]

  TCanvas cv("","",800,450);
  cv.Print(Form("%s/%s[",folder_name.Data(),pdf_file.Data()));
  //antimatter and matter
  for(int iMat=0; iMat<2; iMat++){
    //all the resolutions 
    for(int iRes=0; iRes<8; iRes++){
      //loop on the measurements: mean and stddev
      for(int iMeas=0; iMeas<2; iMeas++){          
        //Kalman filter, O2 
        for(int iFile=0; iFile<2; iFile++){
          std::cout<<Form("%s/f%s_%s_%c",lDir[0],lMeasure[0],lRes[1],lAM[0])<<"\n";
          fHist[iMeas][iFile] = (TH1D*) input_file[iFile]->Get(Form("%s/f%s_%s_%c",lDir[iMeas],lMeasure[iMeas],lRes[iRes],lAM[iMat]));
          set_style(fHist[iMeas][iFile],iFile);
          fHist[iMeas][iFile]->SetTitle(lTitle[iFile]);
        }
        set_range(fHist[iMeas],5);
        rainbow_plot(fHist[iMeas],cv,Form("RainbowPlot_%s_%s_%c",lMeasure[iMeas],lRes[iRes],lAM[iMat]),true);
        output_file.cd();
        cv.Write();          
        cv.Print(Form("%s/%s",folder_name.Data(),pdf_file.Data())); 
        cv.Clear();
      }
      multirainbow_plot(fHist,cv,Form("MultiRainbowPlot_%s_%c",lRes[iRes],lAM[iMat]),true);
      cv.Write();          
      cv.Print(Form("%s/%s",folder_name.Data(),pdf_file.Data())); 
      cv.Clear();    
    }
  }
  cv.Print(Form("%s/%s]",folder_name.Data(),pdf_file.Data()));
}


void compare_distribution(TString Std_name="MassResolutionStd.root", TString KF_name="MassResolutionKF.root", TString O2_name="MassResolutionO2.root", TString output_name="distr_comparison.root", TString pdf_file="distr_comparison.pdf", TString folder_name="nome")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadRightMargin(0.1); 
  gStyle->SetPadLeftMargin(0.15);
  
  //gSystem->Exec(Form("mkdir %s",folder_name));

  const char  lAM[2] = {'A','M'};
  const char* lTitle[3] = {"Standard vertexer","O^{2} vertexer","Kalman Filter"};
  const char* lRes[8] = {"fMass_pT_Int","fPt_pT_Int","fP_p_Int","fCt_ct_Int","fX_X_Int","fY_Y_Int","fZ_Z_Int"};

  TFile* input_file[3];
  input_file[0] = new TFile(Std_name.Data());
  input_file[1] = new TFile(O2_name.Data());
  input_file[2] = new TFile(KF_name.Data());
  TFile output_file(output_name.Data(),"RECREATE");
  TH1D* fHist[3][3];//[mean/stddev][vertexer]

  TCanvas cv("","",800,450);
  cv.Print(Form("%s[",pdf_file.Data()));
  TPaveLabel *titlea = new TPaveLabel(.11,.95,.35,.99,"Antimatter","brndc");
  TPaveLabel *titlem = new TPaveLabel(.11,.95,.35,.99,"Antimatter","brndc");
  
  TLine* line;
  //antimatter and matter
  for(int iMat=0; iMat<2; iMat++){
    //all the resolutions (mass,pt,ct)vs(pt,ct)
    for(int iDist=0; iDist<3; iDist++){
      //Kalman filter, O2 ,standard vertexer
      for(int iFile=0; iFile<3; iFile++){
        fHist[iDist][iFile] = (TH1D*) input_file[iFile]->Get(Form("int_distr/%s_%c",lRes[iDist],lAM[iMat]));
        fHist[iDist][iFile]->Scale(1./fHist[iDist][iFile]->GetMaximum());
        set_style(fHist[iDist][iFile],iFile);
        fHist[iDist][iFile]->GetYaxis()->SetTitle("#frac{counts}{max counts}");
        fHist[iDist][iFile]->SetTitle(lTitle[iFile]);
      }
      set_range(fHist[iDist]);
      
      rainbow_plot(fHist[iDist],cv,Form("MultiRainbowPlot_%c",lAM[iMat]),true);
      cv.cd();
      line = new TLine(0.,0.,0.,1.0);
      line->SetLineColor(2);
      line->SetLineStyle(10);
      line->Draw("SAME");

      if(iMat==0)
        cv.SetTitle("Antimatter");
      else
        cv.SetTitle("Matter");
      
      cv.Write();          
      cv.Print(Form("%s",pdf_file.Data())); 
      cv.Clear(); 
    }   
  }

  //range ct 10 - 20 cm

  for(int iFile=0; iFile<3; iFile++){
    fHist[0][iFile] = (TH1D*) input_file[iFile]->Get(Form("fSlice_10_20_%c",lAM[0]));
    fHist[0][iFile]->Scale(1./fHist[0][iFile]->GetMaximum());
    set_style(fHist[0][iFile],iFile);
    fHist[0][iFile]->GetYaxis()->SetTitle("#frac{counts}{max counts}");
    fHist[0][iFile]->SetTitle(lTitle[iFile]);
  }
  set_range(fHist[0]);
  cv.Print(Form("%s]",pdf_file.Data()));

  rainbow_plot(fHist[0],cv,Form("MultiRainbowPlotSlice_%c",lAM[0]),true);
  cv.cd();
  line = new TLine(0.,0.,0.,1.0);
  line->SetLineColor(2);
  line->SetLineStyle(10);
  line->Draw("SAME");
  cv.Write();          
  cv.Print(Form("slice10_20.pdf"));
  cv.Clear(); 

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

//function to compare efficiencies
void compare_trees(TString findable_name = "efficiencyKF.root",TString task_name = "KFResults.root",TString pdf_file = "compare_trees.pdf")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TFile findable_file(findable_name.Data());
  TFile task_file(task_name.Data());
  const char lAM[3]{"AM"};
  const char* lName[2]={"Antimatter","Matter"};
  TCanvas cv("","",800*2,450);
  cv.Print(Form("%s[",pdf_file.Data()));
  TH1D* histos[2][2];//[eff (ct-pt)][tree]
  //matter and antimatter
  for(int iMat=0;iMat<2;iMat++){
    histos[0][0] = (TH1D*) task_file.Get(Form("efficiency_chi2/fHistEff_Pt_Chi2_deuprot_50.0_%c",lAM[iMat]));
    histos[0][1] = (TH1D*) findable_file.Get(Form("NsigmaTPC/Efficiency/pT/fHistRec_NsigmaTPC_pT_cut_3.0_%c",lAM[iMat]));
    histos[1][0] = (TH1D*) task_file.Get(Form("efficiency_chi2/fHistEff_Ct_Chi2_deuprot_50.0_%c",lAM[iMat]));
    histos[1][1] = (TH1D*) findable_file.Get(Form("NsigmaTPC/Efficiency/ct/fHistRec_NsigmaTPC_ct_cut_3.0_%c",lAM[iMat]));
    histos[0][0]->SetTitle("Task 3KF");
    histos[0][1]->SetTitle("Findable");
    set_range(histos[0]);
    set_range(histos[1]);
    multirainbow_plot(histos,cv,Form("efficiency_%c",lAM[iMat]),true,"PLC PMC");
    cv.Print(Form("%s",pdf_file.Data()));
    cv.Clear();
  }
  cv.Print(Form("%s]",pdf_file.Data()));
}

void compare_times(TString Std_name="selector_resultsStd.root", TString KF_name1="selector_resultsKFChi2_50.root", TString KF_name2="selector_resultsKFChi2_2.root", TString O2_name="selector_resultsO2.root", TString output_name="time_comparison.root", TString pdf_file="time_comparison.pdf", TString folder_name="nome")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  gSystem->Exec(Form("mkdir %s",folder_name.Data()));

  const char* lTitle[4] = {"Standard vertexer","Kalman Filter #chi^{2}/NDF<50","Kalman Filter #chi^{2}/NDF<2","O^{2} vertexer"};

  TFile* input_file[4];
  input_file[0] = new TFile(Std_name.Data());
  input_file[1] = new TFile(KF_name1.Data());
  input_file[2] = new TFile(KF_name2.Data());
  input_file[3] = new TFile(O2_name.Data());

  TLine* line = new TLine(0.,1.,4.,1.);
  line->SetLineColor(2);
  line->SetLineStyle(10);

  TH1D* histTime = nullptr;
  TH1D* histComparison = new TH1D("histComparison",";;#frac{process time}{std vtx time}",4,0,4);

  TCanvas cv("","",800,450);
  for(int iVert=0; iVert<4; iVert++){
    histTime = (TH1D*) input_file[iVert]->Get("fTotTime");
    histComparison->SetBinContent(iVert+1,histTime->GetBinContent(1));
    histComparison->SetBinError(iVert+1,0);
    histComparison->GetXaxis()->SetBinLabel(iVert+1,lTitle[iVert]);
  }

  TGaxis*axis = new TGaxis(4,0.8,4,1.1,0.8*histComparison->GetBinContent(1),1.1*histComparison->GetBinContent(1),510,"+L");
  axis->SetTitle("process time (s)");
  axis->SetLabelFont(0);
  axis->SetTitleFont(0);
  
  histComparison->SetMarkerStyle(8);
  histComparison->SetMarkerColor(4);
  //nomalized to the std vertexer
  histComparison->Scale(1./histComparison->GetBinContent(1));
  histComparison->GetYaxis()->SetRangeUser(0.8,1.1);
  histComparison->GetXaxis()->SetLabelSize(0.05);
  histComparison->Draw("P");

   //draw an axis on the right side
  axis->Draw("same");
  line->Draw("same");
  cv.Print(Form("%s",pdf_file.Data()));
}