#include <TF1.h>
#include <TMarker.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TVirtualPad.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TChain.h>
#include <TString.h>
#include <TImage.h>
#include <TText.h>
#include <TEventList.h>
#include <TStreamerInfo.h>
#include <TList.h>
#include <TFile.h>
#include <TFrame.h>
#include <TLatex.h>
#include "TNtupleD.h"

#include "MStatusDisplay.h"
#include "MStatusArray.h"
#include "MLog.h"
#include "MPointingPos.h"
#include "MImageParDisp.h"
#include "MStereoPar.h"

#include "MTime.h"
#include "MGeomCam.h"
#include "MGeomCamMagicTwo.h"
#include "MHCamera.h"
#include "MArgs.h"
#include "MArray.h"
#include "MObservatory.h"
#include "MStarCamTrans.h"
#include "MParContainer.h"
#include "MMath.h"
#include "MEnv.h"
#include "MGraphicsWizard.h"
#include "MRawEvtHeader.h"
#include "MHadronness.h"
#include "MEnergyEst.h"
#include "MMcEvtBasic.h"
#include "MSrcPosCam.h"
#include "/data/magic/users-ifae/edosouto/MAGIC/BigSearch/Code/MethodExpect/EdgesandBins.h"

#include <iostream>
#include <cstdio>
#include <cstring>
#include <vector>
#include <iomanip>
#include <fstream>
#include <TPaveLabel.h>

int RingCollectionArea_LZd(const char *period){
//void CollectionArea(const Int_t nFineTrue, const Double_t minE, const Double_t maxE){

    const Int_t nFineTrue = 100;
    const Double_t minE = 10.;
    const Double_t maxE = 60000.;

    Double_t PI = TMath::Pi();
    TH1::SetDefaultSumw2(1);
    gStyle->SetOptStat(0);
    
    char *direction = Form("/data/magic/users-ifae/edosouto/MAGIC/BigSearch/Ringwobble_MC/%s/za05to35", period);  

    // Check that the file exists
    //if(gSystem->AccessPathName(Form("%s/GA_*.root",direction))){
    if(gSystem->AccessPathName(Form("%s/",direction))){
        cout << "ERROR: The file doesn't exist\n";
        return -1;
    } else cout << "The file exists\n";

    TChain* data1 = new TChain("Events");
    data1->Add(Form("%s/GA_*.root", direction));

    // Original Montecarlo events data. This container has the energy and telescope's phi and theta of the event
    // This is for the Nsimulated of the Aeff
    TChain* Evtc = new TChain("OriginalMC");
    Evtc->Add(Form("%s/GA_*.root", direction));
    
    MHadronness *hadro1          = new MHadronness();
    MMcEvt *McEv                 = new MMcEvt();
    MMcEvtBasic *Original        = new MMcEvtBasic();
    MStereoPar *Disp             = new MStereoPar();
    MPointingPos *poPos          = new MPointingPos();
    MEnergyEst *EnEst            = new MEnergyEst();


    data1->SetBranchStatus("*", 0);
    Evtc->SetBranchStatus("*", 0);
    data1->SetBranchStatus("MHadronness.fHadronness",1);
    data1->SetBranchStatus("MStereoParDisp.*",1);
    data1->SetBranchStatus("MMcEvt_1.*",1);
    data1->SetBranchStatus("MPointingPos_1.*",1);
    data1->SetBranchStatus("MEnergyEst.*",1);

    Evtc->SetBranchStatus("MMcEvtBasic_1.*", 1);
    data1->SetBranchAddress("MHadronness.", &hadro1);
    data1->SetBranchAddress("MPointingPos_1.", &poPos);
    data1->SetBranchAddress("MMcEvt_1.", &McEv);
    data1->SetBranchAddress("MStereoParDisp.", &Disp);
    Evtc->SetBranchAddress("MMcEvtBasic_1.", &Original);
    data1->SetBranchAddress("MEnergyEst.", &EnEst);
   
    TArrayD edgesEnt = SetEdgesLog(nFineTrue, minE, maxE);
    
    Double_t    Nsim[nFineTrue];
    Double_t    Nfinal[nFineTrue];
    Double_t    Aeff[nFineTrue];
    Double_t    maxr = 0.;
    
    for (int mm=0; mm<nFineTrue; mm++){
        Nsim[mm] = 0.0;
        Nfinal[mm] = 0.0;
        Aeff[mm] = 0.0;
    }

    Double_t  hadCut  = 0.3;
    Double_t  had = 0., zde = 0., aze = 0., ee = 0., et = 0.;                //Containers to check each event's hadronness, Zd angle, Az angle and Energy.
    Double_t  zdd = 0.0, azd = 0.0, enor = 0.0;
    Double_t  etrue = 1.0, eest = 1.0;
    
    std::string period_flat = period;
    std::string period_str = period;
    period_flat.erase(std::remove(period_flat.begin(), period_flat.end(), '.'), period_flat.end());
    
    // Opening the merged results file of the period to add the effective area from the ringwobble MC
    //TFile *f = new TFile(Form("/data/magic/users-ifae/edosouto/MAGIC/BigSearch/Code/MARS3version/%s/merged_%s_100s.root", period, period_flat), "UPDATE");
    std::string filename = "/data/magic/users-ifae/edosouto/MAGIC/BigSearch/Code/MARS3version/" + period_str + "/merged_"+ period_flat + "_100s.root";

    TFile *f = new TFile(filename.c_str(), "UPDATE");
    
    TH1D* hAeff0 = new TH1D("hAeff0","hAeff0", nFineTrue, log10(minE),log10(maxE));
       
    // The original simulation info 
    Double_t impactr = 0.;
    for (int ii = 0; ii<Evtc->GetEntries(); ii++){
         Evtc->GetEntry(ii);
         
         enor = Original->GetEnergy();
         impactr = Original->GetImpact();
         //if (enor > maxE || enor < minE) continue;      // Pongo corte en los eventos originales por coherencia con los mapas de aceptancia

                if (impactr > maxr) maxr = impactr;         // Get the maximum impact radius simulated 
                for (int mm=0; mm<nFineTrue; mm++){
                    if (enor > edgesEnt[mm] && enor <= edgesEnt[mm+1]){
                        
                        Nsim[mm] +=1.;
                         if (ii%100000 == 0) {
                         std::cout << edgesEnt[mm] << "\t" << edgesEnt[mm+1] << std::endl;
                     }
                        break;
                    }
                }
    }
    std::cout << "Maximum impact radius: " << maxr << std::endl;
    Double_t xcam = 0., ycam = 0.;

    Double_t dispCut = 1.4;
    //Double_t dispCut = 0.15;
    std::cout << data1->GetEntries() << std::endl;
    
    for (int ii = 0; ii < data1->GetEntries(); ii++) {
        data1->GetEntry(ii);
        
        had = hadro1->GetHadronness();
        if(had>hadCut) continue;        // Maybe no hadronness cut for the paper test
        
        etrue = McEv->GetEnergy();
        //if (etrue > maxE || etrue < minE) continue;  // I don't think the cuts on energy are necessary
        eest = EnEst->GetEnergy();
        if (eest > 10000 || eest < 70) continue;
        /// Generation radius cut
        Double_t xDisp = Disp->GetDirectionX();
        Double_t yDisp = Disp->GetDirectionY();
        if (xDisp*xDisp+yDisp*yDisp > dispCut*dispCut) continue;     // No dispcut for the paper test
        
        zdd = poPos->GetCorrZd();    // Zd and Az angle of the original generated MC events, now in degrees
        
                for (int mm=0; mm<nFineTrue; mm++){
                    if (etrue > edgesEnt[mm] && etrue <= edgesEnt[mm+1]){
                        
                        Nfinal[mm] +=1.;
                        break;
                    }
                }
    }

        // Another way of saving the maps
    ofstream outputFile2;
    outputFile2.open(Form("/data/magic/users-ifae/edosouto/MAGIC/BigSearch/Melibea/%s/05to35/AEff_%.0fto%.0fGeVin%ibinsza05to35EestCut.txt", period, minE, maxE, nFineTrue));
    outputFile2 << "# TrueEnergy"; 

    outputFile2 << " Aeff[05deg to 35deg Zd]";
    outputFile2 << "\n";

//    for (int ii = 0; ii<nCoarseZd; ii++){
//        maxr[ii] = maxr[ii]/100.;
//        outputFile << "\n["; 
//        for (int mm = 0; mm<nFineTrue; mm++){
//            if (Nsim[ii][mm] == 0.0) Aeff[ii][mm] = 0.0;
//            else Aeff[ii][mm] = (Nfinal[ii][mm]/Nsim[ii][mm])*PI*maxr[ii]*maxr[ii];
//            hAeff0[ii]->SetBinContent(mm+1, Aeff[ii][mm]);
////             std::cout << edgesZd[ii+1] << ", " << edgesAz[jj+1] << ", " << edgesEnt[mm+1] << ", " << Nfinal[ii][jj][mm]  << ", " << Nsim[ii][jj][mm] << ", " << Aeff[ii][jj][mm] << std::endl;
//            outputFile << "[" << edgesEnt[mm+1] << "," << Aeff[ii][mm] << "],\n";
////                 std::cout << edgesZd[ii+1] << "," << edgesAz[jj+1] << "," << edgesEnt[mm+1] << "," << h1[ii][jj]->GetBinContent(mm+1) << "\n";
//            std::cout << edgesZd[ii+1] << "," << edgesEnt[mm+1] << "," << Aeff[ii][mm] << "\n";
//        }
//        
//    }    
    
    for (int mm = 0; mm<nFineTrue; mm++){
        // All zenith ranges have the same energy bins, so they go in the first column.
        outputFile2 << edgesEnt[mm+1] << " "; 
            if (Nsim[mm] == 0.0) Aeff[mm] = 0.0;
            else Aeff[mm] = (Nfinal[mm]/Nsim[mm])*PI*maxr*maxr/10000.;
            
            // We fill the histogram for visualization
            hAeff0->SetBinContent(mm+1, Aeff[mm]);
            std::cout << edgesEnt[mm+1] << ", " << Nfinal[mm]  << ", " << Nsim[mm] << ", " << Aeff[mm] << std::endl;
            // We save the area in the txt
            outputFile2 << Aeff[mm] << "\n";
    }
       
//    outputFile.close();
    outputFile2.close();
/*
    TF1 *fc = new TF1("fc","[0]*TMath::Power(x,[1])", 50., 15000);
    fc->SetParameter(0, 2.0e8);
    fc->FixParameter(1, -2.7);
    fc->SetParLimits(0, 1.0e+06, 1.0e+10);    
    fc->SetParLimits(1, -2.8, -2.6);
    fc->SetLineColor(kRed);
    hw[0][0]->Fit(fc,"R");
  */  
        TCanvas *cmaster = new TCanvas("cmaster", "cmaster", 900, 900);
        cmaster->SetLogy(); 

            hAeff0->SetTitle("05deg to 35deg");
            hAeff0->GetXaxis()->SetTitle("log(E/GeV)");
            hAeff0->GetYaxis()->SetTitle("\\mbox{Collection area } [m^{2}]"); 

            hAeff0->SetLineColor(kPink);
            hAeff0->Draw();
       cmaster->SaveAs(Form("/data/magic/users-ifae/edosouto/MAGIC/BigSearch/Melibea/%s/05to35/AEff_%.0fto%.0fGeVin%ibinsza05to35EestCut.png", period, minE, maxE, nFineTrue)); 
       
       //new TFile("CollArr05to35.root", "RECREATE");
       // Saving the effective area as a histogram in the merged results file
       hAeff0->Write();
       f->Close();
       return 0;
}
