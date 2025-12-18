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
#include "MMcEvtBasic.h"
#include "MSrcPosCam.h"
#include "MEnergyEst.h"
#include "/data/magic/users-ifae/edosouto/MAGIC/BigSearch/Code/MARS3version/EdgesandBins.h"

#include <iostream>
#include <cstdio>
#include <cstring>
#include <vector>
#include <iomanip>
#include <fstream>

///////////
// Macro to generate the acceptance maps required for the BigSearch pipeline (Big MAGIC Search)
// BMAS
// It works on two large Zd bins, 5to35 and 35to50
// TODO add hadronness cuts from efficiency
// The energy threshold in medium Zd range is higher
// /////////////////

void AcceptanceMapsGen(const char *period){
    
    const Int_t nbinsZd = 2;
    const Double_t dispCut = 1.4;
    //const Double_t dispCut = 2.2;
    
    Double_t PI = TMath::Pi();
    TH1::SetDefaultSumw2(1);
    gStyle->SetOptFit(1011);

    TChain* data1 = new TChain("Events");
    data1->Add(Form("/data/magic/users-ifae/edosouto/MAGIC/BigSearch/%s_diffuse/za05to35/GA_5.0to35.0_Q_join_2.*.root", period));
    data1->Add(Form("/data/magic/users-ifae/edosouto/MAGIC/BigSearch/%s_diffuse/za35to50/GA_35.0to50.0_Q_join_2.*.root", period));
    
    MStereoPar *ParDisp1         = new MStereoPar();
    MHadronness *hadro1          = new MHadronness();
    MMcEvt *McEv                 = new MMcEvt();
    MPointingPos *poPos          = new MPointingPos();
    MEnergyEst *EnEst            = new MEnergyEst();

    data1->SetBranchStatus("*", 0);
    data1->SetBranchStatus("MMcEvt_1.*",1);
    data1->SetBranchStatus("MStereoParDisp.*",1);
    data1->SetBranchStatus("MPointingPos_1.*",1);
    data1->SetBranchStatus("MHadronness.fHadronness",1);
    data1->SetBranchStatus("MEnergyEst.fEnergy",1);
    
    data1->SetBranchAddress("MStereoParDisp.", &ParDisp1);
    data1->SetBranchAddress("MHadronness.", &hadro1);
    data1->SetBranchAddress("MMcEvt_1.", &McEv);
    data1->SetBranchAddress("MPointingPos_1.", &poPos);
    data1->SetBranchAddress("MEnergyEst.", &EnEst);

    TArrayD edgesZd = SetEdgesCos(nbinsZd, 4.546, 50.1);      // These limits make it so that both MC files are used separated
    
    Double_t    SumEs[nbinsZd] = {0.};
    Double_t    K[nbinsZd] = {0.};  //(kdata/kMC)*(E^-adata/E^-aMC)
     
    Double_t  hadCut = 0.3;
    Double_t  had = 0., zde = 0., aze = 0.;                //Containers to check each event's hadronness, Zd angle, Az angle and Energy.
    Double_t  xDisp = 0., yDisp = 0.;
    Double_t  etrue = 1.0, edisp = 1.0;
    Double_t  alphadata = 2.7;                          // Real data has this energy slope (Cosmic rays up to knee)
    Double_t  alphamc = 1.6;                            // The MC events were generated with this slope
    
    for (int i=0; i<data1->GetEntries(); i++){
        data1->GetEntry(i);
        
        had = hadro1->GetHadronness();

        // Cut in hadronnes (to be changed TODO)
        if(had>hadCut) continue;

        /// Generation radius cut
        xDisp = ParDisp1->GetDirectionX();
        yDisp = ParDisp1->GetDirectionY();
        if (TMath::IsNaN(xDisp) == 1 || TMath::IsNaN(yDisp) == 1 || xDisp*xDisp+yDisp*yDisp > dispCut*dispCut) continue;

        // Energy cut 
        etrue = McEv->GetEnergy();
        edisp = EnEst->GetEnergy();
        //std::cout<< etrue << " " << edisp << " " << xDisp << " " << yDisp << std::endl; 
        if (etrue < 50 || etrue > 15000) continue;
        if (edisp < 70 || edisp > 10000) continue;
        
        zde = poPos->GetCorrZd();    // Everything in poPos
        aze = poPos->GetCorrAz();     
        for (int k = 0; k< nbinsZd; k++){
            if(zde>edgesZd[k] && zde <=edgesZd[k+1]){
                
                // Second energy cut for higher zenith
                if(k == 1 && etrue < 150) continue;  
                if(k == 1 && edisp < 200) continue;
                
                // This is the energy spectrum slope normalization part   
                SumEs[k] += pow(etrue, -alphadata)/pow(etrue, -alphamc); 
                break;
            }
        }
    }
    
    for (int k = 0; k < nbinsZd; k++){
        // ratio of the ks
        K[k] = 1./SumEs[k];
        std::cout << "K[" << k << "] = " << K[k] << "SumEs = " << SumEs[k] << std::endl; 
    }

    // we create a TNtuple to save the maps as leaves
    //TNtupleD *ExpMap = new TNtupleD("ExpMap","ExpMap","Zd:Az:EnDisp:xDisp:yDisp:Nexpect");
    TNtupleD *ExpMap = new TNtupleD("ExpMap","ExpMap","Zd:EnDisp:xDisp:yDisp:Nexpect");
    TFile *fnew = new TFile(Form("/data/magic/users-ifae/edosouto/MAGIC/BigSearch/Melibea/%s/AcceptanceMaps_%s.root",period,period), "RECREATE");
    std::cout << "Name of the created file:  " << fnew->GetName() << "\n";
    
    Double_t  et = 0., ed = 0.;                //Containers to check each event's estimated variables
    Double_t  zdt = 0., azt = 0.;                              // Containers with the true Zd angle, Az angle and Energy.    
    Double_t phi0 = 0., xp = 0., yp = 0.;
    
//     const double mapBinWidth = 0.03;
    const double mapBinWidth = 0.005;
    TH2D *hex[nbinsZd];
    Int_t numbins = 2*dispCut/(mapBinWidth);
    
    for (int i = 0; i < nbinsZd; i++){
        hex[i] = new TH2D;
        hex[i]->SetBins(numbins,-dispCut,dispCut,numbins,-dispCut,dispCut);
        hex[i]->SetTitle(Form("MC expected map for Zd c (%.0f, %.0f)", edgesZd[i], edgesZd[i+1]));
    }
    
    // Loop to fill the histograms with the montecarlo data
    for (int i=0; i<data1->GetEntries(); i++){
        data1->GetEntry(i);
       
        // TODO efficiency cuts or something 
        had = hadro1->GetHadronness();
        if(had > hadCut) continue;
//        std::cout << "had: " << had << std::endl;

        et = McEv->GetEnergy();
//std::cout << "Energies: " << et << std::endl;
        ed  = EnEst->GetEnergy();
//std::cout << "Energies2: " << ed << std::endl;
        if (et < 50. || et > 15000.) continue;
        if (ed < 70. || ed > 10000.) continue;
//        std::cout << "Energies: " << et << " " << ed << std::endl; 
        xDisp = ParDisp1->GetDirectionX();
        yDisp = ParDisp1->GetDirectionY();
        if (TMath::IsNaN(xDisp) == 1 || TMath::IsNaN(yDisp) == 1 || xDisp*xDisp+yDisp*yDisp > dispCut*dispCut) continue;
        
        azt = poPos->GetCorrAz();       // You can check it in Mars/mpointing/MPointingPosCalc.cc
        zdt = poPos->GetCorrZd();     // The pointing position in Mc files comes from these variables        
        phi0 = azt - 120.;
        phi0 *= PI/180.; // Now in radians
        
        xp = xDisp*cos(phi0) - yDisp*sin(phi0);
        yp = xDisp*sin(phi0) + yDisp*cos(phi0);
        
        Double_t E = pow(et, -alphadata)/pow(et, -alphamc); 
//        std::cout << "i = " << i << ", zdt: " << zdt << " E est: " << ed << std::endl;
        // Loops of conditions to populate each variable interval
        for (int k = 0; k< nbinsZd; k++){
            // Extra conditions for the medium zd range
            if(k == 1 && et < 150) continue;  
            if(k == 1 && ed < 200) continue;
            if(zdt>edgesZd[k] && zdt <=edgesZd[k+1]){

                //ExpMap->Fill(zdt, azt, ed, xp, yp, K[k]*E);       // The error per bin will be computed as sqrt(sum of squares of weight) for each bin.
                ExpMap->Fill(zdt, ed, xp, yp, K[k]*E);       // The error per bin will be computed as sqrt(sum of squares of weight) for each bin.
                hex[k]->Fill(xp, yp, K[k]*E);
                if (i%10000==0){
                    std::cout << "edgesZd[" << k << "] = "<< edgesZd[k] << ", edgesZd[" << k+1 << "] = " << edgesZd[k+1] << std::endl;
                    std::cout << "zdt: " << zdt << " Entrue: " << ed << std::endl;
                }
                break;
            }
        }
    }
    
    ExpMap->Write("", TObject::kOverwrite);
    fnew->Write();
    fnew->WriteObject(hex[0], "Map1");
    fnew->WriteObject(hex[1], "Map2");
    
    fnew->Close();
    delete ExpMap;
    
    std::cout << "Done!\n";

}



