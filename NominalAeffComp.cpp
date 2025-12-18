#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include <TH2.h>
#include <TH3D.h>
#include "TCanvas.h"
#include "TGraph.h"
#include "TF1.h"
#include "TKey.h"

#include <iostream>
#include <cstdio>
#include <cstring>
#include <vector>
#include <iomanip>
#include <fstream>
#include <TLine.h>
#include <stdlib.h>
#include "TLegend.h"
#include <TStyle.h>
#include "../EdgesandBins.h"
#include <TChain.h>
#include <TMath.h>
#include <TSystem.h>
#include <TPaveLabel.h>
#include <TF2.h>
#include <TRandom3.h>
#include <TGraphErrors.h>
#include "MStarCamTrans.h"
#include "MGeomCamMagicTwo.h"
#include "MHCamera.h"

int NominalAeffComp(const char *period){

    Double_t PI = TMath::Pi();
    TH1::SetDefaultSumw2(1);
    gStyle->SetOptFit(1);

    Double_t SearchRadius = 0.15;

    // I open the file with the Montecarlo maps of "pre-expected" number of events
    char *Filename = Form("/data/magic/users-ifae/edosouto/MAGIC/BigSearch/Melibea/%s/AcceptanceMaps_%s.root", period, period);
    if(gSystem->AccessPathName(Form("%s",Filename))){
        cout << "ERROR: The file doesn't exist\n";
        return 1;
    } else{
        cout << "The file exists\n" << endl;
        // I will save the results of this macro in the MAPS file, so that I don't have a thousand small files
        TFile *maps = new TFile(Form("%s",Filename), "UPDATE");

        if (!maps || maps->IsZombie()) {
            cout<<"Cannot open file or file is zombie\n";
        } else {
            cout<<"IsOpen: "<<maps->IsOpen()<<"\n";
            cout<<"IsZombie: "<<maps->IsZombie()<<"\n";
            cout<<"IsWritable: "<<maps->IsWritable()<<"\n";   // very important
            maps->ls();
        }

        // Hard-coded only two Zenith bins
        const int nbinsZd = 2;
        TVectorD v(2);  // Vector to store the epsilon ref of the two Zd bins

        TH2D *hmaps[nbinsZd];
        TArrayD edgesZd = SetEdgesCos(nbinsZd, 4.546, 50.1);
        
        for (int i = 0; i < nbinsZd; i++){
            hmaps[i] = (TH2D*)maps->Get(Form("Map%d", i+1));
            hmaps[i]->SetDirectory(0);
            hmaps[i]->SetTitle(Form("MC map for Zd c(%.0f, %.0f)", edgesZd[i], edgesZd[i+1]));
        }
        
        Int_t xi = 0, yi = 0, zi= 0;       // Containers for the bin position
        Double_t ValMC = 0.;
        
        // I calculate epsilon ref for a couple of cells that are at 0.4deg
        // from the center and then use their mean. The number of NN is kinda random.
        // Not too small, not too big so that this doesn't take that much time. 
        const int NN = 11;
        Double_t y[NN];
        Double_t x[NN];
        
        for (int i = 0; i < NN; i++){
            x[i] = -0.4+i*0.8/(NN-1);
            std::cout <<x[i] << std::endl;
            y[i] = -0.4+i*0.8/(NN-1);
        }
        
        Int_t counter = 0;
        Double_t xp = 0., yp = 0.;
        for (int k = 0; k < nbinsZd; k++){
            Int_t histonum = k;
            for (int i = 0; i < NN; i++){
                 xp = x[i];     
                 for (int j = 0; j < NN; j++){
                     yp = y[j];
         
                     //   xp*xp +yp*yp = 0.4 ? then calculate ValMC
                     if ( xp*xp + yp*yp == 0.4*0.4 ){
                         counter +=1;
                         // Section where we obtain the bins around the given point and their contents
                         Int_t centb = hmaps[histonum]->FindBin(xp, yp);
                         hmaps[histonum]->GetBinXYZ(centb, xi, yi, zi);
                         Double_t bwidth = hmaps[histonum]->GetXaxis()->GetBinWidth(3);
                         Int_t binrel = Int_t(SearchRadius/bwidth)+1;
                         for (int ii = xi - binrel - 1; ii <= xi + binrel + 1; ii++){
                             for (int jj = yi - binrel - 1; jj <= yi + binrel + 1; jj++){
                                 
                                 Double_t xcent = hmaps[histonum]->GetXaxis()->GetBinCenter(ii);
                                 Double_t ycent = hmaps[histonum]->GetYaxis()->GetBinCenter(jj);
                                 if ( (xp-xcent)*(xp-xcent)+(yp-ycent)*(yp-ycent) <= SearchRadius*SearchRadius ) {
                                     ValMC += hmaps[histonum]->GetBinContent(ii,jj);
                                 }
                             }
                         }
                     }
                 }
            }
            std::cout << hmaps[histonum]->GetTitle() << std::endl;
            std::cout << counter << std::endl;
            std::cout << " Reference epsilon value for a 0.4deg cell: " << ValMC/Double_t(counter) << std::endl;
            v[histonum] = ValMC/Double_t(counter);
            counter = 0;
            ValMC = 0.;
        }
    maps->cd();
    v.Write("EpsilonRef", TObject::kOverwrite);
    //maps->GetDirectory("")->WriteTObject(&v, "EpsilonRef", "Overwrite");
    maps->Close();
    }

    return 0;
}
