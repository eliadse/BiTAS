#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include <TH2.h>
#include <TH3D.h>
#include "TCanvas.h"
#include "TGraph.h"
#include "TF1.h"
#include "TKey.h"
#include "TNtupleD.h"
#include "TVectorD.h"

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
#include <TChain.h>
#include <TMath.h>
#include <TSystem.h>
#include <TPaveLabel.h>
#include <TF2.h>
#include <TRandom3.h>
#include <TGraphErrors.h>
#include "MObservatory.h"
#include "MStarCamTrans.h"
#include "MGeomCamMagicTwo.h"
#include "MHCamera.h"
#include "MGeomCam.h"
#include "MMath.h"
#include "MArgs.h"
#include "./EdgesandBins.h"

// This is the analysis macro that adds the normalization factor from rejecting events
// from areas in the FoV that correspond to TeVCat sources.
// The known acceptance from that part of the camera is loaded from the lightcurve file
// N0 = N0'/(1-Int_excl(f(r)dr)) That integral is the sum of the acceptance of the excluded regions
// 20/09/2023
// In this version of analysis, we save obs, exp, Zd, mjd, cellN and hoursobs
// Now we save the value from the acceptance map (depends on the position of the cell in the acceptance map) independently bc we need it to calculate the flux
int AnalysisLZdNFE(const char *filename, const double dispCut, const int binwcoarse, const char *period){
    Double_t PI = TMath::Pi();
    TH1::SetDefaultSumw2(1);
    gStyle->SetOptFit(1);
    
    Int_t binEpsilon = 0;
    Double_t binCenter = 0.;

    if(gSystem->AccessPathName(Form("%s",filename))){
        cout << "ERROR: The file doesn't exist\n";
        return -1;
    } else cout << "The file exists\n";
    
    TString fn = filename; 
    TString od = fn(0,fn.Index("/Light"));
    const char *OutDir = od;
    Double_t SearchRadius = 0.15;
    int const binwfine = 10;    // Finer binning to calculate the number of expected events
    std::cout << "This version of Analysis divides the run in bins of equal length apart from the closest to the end of the run" << std::endl;
    std::cout << OutDir << std::endl;

    const MGeomCamMagicTwo geom;
    const MObservatory observ;
    MStarCamTrans *starcamtrans    = new MStarCamTrans(geom, observ);
    Double_t mm2deg = geom.GetConvMm2Deg();

    // First I open the file containing the cells, the events from the 
    // full camera light curve and the times in mjd that there's a 
    // change in wobble
    TFile *f = new TFile(filename, "READ");
  
    // I extract the name of the source and the night of the observation
    TString fname = filename;                                               // We extract info from the file's name
    TString nightStr = fname(fname.Index("ves")+4,10);                        // All subruns after the stereo upgrade start by 050
    TString sourceStr = fname(fname.Index("Curves")+18, fname.Index(".root")-(fname.Index("Curves")+22));
    const char *night = nightStr;
    const char *source = sourceStr;
    std::cout << "Night:  " << night << "\n";
    std::cout << "Source:  " << source << "\n";
    
    //TString dispStr = fname(fname.Index(".root")-3,3);
    //Double_t dispCut = dispStr.Atof();
    //std::cout << "dispCut = " << dispCut << " and search radius = " << SearchRadius << std::endl;
    // Normalization factor to take into account the exclusion of data from the full camera light curve
    // We exclude the events from the region surrounding the source, but those events are bkg or bkg+sig
    // So in any case we're excluding bkg events from the rates used in the montecarlo to obtain the 
    // expected events 
    // In this final form the normalization factor is added directly to the camera lightcurve as a weight.
    // This weight comes from the acceptance of the camera and was calculated in FullCameraLCCoords.cpp
    // Double_t nf =  1.;

    // Load the vector stored in f and get t0 in mjd
    std::vector<Double_t> *wobtimes;
    f->GetObject("WobbleTimes", wobtimes);
    Double_t t0 = wobtimes->at(0);
    Double_t tf = wobtimes->at(wobtimes->size()-1);
    // Get number of hours observed for displaying later        
    //Double_t hh = (tf-t0)*24.;  
    Double_t hh = 0.;   // hours observed, filled later  
    Int_t preNbins = ((tf-t0)*24.*60.*60./binwfine) + 15;
    
    // WobbleTimes is in mjd. Due to memory/precision it doesn't
    // work with the function "WobbleBins", so we create another
    // vector to store it in seconds after t0
    std::vector<Double_t> secwobs;
    //for (int i = 0; i < wobtimes->size(); i++)  secwobs.push_back((wobtimes->at(i)-t0)*60.*24.*60.);
    for (int i = 0; i < wobtimes->size(); i++){
        std::cout << (wobtimes->at(i)-t0)*60.*24.*60. << " ";
          secwobs.push_back((wobtimes->at(i)-t0)*60.*24.*60.);
    }

    Int_t binsc = 0;
    std::cout << "Number of runs: " << secwobs.size()/2 << std::endl;
    //Double_t tl2 = (secwobs.at(1) - secwobs.at(0))/2;   // Half the size of the first run
    //std::cout << "First run lasts " << tl2*2. << "s" << std::endl;
    
    // To get the livetime of a certain night, we sum the wobble times
    Double_t tl2 = 0.;   // Half the size of the first run
    Int_t tout = 0;
    for (int i = 0; i < secwobs.size(); i+=2){

        if ((secwobs.at(i+1) - secwobs.at(i))/2 > tl2) tl2 = (secwobs.at(i+1) - secwobs.at(i))/2;
        hh += secwobs.at(i+1) - secwobs.at(i);
        std::cout << "Run " << (i + 2)/2 << " lasts " << secwobs.at(i+1) - secwobs.at(i) << "s" << std::endl;
        if ( tl2 > binwcoarse/2) tout = 1;
    }
    if (tout == 0) {
        std::cout << "All runs are too short. Exiting..." << std::endl;
        return -1;
    }
    
    // Convert hh from seconds to hours 
    hh = hh/3600.;
//    if (secwobs.size() < 3 && tl2 < binwcoarse/2) {
//        std::cout << "Observation is too small " << std::endl;
//        return -1;
//    }
//    Int_t normal = 890;     // A run should not last less than 15min/900s +- 10s.
//    if (tl2 < normal/2) {
//        tl2 = (secwobs.at(3) - secwobs.at(2))/2;   // If the first run is busted we use the next one
//        std::cout << "First run is useless, the second one lasts " << tl2*2. << "s" << std::endl;
//    }
    if (binwcoarse > tl2) binsc = (wobtimes->size())*1.5 + 10;   // The time window is bigger than half the run, so the number of bins is fixed.  
    else if (binwcoarse <= tl2 && binwcoarse > 20) binsc = ((tf-t0)*24.*60.*60./binwcoarse)*2;     // Window less than half run
    else if (binwcoarse <= 20) binsc = ((tf-t0)*24.*60.*60./binwcoarse) +15;    // Window smaller than the gaps between runs
    //std::cout << "hey 2" << "\n";

    TH1D *HCam = new TH1D("HCam", "Full camera light curve", preNbins, 0.0, preNbins*binwfine);     // In seconds
    //std::cout << "\nHCam wobble bin change \n /////////////////\n";
    WobbleBins(HCam, secwobs);

    // Load the real events for the lightcurve and fill the hist
    TTree *fullcam = (TTree*)f->Get("FullCamLC");
    Double_t tte = 0., w = 0.;
    fullcam->SetBranchAddress("Time", &tte);
    fullcam->SetBranchAddress("Acc", &w);
    for (int i = 0; i < fullcam->GetEntries(); i++){
        fullcam->GetEntry(i);
        if (w > 0.95) continue;
        Double_t timeevent = (tte - t0)*24.*60.*60.;
        HCam->Fill(timeevent, 1./(1.-w));
    }
    Int_t hbins = HCam->GetNbinsX();
    if (HCam->GetBinLowEdge(hbins+1) < binwcoarse) {
        std::cout << "Observations last " << HCam->GetBinLowEdge(hbins+1) << "s. They are shorter than the time window\n"; 
        return -1;  
    }
 
    // I open the file with the Montecarlo maps of "pre-expected" number of events
    TFile *maps = new TFile(Form("/data/magic/users-ifae/edosouto/MAGIC/BigSearch/Melibea/%s/AcceptanceMaps_%s.root", period, period), "READ");
    // Now we also get the epsilon_ref values for these maps from the acceptance maps file
    // It's a vector. First element is lowZd, second element is mediumZd 
    TVectorD *vec = (TVectorD*)maps->Get("EpsilonRef");
    if (!vec) {
        std::cerr << "Could not find the EpsilonRef!" << std::endl;
        return 1;
    }
    // LowZd
    double epsilonref = (*vec)[0];

    const int nbinsZd = 2;
    TH2D *hmaps[nbinsZd];
    TArrayD edgesZd = SetEdgesCos(nbinsZd, 4.546, 50.1);
    
    for (int i = 0; i < nbinsZd; i++){
        hmaps[i] = (TH2D*)maps->Get(Form("Map%d", i+1));
        hmaps[i]->SetDirectory(0);
        hmaps[i]->SetTitle(Form("MC map for Zd c(%.0f, %.0f)", edgesZd[i], edgesZd[i+1]));
    }
    // I close the montecarlo file to free up some memory
    maps->Close();
    std::cout << "Epsilon_ref LowZd = " << epsilonref << std::endl;

    Int_t histonum = 0;
    
    TH1D *hb = new TH1D("hb","hb", preNbins, 0.0, preNbins*binwfine);
    //std::cout << "\nhb wobble bin change \n /////////////////\n";
    WobbleBins(hb, secwobs);
    Int_t nbhb = hb->GetNbinsX();
    
    // Cloning the previous histogram instead of making it again and running the wobble conversion again
    TH1D* hbE = new TH1D(*hb);
    hbE->SetTitle("hbE");
    hbE->SetName("hbE");
//    std::cout << "Check to see if hbE has replaced hb: " << hbE->GetName() << " " << hb->GetName() << std::endl;
//    std::cout << "bin 13 of hb and hbE: " << hbE->GetBinLowEdge(13) << " " << hb->GetBinLowEdge(13) << std::endl;
    // Histogram to store the cell position (just ValMC so that we don't need a 2D histogram for x and y) in the acceptance map every second
    //TH1D *hbE = new TH1D("hbE","hbE", preNbins, 0.0, preNbins*binwfine);
    //std::cout<< "New wobbleBins" << std::endl;
    //WobbleBins(hbE, secwobs);
   
   // for (int i = 0; i <= hbE->GetNbinsX(); i++){
   //     std::cout << hbE->GetBinLowEdge(i+1) << " ";
   // }
   // std::cout << std::endl;
    
    TH1D *hcell = new TH1D("hcell","hcell", binsc, 0.0, binsc*binwcoarse); 
    hcell->GetXaxis()->SetTitle("Time since start of observation [sec]");
    hcell->GetYaxis()->SetTitle("Events");
    //std::cout << "\nhcell wobble bin change \n /////////////////\n";
    WobbleBins(hcell, secwobs);
    
    Int_t NBINS = hcell->GetNbinsX();
    //std::cout << "NBINS = " << NBINS << std::endl;
    const int nb = NBINS+1; 
    //Double_t xbins[nb] = {0.};
    //TODO This is a fucking problem. This needs to be nb, but it needs to be an array
    //Double_t xbins[100] = {0.};
    std::vector <Double_t> xbinsv;

    //for (int i = 0; i <= NBINS; i++){
    //    xbins[i] = hcell->GetBinLowEdge(i+1);  // These are the edges of hcell's bins, so their size is NBINS+1
    //}
    //std::cout << "xbins into the vector: ";
    for (int i = 0; i <= NBINS; i++){
        xbinsv.push_back(hcell->GetBinLowEdge(i+1));  // These are the edges of hcell's bins, so their size is NBINS+1
       // std::cout << xbinsv.at(i) << " ";
    }
    //std::cout << std::endl;
    //std::cout << "xbin size and nb: "<< xbinsv.size() << " " << nb <<std::endl;
    double* xbins = &xbinsv[0];

    const int nbins = NBINS;
    //Container for the Zd position
    //Double_t ZdArr[nbins] = {0.};           // This container has to be like hcell's bins
    Double_t ZdArr[100] = {0.};           // This container has to be like hcell's bins
    int countz = 0;
    //double binedges[3] = {5.,35.,50.};
    //TH2F *htz = new TH2F("htz", "Zenith in time", NBINS, xbins, 2, binedges);
    // Opening the tree with the camera pointing coordinates at intervals of 1s
    TTree *cam = (TTree*)f->Get("CameraCoords");
    Double_t dec = 0., ha = 0., zdc = 0., azc = 0., lst = 0., mjd = 0.; 

    cam->SetBranchAddress("Dec0", &dec);
    cam->SetBranchAddress("HA0", &ha);
    cam->SetBranchAddress("ZD", &zdc);
    cam->SetBranchAddress("AZ", &azc);
    cam->SetBranchAddress("LST", &lst);
    cam->SetBranchAddress("MJD", &mjd);
    
    Int_t camEntries = cam->GetEntries();
    const int nc = nbhb;
    //Double_t CamCo[nc][6] = {{0.},{0.}};
    
    Double_t **CamCo = (Double_t **)malloc(nc * sizeof(Double_t*));
    for (int i = 0; i < nc; i++) {
        CamCo[i] = (Double_t *)malloc(6*sizeof(Double_t));
    }
    
    int countb = 0;
    for (int i = 0; i < camEntries; i++){
        cam->GetEntry(i);
        Double_t tev = (mjd-t0)*24.*60.*60.;
        //std::cout << "countb = " << countb << "/" << nc-1 << " countz = " << countz << "/" << NBINS- 1 << std::endl;
        //With the information from the camera positioning, we create a histogram-like array of multiple dimensions
        if (tev > hb->GetBinLowEdge(countb+1) && tev < hb->GetBinLowEdge(countb+2)){
            CamCo[countb][0] = dec;
            CamCo[countb][1] = ha;
            CamCo[countb][2] = zdc;
            CamCo[countb][3] = lst;
            CamCo[countb][4] = mjd;
            CamCo[countb][5] = azc;
            countb++;
        }
        else if (tev > hb->GetBinLowEdge(countb+2)){
            countb+=2; 
            // There is a gap in the data and without this fix it gets stuck
            // Need to find the next bin that corresponds to this data point
            while (tev > hb->GetBinLowEdge(countb+1) && countb < nc){
                countb++;
            }
        }
        // We need to create a separate histogram-like array for the zenith that has the hcell shape
        if (tev > xbins[countz] && tev < xbins[countz+1] && countz < NBINS){
            ZdArr[countz] = zdc;        // I don't care to do the average, this is way better than the last approximation already
            countz++;
        }
        else if (tev > xbins[countz+1] && countz <NBINS){
            countz+=2;
            // There is a gap in the data and without this fix it gets stuck
            // Need to find the next bin that corresponds to this data point
            while (tev > xbins[countz] && countz < NBINS){
                countz++;
            }
        }
        //if (countz >= NBINS + 1 ) break;
        if ( countb >= nc + 1 ) break; 
    }

    // Then I get the list of trees, with variable names
    TIter nextit(f->GetListOfKeys());
    TKey *key;
    TString treen;
    std::vector <Int_t> index_read;
    // Read the tree names. Now there's also a histogram, 
    //so I have to check that the current object is a tree
    while ((key=(TKey*)nextit())) 
    {
//         TObject *obj = key->ReadObj();
        treen = key->GetName();
        if (treen[0] == 't'){
            TTree *t = (TTree*)f->Get(treen);
            Int_t numents = t->GetEntries();
     //       if (numents < 100) continue;
            TString cosi = treen(2,treen.Sizeof()-1);
            Int_t index = cosi.Atoi();
            index_read.push_back(index);
        }
        else continue;
    }
    Int_t tnum = 0;
    std::cout << "Bin width selected: " << binwcoarse << " and CamCoord interval: " << binwfine << std::endl; 
    // Instead of analysing everything here, we save the p-values and other relevant info to be merged as needed 
    //TNtupleD *Results = new TNtupleD("Results","Results","Obs:Exp:Zd:MJD:CellN:HoursObs");
    //TNtupleD *Results = new TNtupleD("Results","Results","Obs:Exp:Zd:MJD:CellN:ValMap");
    // Instead of saving ValMap, which is epsilon_cell in my thesis, we load epsilon_ref
    // from the Maps file and apply it to directly save epsilon, their division
    TNtupleD *Results = new TNtupleD("Results","Results","Obs:Exp:Zd:MJD:CellN:Epsilon");
    TFile *fnew = new TFile(Form("%s/Results/Results_%s_%s_%is.root", OutDir, night, source, binwcoarse), "RECREATE");

    Int_t xi, yi, zi;       // Containers for the FetBinXYZ function
    Double_t xp = 0., yp = 0.;  // Containers for the center of the bin in the camera
    Double_t xip = 0., yip = 0.;  // Containers for the center of the bin in the camera
    Double_t phi0 = 0.;

   /////////// From the root forum. Looping over trees
   //    TTree *T; 
   // while ((key=(TKey*)next())) { 
   //     if (strcmp(key->GetClassName(),"TTree")) continue; //do not use keys that are not trees 
   //     T = (TTree*)file->Get(key->GetName()); //fetch the Tree header in memory 
   //     totalEvents = T->GetEntries(); 
   //     printf("%50s\t%i\n",keyname,totalEvents); 
   //     delete T;
   // }
   ///////////
    TTree *tt;
    
    for (int j = 0; j < index_read.size(); j++){
//    for (int j = 0; j < 3; j++){
        tnum = index_read.at(j);
        //TTree *tt = (TTree*)f->Get(Form("t_%i", tnum));
        tt = (TTree*)f->Get(Form("t_%i", tnum));
        Int_t NEntries = tt->GetEntries();
        std::cout << "Tree, entries: " << tnum << ", " << NEntries << std::endl;
        
        hcell->SetTitle(Form("Histogram from cell %d, %s FoV", tnum, source));
        
        Double_t timepoint, RAG, DecG, ZdEv, gridx, gridy;
        //Should be "time:RAGrid:DecGrid:IndexGrid:Zd:Az:ExpectRate:GridY:CamX:CamY:Wangle"
        tt->SetBranchAddress("Time",&timepoint);
        tt->SetBranchAddress("RA",&RAG);
        tt->SetBranchAddress("Dec",&DecG);
        tt->SetBranchAddress("Zd", &ZdEv);         // Pointing pos coords
        tt->SetBranchAddress("GridX", &gridx);     // The gridx or gridx are the expected rates variable 
        tt->SetBranchAddress("GridY", &gridy);
        
        for (int k = 0; k < NEntries; k++){
            tt->GetEntry(k);
            Double_t tevent = timepoint;     
            tevent = (tevent - t0)*24.*60.*60.;       // now tevent is seconds after start 
            hcell->Fill(tevent);
            //htz->Fill(tevent, ZdEv);
        }
        
        tt->GetEntry(0);
        Double_t cellRA = RAG;      // This is a fixed value. One per cell
        Double_t cellDec = DecG;
        for (int l = 1; l <= nbhb; l++){
            Double_t fcont = HCam->GetBinContent(l);
            if(fcont == 0)  continue;        // Bad observation time
            
            Double_t ValMC = 0.;
            starcamtrans->Cel0CelToCam(CamCo[l-1][0], CamCo[l-1][1], cellDec, CamCo[l-1][3]-cellRA, xip, yip);
            xip *= mm2deg;
            yip *= mm2deg;

            // Azimuth rotation part
            phi0 = CamCo[l-1][5] - 120.;
            phi0 *= PI/180.; // Now in radians

            xp = xip*cos(phi0) - yip*sin(phi0);
            yp = xip*sin(phi0) + yip*cos(phi0);

            if (CamCo[l-1][2] >= 35) histonum = 1;
            else histonum = 0;
            if (xp*xp+yp*yp > (dispCut-SearchRadius)*(dispCut-SearchRadius)) continue;
            
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
            //fcont*=nf;
            // Fill the histogram to rebin it later
            hb->SetBinContent(l, fcont*ValMC);
            hbE->SetBinContent(l, ValMC);
        }
        TH1D *hexp = dynamic_cast<TH1D*>(hb->Rebin(NBINS, "hexp", xbins));
        //std::cout << "Bin 13 edges and width: " << hexp->GetBinLowEdge(13) << " " << hexp->GetBinLowEdge(14) << " " << hexp->GetBinWidth(13) << std::endl;
        
        Double_t lastvE = 0.;
        for (int m = 1; m <= hcell->GetNbinsX(); m++){

            //std::cout << "Zd = " << ZdArr[m-1] << " bin = " << m << std::endl;
            // Counting only the bins with the selected time width
            if (hcell->GetBinWidth(m) != binwcoarse) continue;
            Double_t val = hcell->GetBinContent(m);
            Double_t expv = hexp->GetBinContent(m);
            if (expv == 0) continue;

            binCenter = hcell->GetBinCenter(m);
            binEpsilon = hbE->FindBin(binCenter);
            Double_t vE = hbE->GetBinContent(binEpsilon);
            // Enter the loop to find a non zero bin
            Int_t count = 0;
            while (vE == 0){
                count++;
                vE = hbE->GetBinContent(binEpsilon+count);
                if (binEpsilon+count > nbhb) vE = lastvE;
            }
            //std::cout << vE << " " << binEpsilon+count << std::endl;
            Double_t mjdbin = (hcell->GetBinCenter(m)/24./3600.) + t0;
            // ZdArr has the standard 0 to n-1 shape, unlike the 1 to n histogram shape
            // I don't think it's relevant to save hh like this, because it relates to the night and source observed, but after this we execute ResultsPerCell, which erases the connection between cell and source.
            //Results->Fill(val, expv, ZdArr[m-1], mjdbin, tnum, vE);       
            Results->Fill(val, expv, 5.01, mjdbin, tnum, vE/epsilonref);       
            lastvE = vE; 
        }
//        std::cout << "deleting histograms" << std::endl;
        hcell->Reset();
        hexp->Delete();
        hb->Reset();
        hbE->Reset();

        //tt->Reset();
        delete tt;
    }
    /*
    // Appending the observed hours to the file
    std::ofstream outputFile("/data/magic/users-ifae/edosouto/MAGIC/BigSearch/Code/MethodExpect/MediumZdOnly/Livetime.txt", std::ios::app);
    if (outputFile.is_open()) {

        std::cout << "Opening livetime file\n";
        std::cout << "Hours observed = " << hh << std::endl;
        // Source Night Hours Number of runs
        outputFile << source << " " << night << " " << hh << " " << secwobs.size()/2 << std::endl;
        outputFile.close();
    } 
    else {
        std::cout << "Unable to open the file." << std::endl;
    }
    */
    //htz->Draw("zcol");
    // SpaceTrials factor = 0.93
    // TimeTrials factor = 0.87 for 50% overlap. But here there is no overlap
    Results->Write("", TObject::kOverwrite);
    fnew->Write();
    
    fnew->Close();
    delete Results;
    return 0;
}
