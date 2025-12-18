#include <TChain.h>
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
#include <TGraphErrors.h>
#include <TSpline.h>

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

#include <iostream>
#include <cstdio>
#include <cstring>
#include <vector>
#include <iomanip>
#include <fstream>
#include "./EdgesandBins.h"
#include <TFitResult.h>
#include "SearchMC.h"

/////////////////////////////////////////////////////////////////
//
//  This is the macro for loading runs of observation and creating a full 
//  light curve of the camera for the whole night. It also obtainss the 
//  timestamps for the changes of wobble.
//  Both elements are then saved to the corresponding "Lightcurves" file 
//  as a tntuple and a vector, for further processing.
//
//  This step is needed to obtain the number of expected events from the
//  acceptance maps. Run after creating Lightcurves_night_source.root and
//  before running Analysis.cpp
//
/////////////////////////////////////////////////////////////////
Double_t localsiderealT(Double_t GMST)  //GMST is in radians
{
    const MObservatory observ;        // No * for them
    Double_t Long = 0.;
    Double_t LocSidT = 0.;

    Long = observ.GetLongitudeDeg();
    Long *= 1./15.;                     // Change deg to hour
    GMST *=24./TMath::TwoPi();          // Radians to hour

    LocSidT = GMST+Long;
    if (LocSidT<0.) LocSidT+=24.;
    return LocSidT;
}
int FullCameraLCCoordsTeVLZd(const char *runfiles, const char *LCfile, const char *source, const double dispCut, const char *period){
    Double_t PI = TMath::Pi();
    TH1::SetDefaultSumw2(1);
    gStyle->SetOptFit(1);

    TString LCname = LCfile;
//    TString dispStr = LCname(LCname.Index(".root")-3,3);
//    Double_t dispCut = dispStr.Atof();
//     Double_t dispCut = 1.4; 
    Double_t hadCut = 0.3;
    
    TString sourceStr = source;
                                                            // From odie.rc the LE signal cut for theta2 plot is sqrt(0.02)
    const MGeomCamMagicTwo geom;                            // For mm2deg
    const MObservatory observ;        // No * for them


    Double_t mm2deg = geom.GetConvMm2Deg();                 // Conversion factor between mm and degrees in the camera
    
    TChain *data = new TChain("Events");
    // Example of runfiles entry ../Melibea/ST.03.03/B22319+31/2013_08_03/20*
    // DependenceZdCheck needs the night information
    data->Add(runfiles);
   
    Double_t t0 = 0., tf = 0.;
    Double_t xDisp = 0., yDisp = 0.;
    Double_t tevent = 0.;
    Double_t RaTel0 = 0., HaTel0 = 0., DecTel0 = 0.;         // For the original pointing of the telescope
    Double_t LST = 0.;                                       // Local sidereal time
    Double_t HA0_i = 0., Dec0_i = 0.;   

    MTime *mtime                   = new MTime();            // I create the objects to access the methods
    MHadronness *hadro             = new MHadronness();      // and all the Getters to get the info from the tree
    MStereoPar *stereoDisp         = new MStereoPar();
    MPointingPos *poPos            = new MPointingPos();
    MSrcPosCam *SrcPos             = new MSrcPosCam();
    MStarCamTrans *starcamtrans    = new MStarCamTrans(geom, observ);
    
    data->SetBranchStatus("*", 0);
    data->SetBranchStatus("MTime_1.*", 1);
    data->SetBranchStatus("MStereoParDisp.*", 1);
    data->SetBranchStatus("MPointingPos_1.*", 1);
    data->SetBranchStatus("MHadronness.fHadronness", 1);
    data->SetBranchStatus("MSrcPosCam_1.*", 1);        
    
    data->SetBranchAddress("MTime_1.", &mtime);
    data->SetBranchAddress("MStereoParDisp.", &stereoDisp);
    data->SetBranchAddress("MPointingPos_1.", &poPos);
    data->SetBranchAddress("MHadronness.", &hadro);
    data->SetBranchAddress("MSrcPosCam_1.", &SrcPos);
  
    Int_t nEntries = data->GetEntries();
    data->GetEntry(0);
    t0 = mtime->GetMjd();
    double dummyt = 0.;
    std::cout << "first event corresponds to: " << data->GetTreeNumber() << std::endl;
    
    // Original pointing positions to obtain the Src in the FoV
    // The thing is that, since I'm loading all the runs from an observation, It's a bit useless to do this, especially for wobble observations.
    // This part works only as a starting point to load coordinates.
    DecTel0 = poPos->GetCorrDec();                    // Pointing of the camera center in Ha/Dec (already corrected
    HaTel0 = poPos->GetCorrHa();
    
    LST = localsiderealT(mtime->GetGmst());    
    RaTel0 = LST-HaTel0;
    
    Double_t RaTel0_deg = RaTel0*360.0/24.0;        // Now we have the R.A. in degrees, like in the grid file
    if (HaTel0 == 0) {
        RaTel0_deg = poPos->GetRa()*360.0/24.0;
    }
    
    data->GetEntry(nEntries-1);
    tf = mtime->GetMjd();
    
    Double_t ExclReg = 0.2;    // In Crab Odie plots, the excess signal goes up to 0.04theta^2 so just to be sure we use this cut to be more conservative, instead of the one from odie.rc
    
    //     if (sourceStr == "HESSJ1841-055")   ExclReg = 0.41;   // 280.229166666667,-5.55
    //     if (sourceStr == "G24.7+06")        ExclReg = 0.21;   // This SNR was not detected. Ignore the different radius
    
    // Since multiple sources can be in the same FoV, we can't select the cut just by the name of the file. We'll select it based on the coordinates from TevCat
    Double_t TevReg = 0.55;                                // Largest wobble offset from observations other than Crab, to look for TeV sources
    if (sourceStr == "CrabNebula")  TevReg = 1.8;    // There are large offsets in some CrabNebula observations
    
///////////////////////
    //////////////////// Read TevCat coords from file
    TString gfilename="/data/magic/users-ifae/edosouto/MAGIC/BigSearch/Code/MethodExpect/TevCatON/TevCatdeg.txt";     
    int nl = 252;      // The number of lines in TevCatNameTypecomma.txt

    vector<Double_t> RaTev;                                 // Vectors for storing the grid point's coordinates'
    vector<Double_t> DecTev;

    ifstream runs;
    char sep;   
    // If we use another grid map of unknown number of lines
//     if (gfilename != "/data/magic/users-ifae/edosouto/MAGIC/BigSearch/Code/GridCoords_0.csv")    
//     {
//         char dummy_line[1000];
//         Int_t nlines = 0;
//         runs.open(gfilename);
//         while(!runs.eof()){
//             runs.getline(dummy_line, 1000, '\n');
//             nlines+=1;
//         }
//         nlines-=1;
//         runs.close();
//         nll = nlines;
//     }
    Double_t Coords[2];

    // Opening grid file and selecting which points are relevant for this observation
    runs.open(gfilename);
    for(Int_t i=0; i<nl; i++){
        runs>>Coords[0];     // R.A. coordinate in degrees
        runs>>sep;
        runs>>Coords[1];     // Dec coordinate in degrees, from 90 to -90

        // Great-circle distance from wikipedia gives the angle between two coordinates/vectors
        Double_t deltasig = acos(sin(Coords[1]*PI/180.)*sin(DecTel0*PI/180.)+cos(Coords[1]*PI/180.)*cos(DecTel0*PI/180.)*cos((RaTel0_deg-Coords[0])*PI/180.))*180./PI;
//         Double_t condition = 1.1*dispCut+ExclReg;     // I make the condition less stringent to be sure no grid point is left out
        // (1.8+1.4)*2 = 6.4 (We add some leeway)
//         Double_t condition = 6.6;     // The condition has to take into account the wobble changes. And CrabNebula has some 1.8 wobble deviation...
        Double_t condition = (TevReg+dispCut)*2+0.5;     // The condition has to take into account the wobble changes. And CrabNebula has some 1.8 wobble deviation...

        if(deltasig<=condition){
            RaTev.push_back(Coords[0]*24.0/360.);          // Save RA coords as hours for easier use
            DecTev.push_back(Coords[1]);
            std::cout << "TeV source at Ra/Dec = " << Coords[0] << " " << Coords[1] << std::endl; 
        }
    }
    runs.close();
///////////////////////
    
    /// Now we calculate the acceptance of the camera that has to be excised
    // I open the file with the Montecarlo maps of "pre-expected" number of events
    TFile *maps = new TFile(Form("/data/magic/users-ifae/edosouto/MAGIC/BigSearch/Melibea/%s/AcceptanceMaps_%s.root", period, period), "READ");
    
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
    Int_t histonum = 0;
    Double_t ValMC = 0.;
    Double_t phi0 = 0.;
    Double_t xp = 0., yp = 0.;
    std:vector<double> vcoords;

    std::cout.unsetf( std::ios::floatfield );  
    std::cout << std::fixed;
    std::cout.precision(3);
    std::cout << "t0 and tf = " << t0 << " " << tf << std::endl;
    
    TNtupleD *camcoor = new TNtupleD("camcoor", "CameraCoords", "Dec0:HA0:ZD:AZ:MJD:LST");
    TNtupleD *lc = new TNtupleD("lc", "FullCamLC", "Time:Acc");
    std::vector<Double_t> wobtimes;
    int runnum = 0;    
    //wobtimes.push_back(t0);     // We make it start at t0 --> Bad idea, we make so many cuts that this might not match with                                   // our first selected event's time
    
    // Get all the events in the camera (minus cuts) while we get the times where there're wobble 
    // changes
    Double_t tprev = 0;
    Double_t ppAz = 0;
    int flag = 0; // make a flag variable to check if we have to break or not
    for (int i = 0; i < nEntries; i++){
        data->GetEntry(i);
        tevent = mtime->GetMjd();
        ppAz = poPos->GetCorrAz();
        ppAz = (ppAz < 0.0 ? 360.0 + ppAz:ppAz);
        if ( tevent > tprev + 1./(24.*60.*60.)){
            camcoor->Fill(poPos->GetCorrDec(), poPos->GetCorrHa(), poPos->GetCorrZd(), ppAz, tevent, localsiderealT(mtime->GetGmst()));
            tprev = tevent;
        }

        if (hadro->GetHadronness() > hadCut) continue;
        if (poPos->GetCorrZd() >= 35.) continue;
        
        Double_t en = stereoDisp->GetEnergy();  
        if (en < 70 || en > 10000) continue;
       
        xDisp = stereoDisp->GetDirectionX();
        yDisp = stereoDisp->GetDirectionY();
        if ( TMath::IsNaN(xDisp) == 1 || TMath::IsNaN(yDisp) == 1 || xDisp*xDisp+yDisp*yDisp > dispCut*dispCut ) continue;
        
        // To calculate the position of the sources in the camera's FoV
        Double_t xi = 0., yi = 0.;
        HA0_i = poPos->GetCorrHa();
        Dec0_i = poPos->GetCorrDec();
        LST = localsiderealT(mtime->GetGmst());     // We calculate the LST with every event, instead of waiting any number of seconds/minutes
        
        //std::cout << data->GetTreeNumber() <<" "<< wobtimes.size() <<" "<< dummyt << " " << runnum << std::endl;
        // Sometimes the data is bad and the first run doesn't contain the first event that passes the cuts
        // If that happens, the wobble vector has false entries, so we do this instead
        if (wobtimes.size() == 0) {
            wobtimes.push_back(tevent);
            dummyt = tevent; 
            runnum = data->GetTreeNumber();
        }       
        if (runnum != data->GetTreeNumber()){   
            runnum = data->GetTreeNumber();    
            wobtimes.push_back(dummyt);     // We save where the run ended
            wobtimes.push_back(tevent);     // We save where the new run begins
        }
        dummyt = tevent;
        
        // Rotation part for the event's coordinates
        phi0 = ppAz - 120.;
        phi0 *= PI/180.; // Now in radians
        
        // We need to project the TevCat coordinates onto the camera to compare the x/y disp
        flag = 0; // make the flag = 0 to check the closeness
        ValMC = 0.;
        for (int j = 0; j < RaTev.size(); j++)
        {
            starcamtrans->Cel0CelToCam(Dec0_i, HA0_i, DecTev.at(j), LST-RaTev.at(j), xi, yi);     // LST and RaTev are in hours, don't worry
            xi*=mm2deg;            //In stereopardisp x and y are saved in deg, so it's easier to transform the point of the grid and not each event
            yi*=mm2deg;    
           
            Double_t RelPos = (xDisp-xi)*(xDisp-xi)+(yDisp-yi)*(yDisp-yi);  // Distance between the event and the rejected area.
            if (RelPos < ExclReg*ExclReg){        // Event too close to a source. Discarded
                //std::cout << "Point at" << RaTev.at(j)*360./24. << " " << DecTev.at(j) << " " << std::endl;
                flag = 1;       // Setting flag to 1 and use break
                break;          // We don't need to compare the event's position with other coordinates,it's already discarded.
            }
            
            // This is for the MC calculation. We only calculate this when the event is not discarded, because we save them together
            xp = xi*cos(phi0) - yi*sin(phi0);
            yp = xi*sin(phi0) + yi*cos(phi0);
            //if (DecTev.at(j) == -5.5 && RaTev.at(j) < 18.682 && RaTev.at(j) > 18.681) ExclReg = 0.41; // The coordinates of HESSJ1841-055
            //else ExclReg = 0.2;
            vcoords.push_back(xp);
            vcoords.push_back(yp);
            
        }
        if(flag == 1) {         // flag will be 1 only when we have encountered a break in nested loop
            //cout << "element found" << endl;
            vcoords.clear();
            continue;           // We don't want to exit the outer for loop, only to not save this event
        }
        if (RaTev.size() < 2) ValMC = Asearch(hmaps[histonum], xp, yp, ExclReg);
        else ValMC = Asearch(hmaps[histonum], vcoords, ExclReg);
        vcoords.clear();
        //std::cout << i << " " << j << " ValMC:" << ValMC << " ExclReg = " << ExclReg << " " << std::endl;
        lc->Fill(tevent, ValMC);
    }
    wobtimes.push_back(tevent); // For the last event
   
    // Open lightcurves file to add this histogram
    TFile *f = new TFile(LCfile, "UPDATE");
    lc->Write("FullCamLC", TObject::kOverwrite);
    camcoor->Write("CameraCoords", TObject::kOverwrite);
    f->WriteObject(&wobtimes, "WobbleTimes","Overwrite");
    f->Close();
    
    std::cout << "Finished adding the camera light curve to " <<  LCfile << std::endl;
    return 0;
}
