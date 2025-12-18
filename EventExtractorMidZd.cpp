#include <TF1.h>
#include <TMarker.h>
#include <TH1D.h>
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
#include <TGraph.h>
#include "TVirtualFitter.h"
#include "TFitResult.h"
#include "TNtupleD.h"
#include <TLatex.h>
#include <TVector.h>

#include "MEnergyEst.h"
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

#include <utility>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <vector>
#include <iomanip>
#include <cstdlib>
#include <sstream>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <MHadronness.h>

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

int EventExtractorMidZd(const char *filename, const double dispCut, const char *period){
    Double_t PI = TMath::Pi();
    
    TChain* c = new TChain("Events");
    c->Add(filename);
    Int_t numEntries = c->GetEntries();
    if (numEntries < 100) return -1;

    TString fname = filename;                                               // We extract info from the file's name
    TString runStr = fname(fname.Index("_050")+1,8);                        // All subruns after the stereo upgrade start by 050
    int sourcestrnum = fname.Index("_Q_")+3;
    TString sourceStr = fname(sourcestrnum, fname.Index("-W")-sourcestrnum);
    if (sourceStr.Length() == 0) sourceStr = fname(sourcestrnum, fname.Index(".root")-sourcestrnum);
    TString nightStr = fname(fname.Index("_050")-8,4)+"_"+fname(fname.Index("_050")-4,2)+"_"+fname(fname.Index("_050")-2,2);

    const char *run = runStr;
    const char *source = sourceStr;
    const char *night = nightStr;
 
    std::cout << "Source:  " << source << "\n";
    std::cout << "Run:  " << run << "\n";
    
    // Definition and initialization of containers 
    const MGeomCamMagicTwo geom;    // For some reason these two can't be pointers.
    const MObservatory observ;        // No * for them

    MTime *mtime                   = new MTime();                                // I create the objects to access the methods
    MHadronness *hadro             = new MHadronness();                        // and all the Getters to get the info from the tree
    MStereoPar *stereoDisp         = new MStereoPar();
    MPointingPos *poPos            = new MPointingPos();
    MStarCamTrans *starcamtrans    = new MStarCamTrans(geom, observ);
    
    Double_t t0 = 0., tf = 0., tevent = 0.;                  // Declaration of variables that will be used. X and
    Double_t xEvent = 0., yEvent = 0.;                       // Y position in the camera of an event. The time of
    Double_t HA0_i = 0., Dec0_i = 0.;                        // the event in MJD. Its hadronness. The theoretical
    Double_t ppZd = 0., ppAz = 0.;
    Double_t EndEvent = 0.;
    Double_t RaTel0 = 0., HaTel0 = 0., DecTel0 = 0.;         // For the original pointing of the telescope
    Double_t had = 0.;                                       
    Double_t LST = 0.;                                       // Local sidereal time
    
    vector<Double_t> RaGrid;                                 // Vectors for storing the grid point's coordinates'
    vector<Double_t> DecGrid;
    vector<Double_t> IndexGrid;
    
    // Definition of used cuts and conversions
    Double_t mm2deg = geom.GetConvMm2Deg();                 // Conversion factor between mm and degrees in the camera
    //TODO make the correct hadronness cut
    Double_t hadCut = 0.3;
    Double_t SearchRadius = 0.15;                            // [deg] this is the PSF right now    
    
    // Setting status of branches and the addresses to their corresponding objects
    c->SetBranchStatus("*", 0);
    c->SetBranchStatus("MTime_1.*",1);
    c->SetBranchStatus("MStereoParDisp.*",1);
    c->SetBranchStatus("MPointingPos_1.*",1);
    c->SetBranchStatus("MHadronness.fHadronness",1);

    c->SetBranchAddress("MTime_1.", &mtime);
    c->SetBranchAddress("MStereoParDisp.", &stereoDisp);
    c->SetBranchAddress("MPointingPos_1.", &poPos);
    c->SetBranchAddress("MHadronness.", &hadro);

    // Get run's first values
    c->GetEntry(0);
    t0 = mtime->GetMjd();
    std::cout << std::setprecision(18) << "t0 = " <<  t0 << "\n";
    
    // Original pointing positions to obtain the reduced grid
    DecTel0 = poPos->GetCorrDec();                    // Pointing of the camera center in Ha/Dec (already corrected
    HaTel0 = poPos->GetCorrHa();
    
    LST = localsiderealT(mtime->GetGmst());    
    RaTel0 = LST-HaTel0;
    
    Double_t RaTel0_deg = RaTel0*360.0/24.0;        // Now we have the R.A. in degrees, like in the grid file
    if (HaTel0 == 0) {
        RaTel0_deg = poPos->GetRa()*360.0/24.0;
    }
    
    c->GetEntry(numEntries-1);
    tf = mtime->GetMjd();
    std::cout << std::setprecision(18) << "tf = " << tf << "\n";

    //////////////////// Read grid points from file
    //TString gfilename="/data/magic/users-ifae/edosouto/MAGIC/BigSearch/Code/GridCoords_Excluded.csv";
    //int nl = 588879;
    //TString gfilename="/data/magic/users-ifae/edosouto/MAGIC/BigSearch/Code/MARS3version/GridCoords_ExcludedMore.csv";     // MAGIC can see down to -20 Dec
    TString gfilename="/data/magic/users-ifae/edosouto/MAGIC/BigSearch/Code/MARS3version/NoOverlapGrid1.csv";     // MAGIC can see down to -20 Dec
    //int nl = 588080;      // If we're using GridCoords_ExcludedMore.csv this is the number of lines in the file
    int nl = 338005;      // If we're using NoOverlapGrid1.csv this is the number of lines in the file

    std::cout << "Using Grid file: " << gfilename << " If you want to use a new one change in EventExtractor.cpp" << "\n";
    ifstream runs;

    char sep;   
    Double_t Coords[3]; // To store the grid coordinates as we read them

    // Opening grid file and selecting which points are relevant for this observation
    runs.open(gfilename);
    for(Int_t i=0; i<nl; i++){
        runs>>Coords[0];     // Number of grid point
        runs>>sep;
        runs>>Coords[1];     // R.A. coordinate in degrees
        runs>>sep;
        runs>>Coords[2];     // Dec coordinate in degrees, from 90 to -90

        // Great-circle distance from wikipedia gives the angle between two coordinates/vectors
        // deltasigma = arccos(sin(phi1)*sin(phi2)+cos(phi1)*cos(phi2)*cos(deltalambda))
        
        Double_t deltasig = acos(sin(Coords[2]*PI/180.)*sin(DecTel0*PI/180.)+cos(Coords[2]*PI/180.)*cos(DecTel0*PI/180.)*cos((RaTel0_deg-Coords[1])*PI/180.))*180./PI;
        Double_t condition = 2.5*dispCut;     // I make the condition less stringent to be sure no grid point is left out

        if(deltasig<=condition){
            RaGrid.push_back(Coords[1]*24.0/360);          // Save RA coords as hours for easier use
            DecGrid.push_back(Coords[2]);
            IndexGrid.push_back(int(Coords[0]));
        }
    }
    runs.close();
    std::cout << "Ra grid size: " << RaGrid.size() << std::endl;

    TNtupleD *data = new TNtupleD("data","data","time:RAGrid:DecGrid:IndexGrid:Zd:GridX:GridY");
    bool filled = false;

    Double_t phi0 = 0., xp = 0., yp = 0.;
    
    for (int i = 0; i < numEntries; i++){
        c->GetEntry(i);

        // Apply hadronness cut
        had = hadro->GetHadronness();
        if(had > hadCut)     continue;
  
        xEvent = stereoDisp->GetDirectionX();
        yEvent = stereoDisp->GetDirectionY();
        EndEvent = stereoDisp->GetEnergy();
        
        // Apply disp cut
        if ( TMath::IsNaN(xEvent) == 1 || TMath::IsNaN(yEvent) == 1 || xEvent*xEvent+yEvent*yEvent > dispCut*dispCut) continue;        
        if (EndEvent < 200. || EndEvent > 10000.) continue;
//        if (EndEvent > 10000.) continue; // This Energy cut applies to low and medium Zd ranges, so we put it here
//        if (EndEvent < 200.) continue;

        // Apply estimated energy cut      

        tevent = mtime->GetMjd();        
        HA0_i = poPos->GetCorrHa();
        Dec0_i = poPos->GetCorrDec();
        ppZd = poPos->GetCorrZd();
        ppAz = poPos->GetCorrAz();
        ppAz = (ppAz < 0.0 ? 360.0 + ppAz:ppAz);    // Sometimes the Az is negative, so we fix it here already
        
        if (ppZd < 35) continue;
        if (ppZd > 50) continue;
        // Rotation part for the event's coordinates
        phi0 = ppAz - 120.;
        phi0 *= PI/180.; // Now in radians
        
        xp = xEvent*cos(phi0) - yEvent*sin(phi0);
        yp = xEvent*sin(phi0) + yEvent*cos(phi0);
        
        Double_t xi = 0., xip = 0.;
        Double_t yi = 0., yip = 0.;
        LST = localsiderealT(mtime->GetGmst());     // We calculate the LST with every event, instead of waiting any number of seconds/minutes
        
        // Projection of the grid back onto the camera
        for (int j = 0; j < RaGrid.size(); j++)
        {
            starcamtrans->Cel0CelToCam(Dec0_i, HA0_i, DecGrid.at(j), LST-RaGrid.at(j), xi, yi);     // LST and RaGrid are in hours, don't worry'
            xi*=mm2deg;            //In stereopardisp x and y are saved in deg, so it's easier to transform the point of the grid and not each event
            yi*=mm2deg;
            
            // Rotation part for the projection of the RA/Dec point on the camera
            xip = xi*cos(phi0) - yi*sin(phi0);
            yip = xi*sin(phi0) + yi*cos(phi0);
            
            if (xip*xip+yip*yip > (dispCut-SearchRadius)*(dispCut-SearchRadius)) continue;        // If this projected point is further or equal to dispCut from the center we don't take it
            
            Double_t RelPos = (xp-xip)*(xp-xip)+(yp-yip)*(yp-yip);
            
            if (RelPos < SearchRadius*SearchRadius){        // Less or Less equal? Not sure
                
                data->Fill(tevent, RaGrid.at(j), DecGrid.at(j), IndexGrid.at(j), ppZd, xip, yip);
                filled = true;
            }
            else continue;
        }
        
       if (i%100==0) std::cout << "Still alive  " << i <<  "\n";
    }
if (filled) {
    TFile *f0m = new TFile(Form("/data/magic/users-ifae/edosouto/MAGIC/BigSearch/Code/MARS3version/%s/35to50/EL_%s_%s_%s_%.1f.root", period, night, source, run, dispCut),"recreate");
    std::cout << "f0 medium name:  " << f0m->GetName() << "\n";
    data->Write("", TObject::kOverwrite);
    f0m->Write();
    f0m->Close();
    delete f0m;
}
    delete data;

    return 0;
}
