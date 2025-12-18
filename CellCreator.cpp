#include <TCanvas.h>
#include <TChain.h>
#include <TSystem.h>
#include <TFile.h>
#include <TGraph.h>
#include <TKey.h>

#include <string>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <iomanip>      // std::setprecision
#include <cstdlib>
#include <sstream>
#include <stdlib.h>

int CellCreator(const char *Filename, const char *Night, const char *Source, const double dispCut){
   
    TString fn = Filename; 
    TString od = fn(0,fn.Index("/EL_"));
    const char *OutDir = od;
    std::cout << " Opening " << Filename << "\n";
    std::cout << "Out directory : " << OutDir << std::endl;

    TChain *c = new TChain("data");
    c->Add(Filename);

    Double_t timepoint, RAG, DecG, IndG, ZdEv;
    Double_t gridx, gridy;
   
    //Should be "time:RAGrid:DecGrid:IndexGrid:Zd:GridX:GridY:Wangle"
    c->SetBranchAddress("time",&timepoint);
    c->SetBranchAddress("RAGrid",&RAG);
    c->SetBranchAddress("DecGrid",&DecG);
    c->SetBranchAddress("IndexGrid",&IndG);
    c->SetBranchAddress("Zd", &ZdEv);   // Pointing pos coords
    c->SetBranchAddress("GridX", &gridx);     // The gridx or gridx are the expected rates variable 
    c->SetBranchAddress("GridY", &gridy);

    Int_t nEntries = c->GetEntries();
    if (nEntries < 100) {
        std::cout<< "The file is too small" <<std::endl;
        return -1;
    }
    
    //  Here we extract the grid-points used in this observation //
    std::vector<int> vec;
    for(int i = 0; i < nEntries; i++){
         c->GetEntry(i);
         vec.push_back(int(IndG));
    }

    //  Here we sort them in ascending order and remove the duplicated numbers so we can loop over them later //
    sort( vec.begin(), vec.end() );
    vec.erase( unique( vec.begin(), vec.end() ), vec.end() );
    
    //  Now vec contains the index numbers in an unique and orderly manner
    const int nPoints= vec.size();
    
    std::vector <Double_t> v_time[nPoints];
    std::vector <Double_t> v_RAEvent[nPoints];
    std::vector <Double_t> v_DecEvent[nPoints];
    std::vector <Double_t> v_Zd[nPoints];
    std::vector <Double_t> v_gridx[nPoints];
    std::vector <Double_t> v_gridy[nPoints];
 
    for(int k=0; k<nEntries; k++){
        c->GetEntry(k);
        for(int l=0; l < nPoints; l++){
            if (IndG == vec.at(l)) {
                //std::cout << IndG << "\n";
                v_time[l].push_back(timepoint);
                v_RAEvent[l].push_back(RAG);
                v_DecEvent[l].push_back(DecG);
                v_Zd[l].push_back(ZdEv);
                v_gridx[l].push_back(gridx);
                v_gridy[l].push_back(gridy);
            }
        }
    }
 
    Double_t timen, Zdn, RAn, Decn, gridxn, gridyn;
    Double_t timer, Zdr, RAr, Decr, gridxr, gridyr;
    Double_t timew, Zdw, RAw, Decw, gridxw, gridyw;
    
    if (gSystem->AccessPathName(Form("%s/LightCurves_%s_%s_%.1f.root",OutDir, Night, Source, dispCut))){
        
        TFile *f = new TFile(Form("%s/LightCurves_%s_%s_%.1f.root",OutDir, Night, Source, dispCut), "NEW");
        if ( f->IsOpen() ) printf("File opened successfully\n");
        std::cout << "We just created a file with name " << f->GetName() << "\n";
        
        TTree *t[nPoints];
        for(int ii = 0; ii<nPoints; ii++){
            t[ii] = new TTree;
            t[ii]->SetName(Form("t_%i", vec.at(ii)));
            t[ii]->SetTitle(Form("Tree containing Cell%i", vec.at(ii)));
                
            t[ii]->Branch("Time", &timen,"timen/D");
            t[ii]->Branch("RA", &RAn, "RAn/D");
            t[ii]->Branch("Dec", &Decn, "Decn/D"); 
            t[ii]->Branch("Zd", &Zdn, "Zdn/D");
            t[ii]->Branch("GridX", &gridxn, "gridxn/D"); 
            t[ii]->Branch("GridY", &gridyn, "GridYn/D"); 
                       
            for (unsigned int jj = 0; jj < v_time[ii].size(); jj++){
                timen = v_time[ii].at(jj);
                Decn = v_DecEvent[ii].at(jj);
                RAn = v_RAEvent[ii].at(jj);
                Zdn = v_Zd[ii].at(jj);
                gridxn = v_gridx[ii].at(jj);
                gridyn = v_gridy[ii].at(jj);
                t[ii]->Fill();
            }
            t[ii]->Write();
        }
        f->Write();
        f->Close();
    }
    else {      
        TFile *f = new TFile(Form("%s/LightCurves_%s_%s_%.1f.root",OutDir, Night, Source, dispCut), "UPDATE");
        std::vector <Int_t> index_read;
        
        // Black fucking box where we get the key names
        TIter next(f->GetListOfKeys());
        TKey *key;
        while ((key=(TKey*)next())) {
            //std::cout << "File contains " << key->GetClassName() << ": " << key->GetName() << " at " <<key->GetSeekKey() <<"\n";
            TString name = key->GetName();
            TString cosi = name(2,name.Sizeof()-1);
            Int_t index = atoi(cosi);
            index_read.push_back(index);
        } 
        // When we add new events a new cycle is created (ex. t_1;4 to t_1;5)
        // When we read the keys we only obtain the first part (t_1) but as many 
        // times as cycles it has so we need to remove the duplicated entries in index_read
        index_read.erase( unique( index_read.begin(), index_read.end() ), index_read.end() );
        
        for (unsigned int l = 0; l < vec.size(); l++){
            std::vector<int>::iterator it;

            it = find (index_read.begin(), index_read.end(), vec.at(l));
            std::cout << "it = "<< *it << "vec = " << vec.at(l) << "\n";
            
            if (it != index_read.end()){
                std::cout << "File already contains a tree for Cell " << *it << '\n';
                TTree *t = (TTree*)f->Get(Form("t_%i",*it));
                
                t->SetBranchAddress("Time", &timer);
                t->SetBranchAddress("Zd", &Zdr);
                t->SetBranchAddress("RA", &RAr);
                t->SetBranchAddress("Dec", &Decr);
                t->SetBranchAddress("GridX", &gridxr);
                t->SetBranchAddress("GridY", &gridyr);
                
                for (unsigned int m = 0; m < v_time[l].size(); m ++){
                    //std::cout << std::setprecision(9) << "time: " << v_time[l].at(m) << "\n";
                    timer = v_time[l].at(m);
                    Decr = v_DecEvent[l].at(m);
                    RAr = v_RAEvent[l].at(m);
                    Zdr = v_Zd[l].at(m);
                    gridxr = v_gridx[l].at(m);
                    gridyr = v_gridy[l].at(m);
                    t->Fill();
                }
                t->Write(0,TObject::kOverwrite);                
            }
            else{
                std::cout << "No tree found, creating a new one\n";
                TTree *t = new TTree;
                t->SetName(Form("t_%i", vec.at(l)));
                t->SetTitle(Form("Tree containing Cell%i", vec.at(l)));
                
                t->Branch("Time", &timew,"timew/D");
                t->Branch("Zd", &Zdw, "Zdw/D");
                t->Branch("RA", &RAw, "RAw/D");
                t->Branch("Dec", &Decw, "Decw/D");    
                t->Branch("GridX", &gridxw, "gridxw/D");      
                t->Branch("GridY", &gridyw, "gridyw/D");    

                for (unsigned int n = 0; n < v_time[l].size(); n ++){
                    timew = v_time[l].at(n);
                    Decw = v_DecEvent[l].at(n);
                    RAw = v_RAEvent[l].at(n);
                    Zdw = v_Zd[l].at(n);
                    gridxw = v_gridx[l].at(n);
                    gridyw = v_gridy[l].at(n);
                    t->Fill();
                }
                t->Write(0,TObject::kOverwrite);
            }
        }
    }
    return 0;
}
               
