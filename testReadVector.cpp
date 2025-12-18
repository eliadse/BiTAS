#include "TFile.h"
#include "TVectorD.h"
#include <iostream>
#include <stdlib.h>

int testReadVector(){

    TFile *f = new TFile("../../Melibea/ST.03.03/AcceptanceMaps_ST.03.03.root", "READ");
    //TFile *f = new TFile("../../Melibea/ST.03.03/test.root", "READ");
    TVectorD *vres = (TVectorD*)f->Get("EpsilonRef");
    //TVectorD *vres = (TVectorD*)f->Get("NewVec");
    if (vres) {
        for (int i = 0; i < vres->GetNoElements(); i++) {
            std::cout << "v[" << i << "] = " << (*vres)[i] << std::endl;
        }
    }
    f->Close();
    return 0;
}
