#include "TFile.h"
int testWriteVector(){

    //TFile *f = new TFile("../../Melibea/ST.03.03/test.root","RECREATE");
    TFile *f = new TFile("../../Melibea/ST.03.03/test.root","UPDATE");
    TVectorD v(2); v[0]=1; v[1]=2;
    f->cd();
    //v.Write("NewVec", TObject::kOverwrite);
    v.Write("MyVec", TObject::kOverwrite);
    f->Close();
    return 0;
}
