#include "TFile.h"
#include "TH1.h"
#include "fstream"
#include "iostream"

void take()
{
    TFile *f = new TFile("makemass_4180_data.root");
    TH1F *h = (TH1F *)f->Get("Ds[2]");
    std::ofstream out("dataSample.dat", std::ios::app);
    for (int i = 1; i <= h->GetNbinsX(); i++)
    {
        out << i << "\t" << h->GetBinCenter(i) << "\t" << h->GetBinContent(i) << "\t" << h->GetBinError(i) << std::endl;
    }
    out.close();
}