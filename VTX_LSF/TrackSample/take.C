#include "TROOT.h"
#include "TChain.h"
#include "TTree.h"
#include "iostream"
#include "fstream"

void take()
{
    TChain *c = new TChain("TrackInfo");
    c->Add("/mnt/f/ALG_Learning/VTX_LSF/TrackSample/*root");

    // Declaration of leaf types
    Int_t event;
    Int_t run;
    Double_t truth_Ks_x;
    Double_t truth_Ks_y;
    Double_t truth_Ks_z;
    Double_t V_truth_Ks[4];
    Int_t index;
    Double_t HelixInfoPionp[10]; //[index]
    Double_t HelixE1p[10];       //[index]
    Double_t HelixE2p[10];       //[index]
    Double_t HelixE3p[10];       //[index]
    Double_t HelixE4p[10];       //[index]
    Double_t HelixE5p[10];       //[index]
    Double_t HelixInfoPionm[10]; //[index]
    Double_t HelixE1m[10];       //[index]
    Double_t HelixE2m[10];       //[index]
    Double_t HelixE3m[10];       //[index]
    Double_t HelixE4m[10];       //[index]
    Double_t HelixE5m[10];       //[index]

    // List of branches
    TBranch *b_event;          //!
    TBranch *b_run;            //!
    TBranch *b_truth_Ks_x;     //!
    TBranch *b_truth_Ks_y;     //!
    TBranch *b_truth_Ks_z;     //!
    TBranch *b_V_truth_Ks;     //!
    TBranch *b_index;          //!
    TBranch *b_HelixInfoPionp; //!
    TBranch *b_HelixE1p;       //!
    TBranch *b_HelixE2p;       //!
    TBranch *b_HelixE3p;       //!
    TBranch *b_HelixE4p;       //!
    TBranch *b_HelixE5p;       //!
    TBranch *b_HelixInfoPionm; //!
    TBranch *b_HelixE1m;       //!
    TBranch *b_HelixE2m;       //!
    TBranch *b_HelixE3m;       //!
    TBranch *b_HelixE4m;       //!
    TBranch *b_HelixE5m;       //!

    c->SetBranchAddress("event", &event, &b_event);
    c->SetBranchAddress("run", &run, &b_run);
    c->SetBranchAddress("truth_Ks_x", &truth_Ks_x, &b_truth_Ks_x);
    c->SetBranchAddress("truth_Ks_y", &truth_Ks_y, &b_truth_Ks_y);
    c->SetBranchAddress("truth_Ks_z", &truth_Ks_z, &b_truth_Ks_z);
    c->SetBranchAddress("V_truth_Ks", V_truth_Ks, &b_V_truth_Ks);
    c->SetBranchAddress("index", &index, &b_index);
    c->SetBranchAddress("HelixInfoPionp", HelixInfoPionp, &b_HelixInfoPionp);
    c->SetBranchAddress("HelixE1p", HelixE1p, &b_HelixE1p);
    c->SetBranchAddress("HelixE2p", HelixE2p, &b_HelixE2p);
    c->SetBranchAddress("HelixE3p", HelixE3p, &b_HelixE3p);
    c->SetBranchAddress("HelixE4p", HelixE4p, &b_HelixE4p);
    c->SetBranchAddress("HelixE5p", HelixE5p, &b_HelixE5p);
    c->SetBranchAddress("HelixInfoPionm", HelixInfoPionm, &b_HelixInfoPionm);
    c->SetBranchAddress("HelixE1m", HelixE1m, &b_HelixE1m);
    c->SetBranchAddress("HelixE2m", HelixE2m, &b_HelixE2m);
    c->SetBranchAddress("HelixE3m", HelixE3m, &b_HelixE3m);
    c->SetBranchAddress("HelixE4m", HelixE4m, &b_HelixE4m);
    c->SetBranchAddress("HelixE5m", HelixE5m, &b_HelixE5m);

    std::ofstream out1("TrackP.dat", std::ios::app);
    std::ofstream out2("TrackM.dat", std::ios::app);
    std::ofstream out3("TrackEP.dat", std::ios::app);
    std::ofstream out4("TrackEM.dat", std::ios::app);
    for (int i = 0; i < c->GetEntries(); i++)
    {
        c->GetEntry(i);
        out1 << i << "\t";
        for (int j = 0; j < 5; j++)
        {
            out1 << HelixInfoPionp[j] << "\t";
        }
        out1 << std::endl;

        out2 << i << "\t";
        for (int j = 0; j < 5; j++)
        {
            out2 << HelixInfoPionm[j] << "\t";
        }
        out2 << std::endl;

        // Helix Ep
        out3 << i << "\t";
        for (int j = 0; j < 5; j++)
        {
            out3 << HelixE1p[j] << "\t";
        }
        out3 << std::endl;

        out3 << i << "\t";
        for (int j = 0; j < 5; j++)
        {
            out3 << HelixE2p[j] << "\t";
        }
        out3 << std::endl;

        out3 << i << "\t";
        for (int j = 0; j < 5; j++)
        {
            out3 << HelixE3p[j] << "\t";
        }
        out3 << std::endl;

        out3 << i << "\t";
        for (int j = 0; j < 5; j++)
        {
            out3 << HelixE4p[j] << "\t";
        }
        out3 << std::endl;

        out3 << i << "\t";
        for (int j = 0; j < 5; j++)
        {
            out3 << HelixE5p[j] << "\t";
        }
        out3 << std::endl;

        // Helix Em
        out4 << i << "\t";
        for (int j = 0; j < 5; j++)
        {
            out4 << HelixE1m[j] << "\t";
        }
        out4 << std::endl;

        out4 << i << "\t";
        for (int j = 0; j < 5; j++)
        {
            out4 << HelixE2m[j] << "\t";
        }
        out4 << std::endl;

        out4 << i << "\t";
        for (int j = 0; j < 5; j++)
        {
            out4 << HelixE3m[j] << "\t";
        }
        out4 << std::endl;

        out4 << i << "\t";
        for (int j = 0; j < 5; j++)
        {
            out4 << HelixE4m[j] << "\t";
        }
        out4 << std::endl;

        out4 << i << "\t";
        for (int j = 0; j < 5; j++)
        {
            out4 << HelixE5m[j] << "\t";
        }
        out4 << std::endl;
    }
    out1.close();
    out2.close();
    out3.close();
    out4.close();
}