#include "GlauberEvent.cxx"
#include "TFile.h"
void RunGlauberSimulation(int N, const int A, double_t a, double_t NN_crossSection, const char *fName, Option_t *option){
    TFile *f = new TFile(fName, "UPDATE");

    GlauberEvent *g = new GlauberEvent(A, a, NN_crossSection, option);
    TString opt = option;
    opt.ToLower();
    Bool_t progress = kFALSE;
    Bool_t hists = kFALSE;
    Bool_t trees = kFALSE;
    if(opt.Contains("progress")){progress = kTRUE;}
    if(opt.Contains("hists")){hists = kTRUE;}
    if(opt.Contains("trees")){trees = kTRUE;}
    for(int i = 0; i < N; i++){
        g->Collide();
        if(progress){
            if(i%(N/100) == 0) {cout << i/(N/100)<< "%" << endl;
        }
        }
    }
    if(progress){
        cout << "100%" << endl;
        cout << N << " collisions done!" << endl;
    }
    if(hists){
    TH1D *bHist = new TH1D();
    bHist = (TH1D*)g->GetbHist();
    TH1D *Npart = new TH1D();
    Npart = (TH1D*)g->GetNpartHist();
    TH1D *Ncoll = new TH1D();
    Ncoll = (TH1D*)g->GetNcollHist();
    }
    if(trees){
        TTree *GlauberTree = new TTree();
        GlauberTree = (TTree*)g->GetTree();
    }
    f->Write();
    f->Close();
}