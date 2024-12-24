#include "GlauberEvent.cxx"

void RunGlauberSimulation(int N, const int A, double_t a, double_t NN_crossSection, Option_t *option, Option_t *optiona){
    GlauberEvent *g = new GlauberEvent(A, a, NN_crossSection, "makehists");
    TString opt = optiona;
    opt.ToLower();
    Bool_t progress = kFALSE;
    if(opt.Contains("progress")){
        progress = kTRUE;
    }
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
    TH1D *bHist = new TH1D();
    bHist = (TH1D*)g->GetbHist();
    TH1D *Npart = new TH1D();
    Npart = (TH1D*)g->GetNpartHist();
    TH1D *Ncoll = new TH1D();
    Ncoll = (TH1D*)g->GetNcollHist();
}