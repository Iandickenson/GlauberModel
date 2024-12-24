#include "GlauberEvent.cxx"

void RunGlauberSimulation(int N, const int A, double_t a, double_t NN_crossSection, Option_t *option){
    GlauberEvent *g = new GlauberEvent(A, a, NN_crossSection, "makehists");
    for(int i = 0; i < N; i++){
        g->Collide();
    }
    TH1D *bHist = new TH1D();
    bHist = (TH1D*)g->GetbHist();
    TH1D *Npart = new TH1D();
    Npart = (TH1D*)g->GetNpartHist();
    TH1D *Ncoll = new TH1D();
    Ncoll = (TH1D*)g->GetNcollHist();
}