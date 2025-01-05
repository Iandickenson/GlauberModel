#include <cmath>
#include <stdio.h>
#include <iostream>
#include "TH1D.h"
#include "TF1.h"
#include "../inc/GlauberEvent.h"
#include "../inc/GlauberSimulation.h"
using namespace std;
//Default constructor 
GlauberSimulation::GlauberSimulation(){
    gMu = 0;
    gK = 0;

    g_bHist = NULL;
    g_Ncoll = NULL;
    g_Npart = NULL;

    gMultiplicityHist = NULL;
}
//Destructor
GlauberSimulation::~GlauberSimulation(){

}
//Optional Constructor
GlauberSimulation::GlauberSimulation(double_t mu, double_t k, TH1D *ncoll, TH1D *npart, TH1D *bhist){
    gMu = mu;
    gK = k;

    g_bHist = bhist;
    g_Ncoll = ncoll;
    g_Npart = npart;

    gMultiplicityHist = NULL;
    gNegativeBinomial = NULL;
    /*if ( ((TArrayD *)dataHisto->GetSumw2())->GetSize() == 0){
    dataHisto->Sumw2();
    }*/
}
//Optional Constructor - pass GlauberEvent instead of individual hists: this will
// help keep track of what experimental settings are being simulated 
GlauberSimulation::GlauberSimulation(double_t mu, double_t k, GlauberEvent *glauberEvent){
    gMu = mu;
    gK = k;

    g_bHist = glauberEvent->GetbHist();
    g_Ncoll = glauberEvent->GetNcollHist();
    g_Npart = glauberEvent->GetNpartHist();

    gMultiplicityHist = new TH1D("gMultiplicityHist", "Multiplicity Distribution", 1500, 0, 1500);

    gNegativeBinomial = new TF1("NBD", "(TMath::Gamma(x + [0]) / (TMath::Gamma(x + 1)*TMath::Gamma([0]))) * ((TMath::Power([1]/[0], x))/(TMath::Power(([1]/[0]) + 1, x + [0])))", 0, 10);
    gNegativeBinomial->SetParameter(0, gK);
    gNegativeBinomial->SetParameter(1, gMu);
    gNegativeBinomial->SetNpx(500);
}
void GlauberSimulation::SetNegativeBinomialParameters(double_t mu, double_t k){
    gMu = mu;
    gK = k;
    gNegativeBinomial->SetParameter(0, gK);
    gNegativeBinomial->SetParameter(1, gMu);
}
//Generate multiplicity by drawing from negative binomial distribution with current
// mu and k parameters
void GlauberSimulation::GenerateMultiplicity(){
    maxNcoll = g_Ncoll->GetNbinsX();
    for(int ncollIterator = 1; ncollIterator < maxNcoll; ncollIterator++){
        // Remember that the histogram bins are offset by 1 so bin 1 corresponds
        // to a value of 0 collisions not 1 collision
        collisions = ncollIterator - 1;
        binContent = g_Ncoll->GetBinContent(ncollIterator);
        for(int binIterator = 0; binIterator < binContent; binIterator++){
            thisMultiplicity = 0;
            for(int collisionIterator = 0; collisionIterator < collisions; collisionIterator++){
                thisMultiplicity += gNegativeBinomial->GetRandom();
            }
            gMultiplicityHist->Fill(thisMultiplicity);
        }
    }
}
//Option to generate multiplicity with new mu and k parameters
void GlauberSimulation::GenerateMultiplicity(double_t mu, double_t k){
    SetNegativeBinomialParameters(mu, k);
    GenerateMultiplicity();
}
//Normalizing the multiplicity hist
void GlauberSimulation::Normalize(){

}
double_t GlauberSimulation::ChiSquared(TH1D *dataHist){
    return 5.;
}