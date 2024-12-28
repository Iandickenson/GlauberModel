#include <cmath>
#include <stdio.h>
#include <iostream>
#include "TH1D.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TTree.h"
#include "TGraph.h"
#include "../inc/GlauberEvent.h"
using namespace std;
//Default constructor
GlauberEvent::GlauberEvent(){
    g_nNucleons = 0;
    gNN_crossSection = 0;
    gSkinDepth = 0;
    gCollision_d = 0;
    gNuclearRadius = 0;
    

    gRandom = new TRandom3();

    makehists = kFALSE;

    g_bHist = NULL;
    g_Npart = NULL;
    g_Ncoll = NULL;
    gTree = NULL;
}

//Destructor
GlauberEvent::~GlauberEvent(){
    delete gRandom;
    if(makehists){
        delete g_bHist;
        delete g_Npart;
        delete g_Ncoll;
    }
    if(maketrees){
        delete gTree;
    }
    // delete all pointers allocated with new
}

//Optional Constructor
GlauberEvent::GlauberEvent(const int A, double_t a, double_t NN_crossSection, Option_t *option){
    // Nucleus quantities
    g_nNucleons = A;
    gNN_crossSection = NN_crossSection;
    gSkinDepth = a;

    gCollision_d = pow(gNN_crossSection/(10*M_PI), 0.5); // Condition for collisions
                                                         // between nucleons (factor
                                                         // of 10 in denominator is the
                                                         // conversion from mb to fm^2)

    // Positions of nucleons in nuclei
    g_xNucleusA.resize(g_nNucleons);
    g_xNucleusB.resize(g_nNucleons);
    g_yNucleusA.resize(g_nNucleons);
    g_yNucleusB.resize(g_nNucleons);
    // Placeholder spherical coordinates
    gNucleons_r.resize(g_nNucleons);
    g_phi.resize(g_nNucleons);
    g_theta.resize(g_nNucleons);
    // Random number generator
    gRandom = new TRandom3();
    gRandom->SetSeed(0);

    // Distributions
    gNuclearRadius = 1.07 * std::pow(A, 1./3.);
    gWoodsSaxonPDF = new TF1("gWoodSaxonPDF","x*x / (1.+ TMath::Exp((x-[0])/[1]))", 0 , 10);
    gWoodsSaxonPDF->SetParameter(0, gNuclearRadius);
    gWoodsSaxonPDF->SetParameter(1, gSkinDepth);

    gImpactParameterPDF = new TF1("gImpactParameterPDF", "2*TMath::Pi()*x", 0, 18);

    // Variables to be modified by other functions (i.e. )
    gParticipants = 0;
    gCollisions = 0;
    gImpactParameter = 0;
    gIsParticipant = false;
    gIsCollision = false;

    // Option to write histograms 
    TString opt = option;
    opt.ToLower();
    makehists = kFALSE;
    maketrees = kFALSE;
    if(opt.Contains("hists")){makehists = kTRUE;}
    if(opt.Contains("trees")){maketrees = kTRUE;}
    g_bHist = NULL;
    g_Npart = NULL;
    g_Ncoll = NULL;
    gTree = NULL;
    if(makehists){
        g_bHist = new TH1D("bhist", "Impact Parameter", 40, 0, 18);
        Int_t NpartRange = g_nNucleons*2 + 50;
        g_Npart = new TH1D("Npart", "Number of Participants", NpartRange, 0, NpartRange);
        g_Ncoll = new TH1D("Ncoll", "Number of Collisions", 2500, 0, 2500);

        g_bHist->SetTitle("Impact Parameter");
        g_bHist->GetXaxis()->SetTitle("b (fm)");
        g_bHist->GetYaxis()->SetTitle("Counts");

        g_Ncoll->SetTitle("Number of Collisions");
        g_Ncoll->GetXaxis()->SetTitle("Ncoll");
        g_Ncoll->GetYaxis()->SetTitle("Counts");

        g_Npart->SetTitle("Number of Participants");
        g_Npart->GetXaxis()->SetTitle("Npart");
        g_Npart->GetYaxis()->SetTitle("Counts");
    }
    if(maketrees){
        gTree = new TTree("gTree", "Glauber Tree");
        gTree->Branch("b", &gImpactParameter, "b/D");
        gTree->Branch("Participants", &gParticipants, "Participants/I");
        gTree->Branch("Collisions", &gCollisions, "Collisions/I");
    }
}

//GlauberEvent::WoodsSaxonsPDF(){

//}

void GlauberEvent::Collide(){
    // Randomly distribute the nucleons in each pair of nuclei and then loop through 
    // all nucleons to check collisions, participants, and fill impact parameter.
    // Sets values in member variables of GlauberEvent
    
    CollisionReset();
    
    vector<Bool_t> checksA;                                                   // ensure that we do not double 
    vector<Bool_t> checksB;                                                   // count participants 
    checksA.resize(g_nNucleons);
    checksB.resize(g_nNucleons);
    
    
    // Setting impact parameter 
    gImpactParameter = gImpactParameterPDF->GetRandom();

    // Filling each nucleus 
      // do some operations on g_xNucleusA, g_yNucleusB,...
      
    FillNucleus(g_xNucleusA, g_yNucleusA, "nuca");
    FillNucleus(g_xNucleusB, g_yNucleusB, "nucb");
    // Set impact parameter randomly and adjust nuclear positions
    
    for(int iteratorA = 0; iteratorA < g_nNucleons; iteratorA++){

        for(int iteratorB = 0; iteratorB < g_nNucleons; iteratorB++){
            gNucleonDistance = GetNucleonDistance(iteratorA, iteratorB, g_xNucleusA, g_xNucleusB, g_yNucleusA, g_yNucleusB);
            if(gNucleonDistance < gCollision_d){
                gIsCollision = true;
                gCollisions++;
                if(checksA[iteratorA] == false){
                    gParticipants++;
                    checksA[iteratorA] = true;
                }
                if(checksB[iteratorB] == false){
                    gParticipants++;
                    checksB[iteratorB] = true;
                }
            }
        }
    }

    if(makehists){

        if(gIsCollision){
            g_bHist->Fill(gImpactParameter);
        }
        g_Npart->Fill(gParticipants);
        g_Ncoll->Fill(gCollisions);
    }
    if(maketrees){
        gTree->Fill();
    }

}

void GlauberEvent::DrawCollision(){
    TGraph *RedParticipants = new TGraph();
    TGraph *RedSpectators = new TGraph(g_nNucleons, g_xNucleusA.data(), g_yNucleusA.data());

    TGraph *BlueParticipants = new TGraph();
    TGraph *BlueSpectators = new TGraph(g_nNucleons, g_xNucleusB.data(), g_yNucleusB.data());
    
    // Set nucleons as either spectators or participants
    for(int iRed = 0; iRed < g_nNucleons; iRed++){
        for(int jBlue = 0; jBlue < g_nNucleons; jBlue++){
            gNucleonDistance = GetNucleonDistance(iRed, jBlue, g_xNucleusA, g_xNucleusB, g_yNucleusA, g_yNucleusB);
            if(gNucleonDistance < gCollision_d){
                RedParticipants->AddPoint(g_xNucleusA[iRed], g_yNucleusA[iRed]);
                BlueParticipants->AddPoint(g_xNucleusB[jBlue], g_yNucleusB[jBlue]);
            }
        }
    }
    // For simplicity, we initalized the spectators as all nucleons (see above) so we 
    // can show the participants by just drawing over the spectators.
    RedSpectators->SetMarkerStyle(20);
    RedSpectators->SetMarkerSize(2);
    RedSpectators->SetMarkerColor(46);
    RedSpectators->GetXaxis()->CenterTitle();
    RedSpectators->GetXaxis()->SetTitle("x (fm)");
    RedSpectators->SetMinimum(-18);
    RedSpectators->SetMaximum(18);
    RedSpectators->GetXaxis()->SetLimits(-18,18);
    RedSpectators->GetYaxis()->CenterTitle();
    RedSpectators->GetYaxis()->SetTitle("y (fm)");
    RedSpectators->Draw("AP");

    BlueSpectators->SetMarkerStyle(20);
    BlueSpectators->SetMarkerSize(2);
    BlueSpectators->SetMarkerColor(38);
    BlueSpectators->Draw("P");

    RedParticipants->SetMarkerStyle(20);
    RedParticipants->SetMarkerSize(2);
    RedParticipants->SetMarkerColor(2);
    RedParticipants->Draw("P");

    BlueParticipants->SetMarkerStyle(20);
    BlueParticipants->SetMarkerSize(2);
    BlueParticipants->SetMarkerColor(4);
    BlueParticipants->Draw("P");
}

void GlauberEvent::CollisionReset(){
    gParticipants = 0;
    gCollisions = 0;
    gImpactParameter = 0;
    gIsParticipant = false;
    gIsCollision = false;
}

void GlauberEvent::FillNucleus(vector<double_t> &nucleusX, vector<double_t> &nucleusY, Option_t *optionb){
    
    TString whichNuc = optionb;
    whichNuc.ToLower();
    Bool_t isNucA = kFALSE;
    Bool_t isNucB = kFALSE;
    if(whichNuc.Contains("nuca")){
        isNucA = kTRUE;
    }
    else if(whichNuc.Contains("nucb")){
        isNucB = kTRUE;
    }
    else {
        cout <<"ERROR: GlauberEvent::FillNucleus - incorrect nucleus specified" <<endl;
        exit(EXIT_FAILURE);
    }
    
    for(int nucIterator = 0; nucIterator < g_nNucleons; nucIterator++){
        gNucleons_r[nucIterator] = gWoodsSaxonPDF->GetRandom();

        g_phi[nucIterator] = gRandom->Uniform(0, 2.*M_PI);
        g_theta[nucIterator] = acos(gRandom->Uniform(-1,1));

        if(isNucA){
        nucleusX[nucIterator] = gNucleons_r[nucIterator]*sin(g_theta[nucIterator])*cos(g_phi[nucIterator]) + (gImpactParameter/2);
        }
        if(isNucB){
        nucleusX[nucIterator] = gNucleons_r[nucIterator]*sin(g_theta[nucIterator])*cos(g_phi[nucIterator]) - (gImpactParameter/2);
        }
        nucleusY[nucIterator] = gNucleons_r[nucIterator]*sin(g_theta[nucIterator])*sin(g_phi[nucIterator]);
    }
}
/*
TH1D* GetbHist(){
    return g_bHist;
}
TH1D* GetNpartHist(){
    return g_Npart;
}
TH1D* GetNcollHist(){
    return g_Ncoll;
}*/

// Returns number of participants in current Glauber event
Int_t GlauberEvent::GetNparticipants(){
    return gParticipants;
}
// Returns number of collisions in current Glauber event
Int_t GlauberEvent::GetNcollisions(){
    return gCollisions;
}
// Returns impact parameter of current Glauber event
double_t GlauberEvent::GetImpactParameter(){
    return gImpactParameter;
}
double_t GlauberEvent::GetNucleonDistance(int nucAiterator, int nucBiterator, vector<double_t> &xNucleusA, vector<double_t> &xNucleusB, vector<double_t> &yNucleusA, vector<double_t> &yNucleusB){
    return pow((xNucleusB[nucBiterator]-xNucleusA[nucAiterator])*(xNucleusB[nucBiterator]-xNucleusA[nucAiterator]) + (yNucleusB[nucBiterator]-yNucleusA[nucAiterator])*(yNucleusB[nucBiterator]-yNucleusA[nucAiterator]), 0.5);
}
