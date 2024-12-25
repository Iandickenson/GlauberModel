
#ifndef GLAUBEREVENT_H
#define GLAUBEREVENT_H

class TRandom3;
class TF1;
class TH1D;
class TTree;

class GlauberEvent : public TObject{

    private:
        Int_t g_nNucleons;
        double_t gNN_crossSection;
        double_t gSkinDepth;
        double_t gCollision_d;
        double_t gNuclearRadius;

        double_t gImpactParameter;

        Int_t gParticipants;
        Int_t gCollisions;

        Bool_t gIsParticipant;
        Bool_t gIsCollision;

        std::vector<double_t> g_xNucleusA;
        std::vector<double_t> g_xNucleusB;
        std::vector<double_t> g_yNucleusA;
        std::vector<double_t> g_yNucleusB;

        std::vector<double_t> gNucleons_r; //Placeholder nucleon positions 
        std::vector<double_t> g_phi; // used to generate r, phi, theta, 
        std::vector<double_t> g_theta; // then fills x and y with corresponding 
                                       // cartesian coordinates

        TRandom3 *gRandom;

        TF1 *gWoodsSaxonPDF;
        TF1 *gImpactParameterPDF;

        TH1D *g_bHist;
        TH1D *g_Ncoll;
        TH1D *g_Npart;

        TTree *gTree;

        Bool_t makehists;
        Bool_t maketrees;

    public:
        GlauberEvent();
        GlauberEvent(const int A, double_t a, double_t NN_crossSection, Option_t *option);
        ~GlauberEvent();


        void CollisionReset();
        void Collide();
        void FillNucleus(std::vector<double_t> &nucleusX, std::vector<double_t> &nucleusY, Option_t *optionb);

        double_t GetNucleonDistance(int nucAiterator, int nucBiterator, vector<double_t> &xNucleusA, vector<double_t> &xNucleusB, vector<double_t> &yNucleusA, vector<double_t> &yNucleusB);
        double_t GetImpactParameter();
        Int_t GetNcollisions();
        Int_t GetNparticipants();
        TH1D* GetbHist(){return g_bHist;};
        TH1D* GetNpartHist(){return g_Npart;};
        TH1D* GetNcollHist(){return g_Ncoll;};
        TTree* GetTree(){return gTree;}
    ClassDef(GlauberEvent,1);
};

#endif