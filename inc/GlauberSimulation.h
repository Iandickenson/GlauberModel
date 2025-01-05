
#ifndef MULTIPLICITYHIST_H
#define MULTIPLICITYHIST_H

class TH1D;
class TF1;

class GlauberSimulation : public TObject{
    private:
        double_t gMu;
        double_t gK;

        TH1D *g_bHist;
        TH1D *g_Ncoll;
        TH1D *g_Npart;

        TH1D *gMultiplicityHist;

        TF1 *gNegativeBinomial;

        // Auxilliary variables
        double_t binContent;
        double_t maxNcoll;
        double_t thisMultiplicity;
        Int_t collisions;
    public:
        GlauberSimulation();
        ~GlauberSimulation();
        GlauberSimulation(double_t mu, double_t k, TH1D *ncoll, TH1D *npart, TH1D *bhist);
        GlauberSimulation(double_t mu, double_t k, GlauberEvent *g);

        void GenerateMultiplicity();
        void GenerateMultiplicity(double_t mu, double_t k);
        void SetNegativeBinomialParameters(double_t mu, double_t k);

        TH1D* GetMultiplicityHist(){return gMultiplicityHist;};

        double_t ChiSquared(TH1D *dataHist);

    ClassDef(GlauberSimulation,1);
};

#endif