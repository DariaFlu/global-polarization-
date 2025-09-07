#include "set_lambda_parameterization.hpp"


void set_lambda_parameterization(TFile* Lambda_yield, Double_t fBVal, UParticle &ULambda){ //Valeriy's function for properly lambda generation

    TRandom3* rand = new TRandom3(0);  // Seed with 0 for reproducibility

	Double_t centrality = get_centrality(fBVal); //centrality for parameterization (b->centrality for Xe+Xe below)
	Double_t sNN = 2.87; // Energy of the collision in center-of-mass system

	TH2F* h_pt_y;
	if(centrality<10) h_pt_y = (TH2F*) Lambda_yield->Get("h2PartYpT010");
	else if(centrality>10 && centrality <40) h_pt_y = (TH2F*) Lambda_yield->Get("h2PartYpT");
	else h_pt_y = (TH2F*) Lambda_yield->Get("h2PartYpT40100");


    Double_t lambda_pT; // pT of Lambda from pT-y TH2F
	Double_t lambda_y;  // rapidity of Lambda from pT-y TH2F

	h_pt_y->GetRandom2(lambda_y,lambda_pT,rand);

	Double_t v1 = (28.8635/(TMath::Power(sNN,2.89092))) 
                * ((-0.0233*centrality+0.5413* TMath::Power(centrality,1./3) ) 
                * (163.536/18.0188 - 163.536/(lambda_pT +18.0188) )* lambda_y + ( -0.0056*centrality+0.377*TMath::Power(centrality,1./3) ) 
                * (0.6653* lambda_pT - 0.6172 * TMath::Power(lambda_pT,2) +0.1154 * TMath::Power(lambda_pT,3) ) * TMath::Power(lambda_y,3) ); // v1 cent-pT-y func
    // Double_t v2 = 0.07; //uniform distribution
    Double_t v2 = get_V2(sNN, centrality, lambda_pT, lambda_y);
    // std::cout<<"================ "<<v2<<std::endl;
    // if(v2 == 0) v2 =0.07;
    if( v1  > 1 )  v1 =1;
    else if(v1<-1) v1=-1;
    // if( v2  > 1 ) v2 =1;
    // std::cout<<" v2 value =  "<<v2<<std::endl;
    // generate phi according to v1 and v2
    // v1 = 0.0;   //generation WITHOUT flows
    // v2 = 0.0;

    static TF1 f("f", "[0]*(1+2*[1]*TMath::Cos(x)+2*[2]*TMath::Cos(2*x))+[3]", 0,2*TMath::Pi());
    Double_t a1=1+2*v1+2*v2;
    Double_t a2=1-2*v1+2*v2;
    Double_t a3=1-v1*v1/(4*v2)-2*v2;
    Double_t a=0; //
    if (a1<a) a=a1;  // find analytic minimun to shift
    if (a2<a) a=a2;  // find analytic minimun to shift
    if (a3<a) a=a3;  // find analytic minimun to shift
    f.SetParameter(0,1/(2*TMath::Pi()*(1-a))); // norm
    f.SetParameter(1,v1);  // v1
    f.SetParameter(2,v2);  // v2
    f.SetParameter(3,-a/(2*TMath::Pi()*(1-a))); // shift to have probability
    f.SetNpx(10000);  // to get a better result when using TF1::GetRandom
    Double_t phi=f.GetRandom(rand);

    // Double_t fEnergyLambda = ULambda.GetMomentum().E();
    // Calculate total energy (E) properly
    TLorentzVector vec;
    // Calculate transverse mass (m_T)
    Double_t lambda_mass = 1.115683; // GeV/c² (Lambda mass)
    Double_t mT = sqrt(lambda_pT * lambda_pT + lambda_mass * lambda_mass);
    // Convert rapidity (y) to pseudorapidity (η)
    Double_t pz = mT * sinh(lambda_y); // longitudinal momentum
    Double_t lambda_eta = (pz != 0.0) ? TMath::ATanH(pz / sqrt(pz*pz + lambda_pT*lambda_pT)) : 0.0;
    Double_t fEnergyLambda = sqrt(lambda_pT*lambda_pT * cosh(lambda_eta)*cosh(lambda_eta) + lambda_mass*lambda_mass);
    vec.SetPtEtaPhiE(lambda_pT, lambda_eta, phi, fEnergyLambda);
    ULambda.SetMomentum(vec); //ULambda

    delete rand;


}
