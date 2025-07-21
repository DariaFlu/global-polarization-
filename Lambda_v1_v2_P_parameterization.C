void Lambda_v1_v2_P_parameterization(){

	double B; // impact parameter
	double centrality; //centrality for parameterization (b->centrality for Xe+Xe below)
	double sNN = 2.87; // Energy of the collision in center-of-mass system
	if(B<3.44)
                centrality=5;
        else if(B<4.88)
                centrality=15;
        else if(B<5.84)
                centrality=25;
        else if(B<6.64)
                centrality=35;
        else if(B<7.44)
                centrality=45;
        else if(B<8.08)
                centrality=55;
        else if(B<8.72)
                centrality=65;
        else if(B<9.36)
                centrality=75;
        else if(B<9.84)
                centrality=85;
        else
                centrality=95;


	TFile* Lambda_yield = new TFile("/lhep/users/dflusova/lambda/afterburner/v.6/qa_out_xexe.root","read"); // File with TH2F pt-y
	TH2F* h_pt_y;
	if(centrality<10)
		h_pt_y = (TH2F*) Lambda_yield->Get("h2PartYpT010");
	else if(centrality>10 && centrality <40)
		h_pt_y = (TH2F*) Lambda_yield->Get("h2PartYpT");
	else
		h_pt_y = (TH2F*) Lambda_yield->Get("h2PartYpT40100");


	double lambda_pT; // pT of Lambda from pT-y TH2F
	double lambda_y;  // rapidity of Lambda from pT-y TH2F

	h_pt_y->GetRandom2(lambda_y,lambda_pT);


	double v1 = (28.8635/(TMath::Power(sNN,2.89092))) * (  (-0.0233*centrality+0.5413* TMath::Power(centrality,1./3) ) * (163.536/18.0188 - 163.536/(lambda_pT +18.0188) )* lambda_y + ( -0.0056*centrality+0.377*TMath::Power(centrality,1./3) ) * (0.6653* lambda_pT - 0.6172 * TMath::Power(lambda_pT,2) +0.1154 * TMath::Power(lambda_pT,3) ) * TMath::Power(lambda_y,3) ); // v1 cent-pT-y func
	double v2 = 0.07; //uniform distribution

	// generate phi according to v1 and v2
      static TF1 f("f", "[0]*(1+2*[1]*TMath::Cos(x)+2*[2]*TMath::Cos(2*x))+[3]", 0,2*TMath::Pi());
      double a1=1+2*v1+2*v2;
      double a2=1-2*v1+2*v2;
      double a3=1-v1*v1/(4*v2)-2*v2;
      double a=0; //
      if (a1<a) a=a1;  // find analytic minimun to shift
      if (a2<a) a=a2;  // find analytic minimun to shift
      if (a3<a) a=a3;  // find analytic minimun to shift
      f.SetParameter(0,1/(2*TMath::Pi()*(1-a))); // norm
      f.SetParameter(1,v1);  // v1
      f.SetParameter(2,v2);  // v2
      f.SetParameter(3,-a/(2*TMath::Pi()*(1-a))); // shift to have probability
      f.SetNpx(10000);  // to get a better result when using TF1::GetRandom
      double phi=f.GetRandom();


      //Global polarization parameterization

      double mean_P = (2.8569/ TMath::Power(sNN,0.955513) ) * (2.4702 - 0.0461*centrality + 0.0042 * TMath::Power(centrality, 2)); // Value in % 



}
