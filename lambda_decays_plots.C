#ifdef __CLING__
#pragma cling optimize(0)
#endif
void lambda_decays_plots()
{
//=========Macro generated from canvas: c1/Lambda_Decays
//=========  (Mon Mar 30 13:18:21 2026) by ROOT version 6.32.06
   TCanvas *c1 = new TCanvas("c1", "Lambda_Decays",0,0,1200,800);
   c1->SetHighLightColor(2);
   c1->Range(0,0,1,1);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: c1_1
   TPad *c1_1__0 = new TPad("c1_1", "c1_1",0.01,0.51,0.49,0.99);
   c1_1__0->Draw();
   c1_1__0->cd();
   c1_1__0->Range(-0.25,-0.13125,2.25,1.18125);
   c1_1__0->SetFillColor(0);
   c1_1__0->SetBorderMode(0);
   c1_1__0->SetBorderSize(2);
   c1_1__0->SetFrameBorderMode(0);
   c1_1__0->SetFrameBorderMode(0);
   
   TH1D *hLambdaPt__1 = new TH1D("hLambdaPt__1","Lambda pT",100,0,2);
   
   TPaveStats *ptstats = new TPaveStats(0.78,0.775,0.98,0.935,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetTextAlign(12);
   ptstats->SetTextFont(42);
   TText *ptstats_LaTex = ptstats->AddText("hLambdaPt");
   ptstats_LaTex->SetTextSize(0.0368);
   ptstats_LaTex = ptstats->AddText("Entries = 0      ");
   ptstats_LaTex = ptstats->AddText("Mean  =      0");
   ptstats_LaTex = ptstats->AddText("Std Dev   =      0");
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(0);
   ptstats->Draw();
   hLambdaPt__1->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(hLambdaPt__1);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   hLambdaPt__1->SetLineColor(ci);
   hLambdaPt__1->GetXaxis()->SetTitle("pT [GeV/c]");
   hLambdaPt__1->GetXaxis()->SetLabelFont(42);
   hLambdaPt__1->GetXaxis()->SetTitleOffset(1);
   hLambdaPt__1->GetXaxis()->SetTitleFont(42);
   hLambdaPt__1->GetYaxis()->SetTitle("Counts");
   hLambdaPt__1->GetYaxis()->SetLabelFont(42);
   hLambdaPt__1->GetYaxis()->SetTitleFont(42);
   hLambdaPt__1->GetZaxis()->SetLabelFont(42);
   hLambdaPt__1->GetZaxis()->SetTitleOffset(1);
   hLambdaPt__1->GetZaxis()->SetTitleFont(42);
   hLambdaPt__1->Draw("");
   
   TPaveText *pt = new TPaveText(0.4011622,0.9334715,0.5988378,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *pt_LaTex = pt->AddText("Lambda pT");
   pt->Draw();
   c1_1__0->Modified();
   c1->cd();
  
// ------------>Primitives in pad: c1_2
   TPad *c1_2__1 = new TPad("c1_2", "c1_2",0.51,0.51,0.99,0.99);
   c1_2__1->Draw();
   c1_2__1->cd();
   c1_2__1->Range(-1.25,-0.13125,1.25,1.18125);
   c1_2__1->SetFillColor(0);
   c1_2__1->SetBorderMode(0);
   c1_2__1->SetBorderSize(2);
   c1_2__1->SetFrameBorderMode(0);
   c1_2__1->SetFrameBorderMode(0);
   
   TH1D *hCosTheta__2 = new TH1D("hCosTheta__2","Proton cos(#theta)",100,-1,1);
   
   ptstats = new TPaveStats(0.78,0.775,0.98,0.935,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetTextAlign(12);
   ptstats->SetTextFont(42);
   ptstats_LaTex = ptstats->AddText("hCosTheta");
   ptstats_LaTex->SetTextSize(0.0368);
   ptstats_LaTex = ptstats->AddText("Entries = 0      ");
   ptstats_LaTex = ptstats->AddText("Mean  =      0");
   ptstats_LaTex = ptstats->AddText("Std Dev   =      0");
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(0);
   ptstats->Draw();
   hCosTheta__2->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(hCosTheta__2);

   ci = TColor::GetColor("#000099");
   hCosTheta__2->SetLineColor(ci);
   hCosTheta__2->GetXaxis()->SetTitle("cos(#theta)");
   hCosTheta__2->GetXaxis()->SetLabelFont(42);
   hCosTheta__2->GetXaxis()->SetTitleOffset(1);
   hCosTheta__2->GetXaxis()->SetTitleFont(42);
   hCosTheta__2->GetYaxis()->SetTitle("Counts");
   hCosTheta__2->GetYaxis()->SetLabelFont(42);
   hCosTheta__2->GetYaxis()->SetTitleFont(42);
   hCosTheta__2->GetZaxis()->SetLabelFont(42);
   hCosTheta__2->GetZaxis()->SetTitleOffset(1);
   hCosTheta__2->GetZaxis()->SetTitleFont(42);
   hCosTheta__2->Draw("");
   
   pt = new TPaveText(0.384614,0.9367098,0.615386,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   pt_LaTex = pt->AddText("Proton cos(#theta)");
   pt->Draw();
   c1_2__1->Modified();
   c1->cd();
  
// ------------>Primitives in pad: c1_3
   TPad *c1_3__2 = new TPad("c1_3", "c1_3",0.01,0.01,0.49,0.49);
   c1_3__2->Draw();
   c1_3__2->cd();
   c1_3__2->Range(-0.25,-1.25,2.25,1.25);
   c1_3__2->SetFillColor(0);
   c1_3__2->SetBorderMode(0);
   c1_3__2->SetBorderSize(2);
   c1_3__2->SetFrameBorderMode(0);
   c1_3__2->SetFrameBorderMode(0);
   
   TH2D *hAngVsPt = new TH2D("hAngVsPt","cos(#theta) vs Lambda pT",50,0,2,50,-1,1);
   
   ptstats = new TPaveStats(0.78,0.695,0.98,0.935,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetTextAlign(12);
   ptstats->SetTextFont(42);
   ptstats_LaTex = ptstats->AddText("hAngVsPt");
   ptstats_LaTex->SetTextSize(0.0368);
   ptstats_LaTex = ptstats->AddText("Entries = 0      ");
   ptstats_LaTex = ptstats->AddText("Mean x =      0");
   ptstats_LaTex = ptstats->AddText("Mean y =      0");
   ptstats_LaTex = ptstats->AddText("Std Dev x =      0");
   ptstats_LaTex = ptstats->AddText("Std Dev y =      0");
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(0);
   ptstats->Draw();
   hAngVsPt->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(hAngVsPt);

   ci = TColor::GetColor("#000099");
   hAngVsPt->SetLineColor(ci);
   hAngVsPt->GetXaxis()->SetTitle("pT [GeV/c]");
   hAngVsPt->GetXaxis()->SetLabelFont(42);
   hAngVsPt->GetXaxis()->SetTitleOffset(1);
   hAngVsPt->GetXaxis()->SetTitleFont(42);
   hAngVsPt->GetYaxis()->SetTitle("cos(#theta)");
   hAngVsPt->GetYaxis()->SetLabelFont(42);
   hAngVsPt->GetYaxis()->SetTitleFont(42);
   hAngVsPt->GetZaxis()->SetLabelFont(42);
   hAngVsPt->GetZaxis()->SetTitleOffset(1);
   hAngVsPt->GetZaxis()->SetTitleFont(42);
   hAngVsPt->Draw("colz");
   
   pt = new TPaveText(0.3253888,0.9334715,0.6746112,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   pt_LaTex = pt->AddText("cos(#theta) vs Lambda pT");
   pt->Draw();
   c1_3__2->Modified();
   c1->cd();
  
// ------------>Primitives in pad: c1_4
   TPad *c1_4__3 = new TPad("c1_4", "c1_4",0.51,0.01,0.99,0.49);
   c1_4__3->Draw();
   c1_4__3->cd();
   c1_4__3->Range(-0.7853982,-0.13125,7.068584,1.18125);
   c1_4__3->SetFillColor(0);
   c1_4__3->SetBorderMode(0);
   c1_4__3->SetBorderSize(2);
   c1_4__3->SetFrameBorderMode(0);
   c1_4__3->SetFrameBorderMode(0);
   
   TH1D *hLambdaPhi__3 = new TH1D("hLambdaPhi__3","Lambda (#Delta#varphi)",100,0,6.283185);
   
   ptstats = new TPaveStats(0.78,0.775,0.98,0.935,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetTextAlign(12);
   ptstats->SetTextFont(42);
   ptstats_LaTex = ptstats->AddText("hLambdaPhi");
   ptstats_LaTex->SetTextSize(0.0368);
   ptstats_LaTex = ptstats->AddText("Entries = 0      ");
   ptstats_LaTex = ptstats->AddText("Mean  =      0");
   ptstats_LaTex = ptstats->AddText("Std Dev   =      0");
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(0);
   ptstats->Draw();
   hLambdaPhi__3->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(hLambdaPhi__3);

   ci = TColor::GetColor("#000099");
   hLambdaPhi__3->SetLineColor(ci);
   hLambdaPhi__3->GetXaxis()->SetTitle("#Delta#varphi");
   hLambdaPhi__3->GetXaxis()->SetLabelFont(42);
   hLambdaPhi__3->GetXaxis()->SetTitleOffset(1);
   hLambdaPhi__3->GetXaxis()->SetTitleFont(42);
   hLambdaPhi__3->GetYaxis()->SetTitle("Counts");
   hLambdaPhi__3->GetYaxis()->SetLabelFont(42);
   hLambdaPhi__3->GetYaxis()->SetTitleFont(42);
   hLambdaPhi__3->GetZaxis()->SetLabelFont(42);
   hLambdaPhi__3->GetZaxis()->SetTitleOffset(1);
   hLambdaPhi__3->GetZaxis()->SetTitleFont(42);
   hLambdaPhi__3->Draw("");
   
   pt = new TPaveText(0.3907107,0.9334715,0.6092893,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   pt_LaTex = pt->AddText("Lambda (#Delta#varphi)");
   pt->Draw();
   c1_4__3->Modified();
   c1->cd();
   c1->Modified();
   c1->SetSelected(c1);
}
