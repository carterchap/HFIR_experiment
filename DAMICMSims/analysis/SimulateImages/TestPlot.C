{
  gROOT->Reset();
  gStyle->SetOptStat(0);
  TFile f("root/out_testCopper1.root");
  THnSparseF *img_adcu = (THnSparseF*)f.Get("img_adcu_103");
  double off = 10;
  double xmin = 0;
  double ymin = 0;
  double ymax = 4000;
  double xmax = 4000;

  img_adcu->GetAxis(0)->SetRange(xmin, xmax);
  img_adcu->GetAxis(1)->SetRange(ymin, ymax);
  TH2D *hadc0 = img_adcu->Projection(1,0);
  hadc0->SetTitle(0);

  hadc0->GetYaxis()->SetTitleOffset(1.45);
  TCanvas *c_adc0 = new TCanvas("c_adc0","",2);
  c_adc0->SetLogz();
  hadc0->Draw("colz");


}
