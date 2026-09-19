{

	TString vol = "SensPartPV";
	TString iso = "1a1z";

	TChain *c = new TChain("clusters_tree");
	TChain *f = new TChain("finfo");
	
	c->Add("../test_data/"+vol+"/"+vol+iso+".root");
	f->Add("../test_data/"+vol+"/"+vol+iso+".root");
	
	c->SetAlias("ene","charge_total.fVal");
	c->SetAlias("rmsx","charge_rms_x.fVal");
	c->SetAlias("rmsy","charge_rms_y.fVal");
	c->SetAlias("x","charge_mean_x.fVal");
	c->SetAlias("y","charge_mean_y.fVal");
	c->SetAlias("z","simz_mean.fVal");
	c->SetAlias("np","size_npixels.fVal");
	
	c->SetAlias("erange","ene<26 && ene>0.004");

}


