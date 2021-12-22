void ScanAnaVertex(){

	const int nfile = 800;

	int nevents=0;

	float vtx_gen[3];

	int nvertex_reco;
	int nvertex_rave;
	int nvertex_rave_refit;
	int nvertex_acts;
	int nvertex_acts_refit;

	int ntrack_reco;
	int ntrack_reco_mvtx;
	int tot_ntrack_reco=0;

	float vtx_reco[100][3];
	float vtx_reco_err[100][3];
	int vtx_reco_ntrack[100];

	float vtx_rave[100][3];
	float vtx_rave_err[100][3];
	int vtx_rave_ntrack[100];

	float vtx_rave_refit[100][3];
	float vtx_rave_refit_err[100][3];
	int vtx_rave_refit_ntrack[100];

	float vtx_acts[100][3];
	float vtx_acts_err[100][3];
	int vtx_acts_ntrack[100];

	float vtx_acts_refit[100][3];
	float vtx_acts_refit_err[100][3];
	int vtx_acts_refit_ntrack[100];

	int njet04_true;
	float jet04_pt[100];
	float jet04_eta[100];
	float jet04_mass[100];

	int rave_sv_pT10_nvtx[100];
	float rave_sv_pT10_vtx_x[100][30];
	float rave_sv_pT10_vtx_y[100][30];
	float rave_sv_pT10_vtx_z[100][30];

	TH1D *hnvertex_reco = new TH1D("hnvertex_reco","",101,-0.5,100.5);
	TH1D *hnvertex_rave = new TH1D("hnvertex_rave","",101,-0.5,100.5);
	TH1D *hnvertex_rave_refit = new TH1D("hnvertex_rave_refit","",101,-0.5,100.5);
	TH1D *hnvertex_acts = new TH1D("hnvertex_acts","",101,-0.5,100.5);
	TH1D *hnvertex_acts_refit = new TH1D("hnvertex_acts_refit","",101,-0.5,100.5);

	TH2D *hntrk_pv_rave_acts = new TH2D("hntrk_pv_rave_acts","",51,-0.5,50.5,51,-0.5,50.5);

	TH2D *hvtx_res_acts[3];
	TH2D *hvtx_res_acts_refit[3];
	TH2D *hvtx_res_rave[3];
	TH2D *hvtx_res_rave_refit[3];

#if 0//total efficiency
	TH2D *htrack_eta_pt_mvtx0 = new TH2D("htrack_eta_pt_mvtx0","",30,-1.5,1.5,100,0,10);
	TH2D *htrack_eta_pt_mvtx1 = new TH2D("htrack_eta_pt_mvtx1","",30,-1.5,1.5,100,0,10);
	TH2D *htrack_eta_pt_mvtx2 = new TH2D("htrack_eta_pt_mvtx2","",30,-1.5,1.5,100,0,10);
	TH2D *htrack_eta_pt_mvtx3 = new TH2D("htrack_eta_pt_mvtx3","",30,-1.5,1.5,100,0,10);
#endif//total efficiency

#if 1//prompt vs. non-prompt
	TH2D *htrack_eta_pt_p_mvtx0 = new TH2D("htrack_eta_pt_p_mvtx0","",30,-1.5,1.5,100,0,10);
	TH2D *htrack_eta_pt_p_mvtx1 = new TH2D("htrack_eta_pt_p_mvtx1","",30,-1.5,1.5,100,0,10);
	TH2D *htrack_eta_pt_p_mvtx2 = new TH2D("htrack_eta_pt_p_mvtx2","",30,-1.5,1.5,100,0,10);
	TH2D *htrack_eta_pt_p_mvtx3 = new TH2D("htrack_eta_pt_p_mvtx3","",30,-1.5,1.5,100,0,10);

	TH2D *htrack_eta_pt_s_mvtx0 = new TH2D("htrack_eta_pt_s_mvtx0","",30,-1.5,1.5,100,0,10);
	TH2D *htrack_eta_pt_s_mvtx1 = new TH2D("htrack_eta_pt_s_mvtx1","",30,-1.5,1.5,100,0,10);
	TH2D *htrack_eta_pt_s_mvtx2 = new TH2D("htrack_eta_pt_s_mvtx2","",30,-1.5,1.5,100,0,10);
	TH2D *htrack_eta_pt_s_mvtx3 = new TH2D("htrack_eta_pt_s_mvtx3","",30,-1.5,1.5,100,0,10);
#endif//prompt vs. non-prompt

	TH3D *hjet04_mass_eta_pt = new TH3D("hjet04_mass_eta_pt","",100,0,100,30,-1.5,1.5,900,10,100);
	TH1D *hdistance_pv_sv = new TH1D("hdistance_pv_sv","",500,0,5);

	for (int jj=0; jj<3; jj++){
		hvtx_res_rave[jj] = new TH2D(Form("hvtx_res_rave_%d",jj),"",50,0,50,2000,-0.5,0.5);
		hvtx_res_rave_refit[jj] = new TH2D(Form("hvtx_res_rave_refit_%d",jj),"",50,0,50,2000,-0.5,0.5);
		hvtx_res_acts[jj] = new TH2D(Form("hvtx_res_acts_%d",jj),"",50,0,50,2000,-0.5,0.5);
		hvtx_res_acts_refit[jj] = new TH2D(Form("hvtx_res_acts_refit_%d",jj),"",50,0,50,2000,-0.5,0.5);
	}

	TH2D *hntrk_tpc_mvtx = new TH2D("hntrk_tpc_mvtx","hntrk_tpc_mvtx",51,-0.5,50.5,51,-0.5,50.5);

	for (int ii=0; ii<nfile; ii++){

		//TFile *infile = new TFile("AnaSvtxVertex.root","read");
		TFile *infile = new TFile(Form("jobdir_AnaVertex_38/AnaSvtxVertex_%04d.root",ii),"read");
		//TFile *infile = new TFile(Form("test_jobdir_AnaVertex_08/AnaSvtxVertex_%04d.root",ii),"read");
		if(infile->IsOpen())
		{
			cout << "OPEN: " << infile->GetName() << endl;
		}else
		{
			infile->Close();
			delete infile;
			continue;
		}

		TTree *T = (TTree*)infile->Get("events");

		if(!T)
		{
			infile->Close();
			delete infile;
			continue;
		}

		T->SetBranchAddress("vtx_gen",vtx_gen);

		T->SetBranchAddress("ntrack_reco",&ntrack_reco);
		T->SetBranchAddress("ntrack_reco_mvtx",&ntrack_reco_mvtx);

		T->SetBranchAddress("nvertex_reco",&nvertex_reco);
		T->SetBranchAddress("vtx_reco",vtx_reco);
		T->SetBranchAddress("vtx_reco_ntrack",vtx_reco_ntrack);

		T->SetBranchAddress("nvertex_rave",&nvertex_rave);
		T->SetBranchAddress("vtx_rave",vtx_rave);
		T->SetBranchAddress("vtx_rave_ntrack",vtx_rave_ntrack);

		T->SetBranchAddress("nvertex_rave_refit",&nvertex_rave_refit);
		T->SetBranchAddress("vtx_rave_refit",vtx_rave_refit);
		T->SetBranchAddress("vtx_rave_refit_ntrack",vtx_rave_refit_ntrack);

		T->SetBranchAddress("nvertex_acts",&nvertex_acts);
		T->SetBranchAddress("vtx_acts",vtx_acts);
		T->SetBranchAddress("vtx_acts_ntrack",vtx_acts_ntrack);

		T->SetBranchAddress("nvertex_acts_refit",&nvertex_acts_refit);
		T->SetBranchAddress("vtx_acts_refit",vtx_acts_refit);
		T->SetBranchAddress("vtx_acts_refit_ntrack",vtx_acts_refit_ntrack);

		T->SetBranchAddress("njet04_true",&njet04_true);
		T->SetBranchAddress("jet04_pt",jet04_pt);
		T->SetBranchAddress("jet04_eta",jet04_eta);
		T->SetBranchAddress("jet04_mass",jet04_mass);

		T->SetBranchAddress("rave_sv_pT10_nvtx",rave_sv_pT10_nvtx);
		T->SetBranchAddress("rave_sv_pT10_vtx_x",rave_sv_pT10_vtx_x);
		T->SetBranchAddress("rave_sv_pT10_vtx_y",rave_sv_pT10_vtx_y);
		T->SetBranchAddress("rave_sv_pT10_vtx_z",rave_sv_pT10_vtx_z);

#if 0//total efficiency
		TH2D *htmp0 = (TH2D*) infile->Get("h2d_eta_pt_mvtx0");
		TH2D *htmp1 = (TH2D*) infile->Get("h2d_eta_pt_mvtx1");
		TH2D *htmp2 = (TH2D*) infile->Get("h2d_eta_pt_mvtx2");
		TH2D *htmp3 = (TH2D*) infile->Get("h2d_eta_pt_mvtx3");
#endif//total efficiency

#if 1//prompt vs. non-prompt
		TH2D *htmp_p0 = (TH2D*) infile->Get("h2d_eta_pt_p_mvtx0");
		TH2D *htmp_p1 = (TH2D*) infile->Get("h2d_eta_pt_p_mvtx1");
		TH2D *htmp_p2 = (TH2D*) infile->Get("h2d_eta_pt_p_mvtx2");
		TH2D *htmp_p3 = (TH2D*) infile->Get("h2d_eta_pt_p_mvtx3");

		TH2D *htmp_s0 = (TH2D*) infile->Get("h2d_eta_pt_s_mvtx0");
		TH2D *htmp_s1 = (TH2D*) infile->Get("h2d_eta_pt_s_mvtx1");
		TH2D *htmp_s2 = (TH2D*) infile->Get("h2d_eta_pt_s_mvtx2");
		TH2D *htmp_s3 = (TH2D*) infile->Get("h2d_eta_pt_s_mvtx3");
#endif//prompt vs. non-prompt

#if 0//total efficiency
		htrack_eta_pt_mvtx0->Add(htmp0);
		htrack_eta_pt_mvtx1->Add(htmp1);
		htrack_eta_pt_mvtx2->Add(htmp2);
		htrack_eta_pt_mvtx3->Add(htmp3);
#endif//total efficiency

#if 1//prompt vs. non-prompt
		htrack_eta_pt_p_mvtx0->Add(htmp_p0);
		htrack_eta_pt_p_mvtx1->Add(htmp_p1);
		htrack_eta_pt_p_mvtx2->Add(htmp_p2);
		htrack_eta_pt_p_mvtx3->Add(htmp_p3);

		htrack_eta_pt_s_mvtx0->Add(htmp_s0);
		htrack_eta_pt_s_mvtx1->Add(htmp_s1);
		htrack_eta_pt_s_mvtx2->Add(htmp_s2);
		htrack_eta_pt_s_mvtx3->Add(htmp_s3);
#endif//prompt vs. non-prompt

		int nentries = T->GetEntries();
		nevents += nentries;

		for (int ien=0; ien<nentries; ien++){

			T->GetEntry(ien);
			//tot_ntrack_reco += ntrack_reco;

			//hntrk_tpc_mvtx->Fill(ntrack_reco, ntrack_reco_mvtx);

			hnvertex_reco->Fill(nvertex_reco);
			hnvertex_rave->Fill(nvertex_rave);
			hnvertex_rave_refit->Fill(nvertex_rave_refit);
			hnvertex_acts->Fill(nvertex_acts);
			hnvertex_acts_refit->Fill(nvertex_acts_refit);

#if 0//comparing acts vs rave
			if ( nvertex_rave==1 && nvertex_acts==1 ){
				hntrk_pv_rave_acts->Fill(vtx_rave_ntrack[0], vtx_acts_ntrack[0]);
			}
#endif //comparing acts vs rave

			if ( nvertex_rave==1 ){
				for (int jj=0; jj<3; jj++){
					hvtx_res_rave[jj]->Fill(vtx_rave_ntrack[0], vtx_rave[0][jj] - vtx_gen[jj]);
				}//jj
			}

			if ( nvertex_rave_refit==1 ){
				for (int jj=0; jj<3; jj++){
					hvtx_res_rave_refit[jj]->Fill(vtx_rave_refit_ntrack[0], vtx_rave_refit[0][jj] - vtx_gen[jj]);
				}//jj
			}

			if ( nvertex_acts==1 && vtx_acts_ntrack[0]>0 ){
				for (int jj=0; jj<3; jj++){
					hvtx_res_acts[jj]->Fill(vtx_acts_ntrack[0], vtx_acts[0][jj] - vtx_gen[jj]);
				}//jj
			}

			if ( nvertex_acts_refit==1 && vtx_acts_refit_ntrack[0]>0 ){
				for (int jj=0; jj<3; jj++){
					hvtx_res_acts_refit[jj]->Fill(vtx_acts_refit_ntrack[0], vtx_acts_refit[0][jj] - vtx_gen[jj]);
				}//jj
			}

			for (int ijet=0; ijet<njet04_true;ijet++)
			{
				float pt = jet04_pt[ijet];
				float eta = jet04_eta[ijet];
				float mass = jet04_mass[ijet];

				hjet04_mass_eta_pt->Fill(mass,eta,pt);

				for (int isv=0;isv<rave_sv_pT10_nvtx[ijet];isv++)
				{
					float sv_x = rave_sv_pT10_vtx_x[ijet][isv];
					float sv_y = rave_sv_pT10_vtx_y[ijet][isv];
					float sv_z = rave_sv_pT10_vtx_z[ijet][isv];

					float pv_x = vtx_rave[0][0];
					float pv_y = vtx_rave[0][1];
					float pv_z = vtx_rave[0][2];

					float dist = sqrt(pow(pv_x-sv_x,2)+pow(pv_y-sv_y,2)+pow(pv_z-sv_z,2));
					hdistance_pv_sv->Fill(dist);
				}//-sv loop
			}//ijet

		}//ien
		infile->Close();
		delete infile;

	}//ii (file loop)

#if 1//not making outfile
	TFile *outfile = new TFile("outfile_AnaVertex.root","recreate");
	hnvertex_reco->Write();
	hnvertex_rave->Write();
	hnvertex_rave_refit->Write();
	hnvertex_acts->Write();
	hnvertex_acts_refit->Write();

#if 0//total efficiency
	htrack_eta_pt_mvtx0->Write();
	htrack_eta_pt_mvtx1->Write();
	htrack_eta_pt_mvtx2->Write();
	htrack_eta_pt_mvtx3->Write();
#endif//total efficiency

#if 1//prompt vs. non-prompt
	htrack_eta_pt_p_mvtx0->Write();
	htrack_eta_pt_p_mvtx1->Write();
	htrack_eta_pt_p_mvtx2->Write();
	htrack_eta_pt_p_mvtx3->Write();

	htrack_eta_pt_s_mvtx0->Write();
	htrack_eta_pt_s_mvtx1->Write();
	htrack_eta_pt_s_mvtx2->Write();
	htrack_eta_pt_s_mvtx3->Write();
#endif//prompt vs. non-prompt

#if 0//comparing acts vs rave
	hntrk_pv_rave_acts->Write();
#endif //comparing acts vs rave

	hntrk_tpc_mvtx->Write();

	for (int jj=0; jj<3; jj++){
		hvtx_res_rave[jj]->Write();
		hvtx_res_rave_refit[jj]->Write();
		hvtx_res_acts[jj]->Write();
		hvtx_res_acts_refit[jj]->Write();
	}

	hdistance_pv_sv->Write();
	hjet04_mass_eta_pt->Write();

	outfile->Close();
#endif //not making outfile

}
