#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TFile.h"
#include <iostream>

void SV_scan()
{
	const int nfile = 200;

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

	float vtx_rave[100][3];
	float vtx_rave_err[100][3];

	float vtx_rave_refit[100][3];
	float vtx_rave_refit_err[100][3];

	int njet;
	int jet_prop_parton[100];
	int jet_prop_hadron[100];
	float jet_pt[100];
	float jet_mass[100];

	int nsv[100];
	float sv_x[100][30];
	float sv_y[100][30];
	float sv_z[100][30];
	float sv_pt[100][30];
	float sv_mass[100][30];
	float sv_mass_corr[100][30];
	int sv_ntrack[100][30];

	TH1D* hjet_nsv[3];
	TH1D* hjet_sv_dist[3];
	TH1D* hjet_sv_mass[3];
	TH1D* hjet_sv_mass_corr[3];
	TH1D* hjet_sv_ntrack[3];
	TH2D* hjet_pt_mass[3];
	for(int iflav=0; iflav<3; iflav++)
	{
		hjet_nsv[iflav] = new TH1D(Form("hjet_nsv_%d",iflav),"",10,0,10);
		hjet_sv_dist[iflav] = new TH1D(Form("hjet_sv_dist_%d",iflav),"",200,0,2);
		hjet_sv_mass[iflav] = new TH1D(Form("hjet_sv_mass_%d",iflav),"",50,0,10);
		hjet_sv_mass_corr[iflav] = new TH1D(Form("hjet_sv_mass_corr_%d",iflav),"",50,0,10);
		hjet_sv_ntrack[iflav] = new TH1D(Form("hjet_sv_ntrack_%d",iflav),"",10,0,10);
		hjet_pt_mass[iflav] = new TH2D(Form("hjet_pt_mass_%d",iflav),"",900,10,100,100,0,100);
	}//iflav

	for (int ii=0; ii<nfile; ii++){

		//TFile *infile = new TFile("AnaSvtxVertex.root","read");
		TFile *infile = new TFile(Form("jobdir_AnaVertex_40/AnaSvtxVertex_%04d.root",ii),"read");
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
		//T->SetBranchAddress("vtx_reco_ntrack",vtx_reco_ntrack);

		T->SetBranchAddress("nvertex_rave",&nvertex_rave);
		T->SetBranchAddress("vtx_rave",vtx_rave);
		//T->SetBranchAddress("vtx_rave_ntrack",vtx_rave_ntrack);

		T->SetBranchAddress("nvertex_rave_refit",&nvertex_rave_refit);
		T->SetBranchAddress("vtx_rave_refit",vtx_rave_refit);
		//T->SetBranchAddress("vtx_rave_refit_ntrack",vtx_rave_refit_ntrack);

		/*
		T->SetBranchAddress("nvertex_acts",&nvertex_acts);
		T->SetBranchAddress("vtx_acts",vtx_acts);
		T->SetBranchAddress("vtx_acts_ntrack",vtx_acts_ntrack);

		T->SetBranchAddress("nvertex_acts_refit",&nvertex_acts_refit);
		T->SetBranchAddress("vtx_acts_refit",vtx_acts_refit);
		T->SetBranchAddress("vtx_acts_refit_ntrack",vtx_acts_refit_ntrack);
		*/

		T->SetBranchAddress("njet04_true",&njet);
		T->SetBranchAddress("jet04_prop_parton",jet_prop_parton);
		T->SetBranchAddress("jet04_prop_hadron",jet_prop_hadron);
		T->SetBranchAddress("jet04_pt",jet_pt);
		T->SetBranchAddress("jet04_mass",jet_mass);

		T->SetBranchAddress("rave_sv_pT10_nvtx",nsv);
		T->SetBranchAddress("rave_sv_pT10_vtx_x",sv_x);
		T->SetBranchAddress("rave_sv_pT10_vtx_y",sv_y);
		T->SetBranchAddress("rave_sv_pT10_vtx_z",sv_z);
		T->SetBranchAddress("rave_sv_pT10_vtx_pt",sv_pt);
		T->SetBranchAddress("rave_sv_pT10_vtx_mass",sv_mass);
		T->SetBranchAddress("rave_sv_pT10_vtx_mass_corr",sv_mass_corr);
		T->SetBranchAddress("rave_sv_pT10_vtx_ntrack",sv_ntrack);

		int nentries = T->GetEntries();
		nevents += nentries;

		for(int ien=0; ien<nentries; ien++)
		{

			T->GetEntry(ien);

			for(int ijet=0;ijet<njet;ijet++)
			{
				int idx_flv;

				//cout << "jet flavor : " << jet_prop_parton[ijet] << endl;

				/*
				if(abs(jet_prop_hadron[ijet])==0) idx_flv=0;
				else if(abs(jet_prop_hadron[ijet])==4) idx_flv=1;
				else if(abs(jet_prop_hadron[ijet])==5) idx_flv=2;
				*/
				if(abs(jet_prop_parton[ijet])==0) idx_flv=0;
				else if(abs(jet_prop_parton[ijet])==4) idx_flv=1;
				else if(abs(jet_prop_parton[ijet])==5) idx_flv=2;
				else continue;
				
				//cout << "idx_flv :" << idx_flv << endl;

				//cout << "jet pt : " << jet_pt[ijet] << "	jet mass : " << jet_mass[ijet] << endl;
				hjet_pt_mass[idx_flv]->Fill(jet_pt[ijet],jet_mass[ijet]);
				hjet_nsv[idx_flv]->Fill(nsv[ijet]);

				if(nsv[ijet]==0) continue;
				float max_dist = 0;
				float far_sv_mass;
				float far_sv_mass_corr;
				int far_sv_ntrack;
				for(int isv=0;isv<nsv[ijet];isv++)
				{
					//hjet_sv_mass[idx_flv]->Fill(sv_mass[ijet][isv]);

					float dx = sv_x[ijet][isv]-vtx_rave_refit[0][0];
					float dy = sv_y[ijet][isv]-vtx_rave_refit[0][1];
					float dz = sv_z[ijet][isv]-vtx_rave_refit[0][2];

					float dist = sqrt(dx*dx+dy*dy+dz*dz);

					if( dist > 2.0 ) continue;

					if(dist>max_dist)
					{
						max_dist = dist;
						far_sv_mass = sv_mass[ijet][isv];
						far_sv_mass_corr = sv_mass_corr[ijet][isv];
						far_sv_ntrack = sv_ntrack[ijet][isv];
					}

				}//sv loop

				hjet_sv_dist[idx_flv]->Fill(max_dist);
				//
				//if( max_dist > 0.1 && max_dist < 2.0 ) 
				if( max_dist > 0.1 )
				{
					hjet_sv_mass[idx_flv]->Fill(far_sv_mass);
					hjet_sv_mass_corr[idx_flv]->Fill(far_sv_mass_corr);
					hjet_sv_ntrack[idx_flv]->Fill(far_sv_ntrack);
				}//SV flight distance cut

			}//jet loop

		}//ien

		infile->Close();
		delete infile;

	}//ii (file loop)

	TFile* outfile = new TFile("outfile_SV_scan.root","recreate");
	for(int iflav=0; iflav<3; iflav++){
		hjet_nsv[iflav]->Write();
		hjet_sv_dist[iflav]->Write();
		hjet_pt_mass[iflav]->Write();
		hjet_sv_mass[iflav]->Write();
		hjet_sv_mass_corr[iflav]->Write();
		hjet_sv_ntrack[iflav]->Write();
	}//iflav

	outfile->Close();

}//SV_scan()
