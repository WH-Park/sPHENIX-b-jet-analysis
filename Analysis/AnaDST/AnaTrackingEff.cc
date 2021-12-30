////////////////////////////////////////////////////////////////////////////////
//
// This module is desgined to grab svtx tracks and put truth and cluster
// information into a TTree for GenFit testing
//


#include "AnaTrackingEff.h"

#include <phool/phool.h>
#include <phool/getClass.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4VtxPoint.h>
#include <fun4all/PHTFileServer.h>
#include <fun4all/Fun4AllServer.h>

#include <g4hough/SvtxVertexMap.h>
#include <g4hough/SvtxVertex.h>
#include <g4hough/SvtxTrackMap.h>
#include <g4hough/SvtxTrack.h>
#include <g4hough/SvtxClusterMap.h>
#include <g4hough/SvtxCluster.h>
#include <g4hough/SvtxHitMap.h>
#include <g4hough/SvtxHit.h>

#include <g4eval/SvtxEvalStack.h>
#include <g4eval/SvtxTrackEval.h>
#include <g4eval/SvtxClusterEval.h>
#include <g4eval/SvtxTruthEval.h>
#include <g4eval/SvtxVertexEval.h>
#include <g4eval/SvtxHitEval.h>

#include <g4jets/JetMap.h>
#include <g4jets/Jet.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TTree.h>
#include <TVector3.h>
#include <TMath.h>
#include <phhepmc/PHHepMCGenEventMap.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#include <HepMC/GenParticle.h>

#include <HFJetTruthGeneration/HFJetDefs.h>

#include <iostream>
#include <set>

using namespace std;

//----------------------------------------------------------------------------//
//-- Constructor:
//--  simple initialization
//----------------------------------------------------------------------------//
AnaTrackingEff::AnaTrackingEff(const string &name):
  SubsysReco( name ),
  _events( NULL )
{
  //initialize
  _event = 0;
  _outfile = "AnaTrackingEff.root";

	hz_in = NULL;
	hz_out = NULL;
	hz_diff = NULL;

	_truth_container = NULL;
	_hepmc_eventmap = NULL;

	_trkmap = NULL;
	_clustermap = NULL;
	_svtxevalstack = NULL;

	
	_n_maps_layers = 3;
	_n_intt_layers = 4;

}

//----------------------------------------------------------------------------//
//-- Init():
//--   Intialize all histograms, trees, and ntuples
//----------------------------------------------------------------------------//
int AnaTrackingEff::Init(PHCompositeNode *topNode)
{
  cout << PHWHERE << " Openning file " << _outfile << endl;
  PHTFileServer::get().open( _outfile, "RECREATE");

  // create TTree
  _events = new TTree("events", "Svtx Event");

  _events->Branch("event",0,"event/I");

	_events->Branch("vtx_gen",0,"vtx_gen[3]/F");
	_events->Branch("nvertex",0,"nvertex/I");
	_events->Branch("vtx_reco",0,"vtx_reco[nvertex][3]/F");
	//_events->Branch("vtx_reco_err",0,"vtx_reco_err[3]/F");
	_events->Branch("vtx_id",0,"vtx_id[nvertex]/S");

	_events->Branch("npart",0,"npart/I");
	_events->Branch("gen_pid",0,"gen_pid[npart]/S");
	_events->Branch("gen_embed",0,"gen_embed[npart]/S");
	_events->Branch("gen_primary",0,"gen_primary[npart]/S");
	_events->Branch("gen_px",0,"gen_px[npart]/F");
	_events->Branch("gen_py",0,"gen_py[npart]/F");
	_events->Branch("gen_pz",0,"gen_pz[npart]/F");
	_events->Branch("gen_vx",0,"gen_vx[npart]/F");
	_events->Branch("gen_vy",0,"gen_vy[npart]/F");
	_events->Branch("gen_vz",0,"gen_vz[npart]/F");
	_events->Branch("gen_ngmaps",0,"gen_ngmaps[npart]/S");
	_events->Branch("gen_ngintt",0,"gen_ngintt[npart]/S");
	_events->Branch("gen_ngtpc",0,"gen_ngtpc[npart]/S");
	_events->Branch("gen_ngmaps_layer",0,"gen_ngmaps_layer[npart]/S");
	_events->Branch("gen_ngintt_layer",0,"gen_ngintt_layer[npart]/S");
	_events->Branch("gen_ngtpc_layer",0,"gen_ngtpc_layer[npart]/S");

	//_events->Branch("ncluster",0,"ncluster/I");
	//_events->Branch("cluster_layer",0,"cluster_layer[ncluster]/S");
	//_events->Branch("cluster_gprimary",0,"cluster_gprimary[ncluster]/S");
	//_events->Branch("cluster_gembed",0,"cluster_gembed[ncluster]/S");
	//_events->Branch("cluster_gflavor",0,"cluster_gflavor[ncluster]/S");
	//_events->Branch("cluster_gpx",0,"cluster_gpx[ncluster]/F");
	//_events->Branch("cluster_gpy",0,"cluster_gpy[ncluster]/F");
	//_events->Branch("cluster_gvx",0,"cluster_gvx[ncluster]/F");
	//_events->Branch("cluster_gvy",0,"cluster_gvy[ncluster]/F");
	//_events->Branch("cluster_err",0,"cluster_err[ncluster][3]/F");

	_events->Branch("ntrack",0,"ntrack/I");
	_events->Branch("track_pt",0,"track_pt[ntrack]/F");
	_events->Branch("track_eta",0,"track_eta[ntrack]/F");
	_events->Branch("track_chiq",0,"track_chiq[ntrack]/F");
	_events->Branch("track_ndf",0,"track_ndf[ntrack]/S");
	_events->Branch("track_cluster_n_maps",0,"track_cluster_n_maps[ntrack]/S");
	_events->Branch("track_cluster_n_intt",0,"track_cluster_n_intt[ntrack]/S");
	_events->Branch("track_cluster_n_tpc",0,"track_cluster_n_tpc[ntrack]/S");
	_events->Branch("track_cluster_n_layer_maps",0,"track_cluster_n_layer_maps[ntrack]/S");
	_events->Branch("track_cluster_n_layer_intt",0,"track_cluster_n_layer_intt[ntrack]/S");
	_events->Branch("track_cluster_n_layer_tpc",0,"track_cluster_n_layer_tpc[ntrack]/S");
	_events->Branch("track_cluster_n_maps_truth",0,"track_cluster_n_maps_truth[ntrack]/S");
	_events->Branch("track_cluster_n_intt_truth",0,"track_cluster_n_intt_truth[ntrack]/S");
	_events->Branch("track_cluster_n_tpc_truth",0,"track_cluster_n_tpc_truth[ntrack]/S");
	//_events->Branch("track_cluster_maps_pattern",0,"track_cluster_maps_pattern[ntrack]/S");
	//_events->Branch("track_cluster_intt_pattern",0,"track_cluster_intt_pattern[ntrack]/S");
	_events->Branch("track_gprimary",0,"track_gprimary[ntrack]/S");
	_events->Branch("track_gembed",0,"track_gembed[ntrack]/S");
	_events->Branch("track_gflavor",0,"track_gflavor[ntrack]/S");
	_events->Branch("track_parent_id",0,"track_parent_id[ntrack]/S");
	//_events->Branch("track_matched_cluster_n",0,"track_matched_cluster_n[ntrack]/S");
	_events->Branch("track_gvx",0,"track_gvx[ntrack]/F");
	_events->Branch("track_gvy",0,"track_gvy[ntrack]/F");
	_events->Branch("track_gvz",0,"track_gvz[ntrack]/F");
	_events->Branch("track_gpx",0,"track_gpx[ntrack]/F");
	_events->Branch("track_gpy",0,"track_gpy[ntrack]/F");
	_events->Branch("track_gpz",0,"track_gpz[ntrack]/F");
	_events->Branch("track_dca2d",0,"track_dca2d[ntrack]/F");
	_events->Branch("track_dca3d_xy",0,"track_dca3d_xy[ntrack]/F");
	_events->Branch("track_dca3d_z",0,"track_dca3d_z[ntrack]/F");

	hz_in = new TH1D("hz_in","",3000,-150,150);
	hz_out = new TH1D("hz_out","",3000,-150,150);
	hz_diff = new TH1D("hz_diff","",3000,-150,150);
	hn_out = new TH1D("hn_out","",1001,-0.5,1000.5);
	hpT_eta_in = new TH2D("hpT_eta_in","",200,0,20,100,-5,5);
	hpT_eta_in_tpc = new TH2D("hpT_eta_in_tpc","",200,0,20,100,-5,5);
	hpT_eta_in_tpc_mvtx = new TH2D("hpT_eta_in_tpc_mvtx","",200,0,20,100,-5,5);
	hpT_eta_in_tpc_mvtx_intt1 = new TH2D("hpT_eta_in_tpc_mvtx_intt1","",200,0,20,100,-5,5);
	hpT_eta_in_tpc_mvtx_intt2 = new TH2D("hpT_eta_in_tpc_mvtx_intt2","",200,0,20,100,-5,5);
	hpT_eta_in_tpc_mvtx_intt3 = new TH2D("hpT_eta_in_tpc_mvtx_intt3","",200,0,20,100,-5,5);
	hpT_eta_in_tpc_mvtx_intt4 = new TH2D("hpT_eta_in_tpc_mvtx_intt4","",200,0,20,100,-5,5);
	hpT_eta_out = new TH2D("hpT_eta_out","",200,0,20,100,-5,5);
	hpT_eta_out_tpc = new TH2D("hpT_eta_out_tpc","",200,0,20,100,-5,5);
	hpT_eta_out_tpc_mvtx = new TH2D("hpT_eta_out_tpc_mvtx","",200,0,20,100,-5,5);
	hpT_eta_out_tpc_mvtx_intt1 = new TH2D("hpT_eta_out_tpc_mvtx_intt1","",200,0,20,100,-5,5);
	hpT_eta_out_tpc_mvtx_intt2 = new TH2D("hpT_eta_out_tpc_mvtx_intt2","",200,0,20,100,-5,5);
	hpT_eta_out_tpc_mvtx_intt3 = new TH2D("hpT_eta_out_tpc_mvtx_intt3","",200,0,20,100,-5,5);
	hpT_eta_out_tpc_mvtx_intt4 = new TH2D("hpT_eta_out_tpc_mvtx_intt4","",200,0,20,100,-5,5);

	for (int ii=0; ii<4; ii++){
		hclus_intt_zphi_gen[ii] = new TH2D(Form("hclus_intt_zphi_gen_layer%d",ii),"",300,-30,30,300,-TMath::Pi(),TMath::Pi());
		hclus_intt_zphi_reco[ii] = new TH2D(Form("hclus_intt_zphi_reco_layer%d",ii),"",300,-30,30,300,-TMath::Pi(),TMath::Pi());
	}

	hnclus_pT_tpc = new TH2D("hnclus_pT_tpc","",100,0,10,61,-0.5,60.5);
	hnclus_pT_intt = new TH2D("hnclus_pT_intt","",100,0,10,61,-0.5,60.5);
	hnclus_pT_mvtx = new TH2D("hnclus_pT_mvtx","",100,0,10,61,-0.5,60.5);

	hnclus_pT_intt_tpc35 = new TH2D("hnclus_pT_intt_tpc35","",100,0,10,61,-0.5,60.5);
	hnclus_pT_mvtx_tpc35 = new TH2D("hnclus_pT_mvtx_tpc35","",100,0,10,61,-0.5,60.5);

  return 0;
}

//----------------------------------------------------------------------------//
//-- process_event():
//--   Call user instructions for every event.
//--   This function contains the analysis structure.
//----------------------------------------------------------------------------//
int AnaTrackingEff::process_event(PHCompositeNode *topNode)
{
  _event++;
  //if (1)
  if (_event % 1000 == 0)
    cout << PHWHERE << "Events processed: " << _event << endl;

  GetNodes(topNode);

	//cout << "GetNodes" << endl;

  if (!_svtxevalstack) {
    _svtxevalstack = new SvtxEvalStack(topNode);
    _svtxevalstack->set_strict(false);
    //_svtxevalstack->set_verbosity(verbosity + 1);
  } else {
    _svtxevalstack->next_event(topNode);
  }

	//cout << "SvtxEvalStack" << endl;

  fill_tree(topNode);

  return 0;
}

//----------------------------------------------------------------------------//
//-- End():
//--   End method, wrap everything up
//----------------------------------------------------------------------------//
int AnaTrackingEff::End(PHCompositeNode *topNode)
{

	cout << "-----AnaTrackingEff::End------" << endl;

  PHTFileServer::get().cd( _outfile );
  _events->Write();
	hz_in->Write();
	hz_out->Write();
	hz_diff->Write();
	hn_out->Write();
	hpT_eta_in->Write();
	hpT_eta_in_tpc->Write();
	hpT_eta_in_tpc_mvtx->Write();
	hpT_eta_in_tpc_mvtx_intt1->Write();
	hpT_eta_in_tpc_mvtx_intt2->Write();
	hpT_eta_in_tpc_mvtx_intt3->Write();
	hpT_eta_in_tpc_mvtx_intt4->Write();
	hpT_eta_out->Write();
	hpT_eta_out_tpc->Write();
	hpT_eta_out_tpc_mvtx->Write();
	hpT_eta_out_tpc_mvtx_intt1->Write();
	hpT_eta_out_tpc_mvtx_intt2->Write();
	hpT_eta_out_tpc_mvtx_intt3->Write();
	hpT_eta_out_tpc_mvtx_intt4->Write();
  //PHTFileServer::get().close();
	//
	for (int ii=0; ii<4; ii++){
		hclus_intt_zphi_gen[ii]->Write();
		hclus_intt_zphi_reco[ii]->Write();
	}

	hnclus_pT_tpc->Write();
	hnclus_pT_intt->Write();
	hnclus_pT_mvtx->Write();

	hnclus_pT_intt_tpc35->Write();
	hnclus_pT_mvtx_tpc35->Write();

	//if ( trutheval ) trutheval->Delete();
	//if ( trackeval ) trackeval->Delete();

  return 0;
}


//----------------------------------------------------------------------------//
//-- fill_tree():
//--   Fill the trees with truth, track fit, and cluster information
//----------------------------------------------------------------------------//
void AnaTrackingEff::fill_tree(PHCompositeNode *topNode)
{

  // Make sure to reset all the TTree variables before trying to set them.
  reset_variables();

	//cout << "reset_variables" << endl;

	SvtxClusterEval *clustereval = _svtxevalstack->get_cluster_eval();
	SvtxTrackEval *trackeval = _svtxevalstack->get_track_eval();
	SvtxTruthEval *trutheval = _svtxevalstack->get_truth_eval();
	//SvtxVertexEval *vertexeval = _svtxevalstack->get_vertex_eval();

	//cout << "ready" << endl;

	//int ncluster = 0;
	int ntrack = 0;
	int npart = 0;

	cout << "-EVENT: " << _event << "---------------------------------------------------------------------------------------------------" << endl;

	if ( _hepmc_eventmap ){

		float zvertex_in = -9999;

		for (PHHepMCGenEventMap::ConstIter iter=_hepmc_eventmap->begin(); iter!=_hepmc_eventmap->end(); ++iter){
			PHHepMCGenEvent *phhepmc_event = iter->second;

			if ( phhepmc_event->get_embedding_id()!=2 ) continue;

			zvertex_in = phhepmc_event->get_collision_vertex().z();
			hz_in->Fill(zvertex_in);

			/*
			cout 
				<< "HepMC x: " << phhepmc_event->get_collision_vertex().x() << " " 
				<< "HepMC y: " << phhepmc_event->get_collision_vertex().y() << " " 
				<< "HepMC z: " << phhepmc_event->get_collision_vertex().z() << " " 
				<< endl;
				*/

		}//iter

		int n_out = 0;
		for (PHHepMCGenEventMap::ConstIter iter=_hepmc_eventmap->begin(); iter!=_hepmc_eventmap->end(); ++iter){
			PHHepMCGenEvent *phhepmc_event = iter->second;

			if ( phhepmc_event->get_embedding_id()==2 ) continue;

			float zvertex_out = phhepmc_event->get_collision_vertex().z();
			hz_out->Fill(zvertex_out);
			hz_diff->Fill(zvertex_out - zvertex_in);

			n_out++;
		}//iter

		hn_out->Fill(n_out);

	}//

	/*
	//if ( _hepmc_event ){
	if ( 0 ){

		HepMC::GenEvent *hepmc_evt = _hepmc_event->getEvent();

		for ( HepMC::GenEvent::particle_iterator p=hepmc_evt->particles_begin(); p!=hepmc_evt->particles_end(); ++p ){   

			HepMC::GenParticle *part = hepmc_evt->barcode_to_particle((*p)->barcode());

			short tmp_pid = (short)part->pdg_id();
			short tmp_status = (short)part->status();

			if ( tmp_status!=1 ){
				cout << "PID: " << tmp_pid << endl;
				continue;
			}

			//if ( abs(tmp_pid)!=211 && abs(tmp_pid)!=321 && abs(tmp_pid)!=2212 ) continue;

			gen_px[npart] = part->momentum().x();
			gen_py[npart] = part->momentum().y();
			gen_pz[npart] = part->momentum().z();

			TVector3 vec(gen_px[npart], gen_py[npart], gen_pz[npart]);
			if ( vec.Pt()<0.1 || fabs(vec.Eta())>1.0 ) continue;

			gen_x[npart] = part->production_vertex()->position().x();
			gen_y[npart] = part->production_vertex()->position().y();
			gen_z[npart] = part->production_vertex()->position().z();

			gen_pid[npart] = tmp_pid;
			gen_status[npart] = tmp_status;

			npart++;

			if ( npart>=1000 ) break;

		}//particle 
	}//hepmc_event
	*/

	/*
	bool bHF = false;

	//if ( _hepmc_event ){
	if ( 0 ){

		HepMC::GenEvent *hepmc_evt = _hepmc_event->getEvent();

		for ( HepMC::GenEvent::particle_iterator p=hepmc_evt->particles_begin(); p!=hepmc_evt->particles_end(); ++p ){   

			HepMC::GenParticle *part = hepmc_evt->barcode_to_particle((*p)->barcode());

			short tmp_pid = (short)part->pdg_id();
			short tmp_status = (short)part->status();

			if ( tmp_status!=1 ){
				if ( (abs(tmp_pid)>400 && abs(tmp_pid)<600) || (abs(tmp_pid)>4000 && abs(tmp_pid)<6000) ){
					bHF = true;
					break;
				}
			}
		}//for
	}//_hepmc_event

	if ( bHF ){
		cout << "SKIP HF event" << endl;
		return;
	}
	*/

	//cout << "NPART: " << npart << endl;

	if ( _truth_container ){
	//if ( 0 ){
	
		PHG4VtxPoint* first_point = _truth_container->GetPrimaryVtx(_truth_container->GetPrimaryVertexIndex());
		vtx_gen[0] = first_point->get_x();
		vtx_gen[1] = first_point->get_y();
		vtx_gen[2] = first_point->get_z();

		PHG4TruthInfoContainer::ConstRange range =  _truth_container->GetParticleRange();

		npart = 0;
		for (PHG4TruthInfoContainer::ConstIterator iter=range.first; iter!=range.second; ++iter)
		{

			set<int> maps_layer;
			set<int> intt_layer;
			set<int> tpc_layer;
			maps_layer.clear();
			intt_layer.clear();
			tpc_layer.clear();

			PHG4Particle *g4particle = iter->second;
			//PHG4VtxPoint *vtx = trutheval->get_vertex(g4particle);

			if ( !g4particle ) continue;
			if ( trutheval->is_primary(g4particle)==0 ) continue;
			short pid = g4particle->get_pid();
			if ( abs(pid)!=211 && abs(pid)!=321 && abs(pid)!=2212 ) continue;

			int embedding_id = trutheval->get_embed(g4particle);

			short nmaps = 0, nintt = 0, ntpc = 0;
			short nmaps_layer = 0, nintt_layer = 0, ntpc_layer = 0;
			std::set<SvtxCluster*> g4clusters = clustereval->all_clusters_from(g4particle);
			for (const SvtxCluster* g4cluster:g4clusters){
				int layer = g4cluster->get_layer();
				float clus_x = g4cluster->get_x();
				float clus_y = g4cluster->get_y();
				float clus_z = g4cluster->get_z();
				if ( layer<_n_maps_layers ){
					nmaps++;
					maps_layer.insert(layer);
				}else if ( layer<(_n_maps_layers + _n_intt_layers) ){
					/*
					cout << "Truth cluster, layer: " << layer
						<< ", r: " << sqrt(clus_x*clus_x + clus_y*clus_y) 
						<< endl;
					*/
					hclus_intt_zphi_gen[(layer-_n_maps_layers)/2]->Fill(clus_z, atan2(clus_y, clus_x));
					nintt++;
					intt_layer.insert((layer-_n_maps_layers)/2);
				}else{
					ntpc++;
					tpc_layer.insert(layer);
				}
			}

			nmaps_layer = maps_layer.size();
			nintt_layer = intt_layer.size();
			ntpc_layer = tpc_layer.size();

			TVector3 vec(g4particle->get_px(), g4particle->get_py(), g4particle->get_pz());

			if ( embedding_id==2 ){
				hpT_eta_in->Fill(vec.Pt(), vec.Eta());

				if ( fabs(vec.Eta())<1.0 ){
					hnclus_pT_tpc->Fill(vec.Pt(), ntpc);
					hnclus_pT_intt->Fill(vec.Pt(), nintt);
					hnclus_pT_mvtx->Fill(vec.Pt(), nmaps);
				}

				if ( ntpc_layer>=35 ){
					hpT_eta_in_tpc->Fill(vec.Pt(), vec.Eta());

					if ( fabs(vec.Eta())<1.0 ){
						hnclus_pT_intt_tpc35->Fill(vec.Pt(), nintt_layer);
						hnclus_pT_mvtx_tpc35->Fill(vec.Pt(), nmaps_layer);
					}

					if ( nmaps_layer>=2 ){
						hpT_eta_in_tpc_mvtx->Fill(vec.Pt(), vec.Eta());
						if ( nintt_layer>=1 )
							hpT_eta_in_tpc_mvtx_intt1->Fill(vec.Pt(), vec.Eta());
						if ( nintt_layer>=2 )
							hpT_eta_in_tpc_mvtx_intt2->Fill(vec.Pt(), vec.Eta());
						if ( nintt_layer>=3 )
							hpT_eta_in_tpc_mvtx_intt3->Fill(vec.Pt(), vec.Eta());
						if ( nintt_layer>=4 )
							hpT_eta_in_tpc_mvtx_intt4->Fill(vec.Pt(), vec.Eta());
					}
				}
			}else{
				hpT_eta_out->Fill(vec.Pt(), vec.Eta());

				if ( ntpc_layer>=35 ){
					hpT_eta_out_tpc->Fill(vec.Pt(), vec.Eta());
					if ( nmaps_layer>=2 ){
						hpT_eta_out_tpc_mvtx->Fill(vec.Pt(), vec.Eta());
						if ( nintt_layer>=1 )
							hpT_eta_out_tpc_mvtx_intt1->Fill(vec.Pt(), vec.Eta());
						if ( nintt_layer>=2 )
							hpT_eta_out_tpc_mvtx_intt2->Fill(vec.Pt(), vec.Eta());
						if ( nintt_layer>=3 )
							hpT_eta_out_tpc_mvtx_intt3->Fill(vec.Pt(), vec.Eta());
						if ( nintt_layer>=4 )
							hpT_eta_out_tpc_mvtx_intt4->Fill(vec.Pt(), vec.Eta());
					}
				}

				continue;
			}

			//if ( trutheval->get_embed(g4particle)<=1 ) continue;

			gen_pid[npart]  = pid;
			gen_px[npart] = g4particle->get_px();
			gen_py[npart] = g4particle->get_py();
			gen_pz[npart] = g4particle->get_pz();

			PHG4VtxPoint* vtx = trutheval->get_vertex(g4particle);

			if ( vtx ){
				gen_vx[npart] = vtx->get_x();
				gen_vy[npart] = vtx->get_y();
				gen_vz[npart] = vtx->get_z();
			}else{
				gen_vx[npart] = -999;
				gen_vy[npart] = -999;
				gen_vz[npart] = -999;
			}

			gen_embed[npart] = trutheval->get_embed(g4particle);
			gen_primary[npart] = trutheval->is_primary(g4particle);


			gen_ngmaps[npart] = nmaps;
			gen_ngintt[npart] = nintt;
			gen_ngtpc[npart] = ntpc;

			gen_ngmaps_layer[npart] = nmaps_layer;
			gen_ngintt_layer[npart] = nintt_layer;
			gen_ngtpc_layer[npart] = ntpc_layer;

			npart++;
		}
	}//truth_container

	if ( _vtxmap ){
	//if ( 0 ){

		for (SvtxVertexMap::ConstIter iter=_vtxmap->begin(); iter!=_vtxmap->end(); ++iter){
			SvtxVertex *vtx = iter->second;

			vtx_id[nvertex] = vtx->get_id();

			vtx_reco[nvertex][0] = vtx->get_x();
			vtx_reco[nvertex][1] = vtx->get_y();
			vtx_reco[nvertex][2] = vtx->get_z();

			vtx_reco_err[nvertex][0] = vtx->get_error(0,0);
			vtx_reco_err[nvertex][1] = vtx->get_error(1,1);
			vtx_reco_err[nvertex][2] = vtx->get_error(2,2);

			/*
			PHG4VtxPoint* point = vertexeval->max_truth_point_by_ntracks(vtx);

			if ( point ){
				vtx_gen[nvertex][0] = point->get_x();
				vtx_gen[nvertex][1] = point->get_y();
				vtx_gen[nvertex][2] = point->get_z();

				cout 
					<< "Eval x: " << vtx_gen[0] << " " 
					<< "Eval y: " << vtx_gen[1] << " " 
					<< "Eval z: " << vtx_gen[2] << " " 
					<< endl;
			}
			*/

			//cout << "id: " << vtx->get_id() << ", dx: " << vtx_reco[nvertex][0]-vtx_gen[nvertex][0] << ", dy: " << vtx_reco[nvertex][1]-vtx_gen[nvertex][1] << ", dz: " << vtx_reco[nvertex][2]-vtx_gen[nvertex][2] << endl;
			nvertex++;
		}
	}//vtxmap

			/*

	//if ( _clustermap ){
	if ( 0 ){
		//bool bFIRST = true;
		for (SvtxClusterMap::ConstIter iter=_clustermap->begin(); iter!=_clustermap->end(); ++iter){
			SvtxCluster *svtx_cluster = iter->second;

			cluster_layer[ncluster] = svtx_cluster->get_layer();

			//if ( cluster_layer[ncluster]>=3 ) continue;

			for (int ii=0; ii<3; ii++){
				cluster_err[ncluster][ii] = svtx_cluster->get_error(ii,ii);
			}

			PHG4Particle* g4particle = clustereval->max_truth_particle_by_energy(svtx_cluster);

			cluster_gprimary[ncluster] = trutheval->is_primary(g4particle);
			cluster_gflavor[ncluster] = g4particle->get_pid();

			PHG4Hit *g4hit = clustereval->max_truth_hit_by_energy(svtx_cluster);
			if ( g4hit ){
				PHG4Particle *g4particle = trutheval->get_particle(g4hit);

				if ( g4particle ){
					cluster_gembed[ncluster] = trutheval->get_embed(g4particle);
				}
			}

			cluster_gpx[ncluster] = g4particle->get_px();
			cluster_gpy[ncluster] = g4particle->get_py();

			PHG4VtxPoint* vtx = trutheval->get_vertex(g4particle);
			if (vtx){
				cluster_gvx[ncluster] = vtx->get_x();
				cluster_gvy[ncluster] = vtx->get_y();
			} 

			if ( !_hepmc_event && cluster_gprimary[ncluster]==1 && bFIRST ){
				gen_px[npart] = g4particle->get_px();
				gen_py[npart] = g4particle->get_py();
				gen_pz[npart] = g4particle->get_pz();
				gen_pid[npart] = g4particle->get_pid();
				gen_status[npart] = 1;
				npart++;
				bFIRST = false;
			}

			ncluster++;
			if ( ncluster>=5000 ) break;
		}
	}
			*/

	//cout << "CHECK1" << endl;
	if ( _trkmap ){
	//if ( 0 ){
		for (SvtxTrackMap::ConstIter iter=_trkmap->begin(); iter!=_trkmap->end(); ++iter){
			SvtxTrack *svtx_trk = iter->second;

			PHG4Particle* g4particle = trackeval->max_truth_particle_by_nclusters(svtx_trk);

			if ( g4particle ){
				track_gprimary[ntrack] = trutheval->is_primary(g4particle); 
				track_gembed[ntrack] = trutheval->get_embed(g4particle); 
				track_gflavor[ntrack] = g4particle->get_pid();
				track_parent_id[ntrack] = g4particle->get_parent_id();

				if ( track_parent_id[ntrack]!=0 ){
					PHG4Particle* g4particle_p = _truth_container->GetParticle(track_parent_id[ntrack]);
					if ( g4particle_p ){
						track_parent_id[ntrack] = g4particle_p->get_pid();
					}else{
						track_parent_id[ntrack] = -999;
					}
				}

				track_matched_cluster_n[ntrack] = trackeval->get_nclusters_contribution(svtx_trk, g4particle);

				PHG4VtxPoint* vtx = trutheval->get_vertex(g4particle);
				if ( vtx ){
					track_gvx[ntrack] = vtx->get_x();
					track_gvy[ntrack] = vtx->get_y();
					track_gvz[ntrack] = vtx->get_z();
				}else{
					track_gvx[ntrack] = -999;
					track_gvy[ntrack] = -999;
					track_gvz[ntrack] = -999;
				}

				track_gpx[ntrack] = g4particle->get_px();
				track_gpy[ntrack] = g4particle->get_py();
				track_gpz[ntrack] = g4particle->get_pz();
			}//g4particle

			/*
			int trk_cluster_n = 0;
			for (SvtxTrack::ConstClusterIter iter2 = svtx_trk->begin_clusters(); iter2!=svtx_trk->end_clusters(); ++iter2) {
				unsigned int cluster_id = *iter2;
				SvtxCluster* cluster = _clustermap->get(cluster_id);

				track_cluster_layer[ntrack][trk_cluster_n] = cluster->get_layer();

				PHG4Hit *g4hit = clustereval->max_truth_hit_by_energy(cluster);
				if ( g4hit ){
					PHG4Particle *g4particle = trutheval->get_particle(g4hit);

					if ( g4particle ){
						track_cluster_gembed[ntrack][trk_cluster_n] = trutheval->get_embed(g4particle);
					}
				}

				trk_cluster_n++;

				if ( trk_cluster_n>=100 ) break;
			}

			if ( trk_cluster_n!=int(svtx_trk->size_clusters()) ){
				cout << "Cluster sizes are not matched!!" << endl;
			}
			*/

			short nmaps = 0, nintt = 0, ntpc = 0;
			short nmaps_truth = 0, nintt_truth = 0, ntpc_truth = 0;
			short hit_pattern_maps = 0, hit_pattern_intt = 0;

			set<int> maps_layer;
			set<int> intt_layer;
			set<int> tpc_layer;
			maps_layer.clear();
			intt_layer.clear();
			tpc_layer.clear();

			for (SvtxTrack::ConstClusterIter iter2 = svtx_trk->begin_clusters(); iter2!=svtx_trk->end_clusters(); ++iter2) {
				unsigned int cluster_id = *iter2;
				SvtxCluster* cluster = _clustermap->get(cluster_id);
				short layer = cluster->get_layer();
				float clus_x = cluster->get_x();
				float clus_y = cluster->get_y();
				float clus_z = cluster->get_z();

				PHG4Hit *g4hit = clustereval->max_truth_hit_by_energy(cluster);
				PHG4Particle *g4particle = trutheval->get_particle(g4hit);

				short gembed_id = -999;
				if ( g4particle ){
					gembed_id = trutheval->get_embed(g4particle);
				}

				if ( layer<_n_maps_layers ){
					nmaps++;
					if ( gembed_id==2 ) nmaps_truth++;
					hit_pattern_maps += pow(2,layer);

					maps_layer.insert(layer);
				}else if ( layer<(_n_maps_layers + _n_intt_layers) ){
					nintt++;
					if ( gembed_id==2 ) nintt_truth++;
					hit_pattern_intt += pow(2,layer-_n_maps_layers);

					hclus_intt_zphi_reco[(layer-_n_maps_layers)/2]->Fill(clus_z, atan2(clus_y, clus_x));

					intt_layer.insert((layer-_n_maps_layers)/2);
				}else{
					ntpc++;
					if ( gembed_id==2 ) ntpc_truth++;

					tpc_layer.insert(layer);
				}
			}

			track_cluster_n[ntrack] = svtx_trk->size_clusters();
			track_cluster_n_maps[ntrack] = nmaps;
			track_cluster_n_intt[ntrack] = nintt;
			track_cluster_n_tpc[ntrack] = ntpc;
			track_cluster_n_layer_maps[ntrack] = maps_layer.size();
			track_cluster_n_layer_intt[ntrack] = intt_layer.size();
			track_cluster_n_layer_tpc[ntrack] = tpc_layer.size();
			track_cluster_n_maps_truth[ntrack] = nmaps_truth;
			track_cluster_n_intt_truth[ntrack] = nintt_truth;
			track_cluster_n_tpc_truth[ntrack] = ntpc_truth;
			track_cluster_pattern_maps[ntrack] = hit_pattern_maps;
			track_cluster_pattern_intt[ntrack] = hit_pattern_intt;
			//track_cluster_n[ntrack] = trk_cluster_n;
			track_pt[ntrack] = svtx_trk->get_pt();
			track_eta[ntrack] = svtx_trk->get_eta();
			track_chiq[ntrack] = svtx_trk->get_chisq();
			track_ndf[ntrack] = short(svtx_trk->get_ndf());

			track_dca2d[ntrack] = svtx_trk->get_dca2d();
			track_dca3d_xy[ntrack] = svtx_trk->get_dca3d_xy();
			track_dca3d_z[ntrack] = svtx_trk->get_dca3d_z();

			ntrack++;

			if ( ntrack>=2000 ) break;
		}
	}

	//cout << "CHECK2" << endl;

	/*
	if ( _trkmap_refit ){
	//if ( 0 ){
		for (SvtxTrackMap::ConstIter iter=_trkmap_refit->begin(); iter!=_trkmap_refit->end(); ++iter){

			SvtxTrack *trk_refit = iter->second;

			trk_id[ntrack] = trk_refit->get_id();
			trk_chi2pdf[ntrack] = trk_refit->get_chisq()/trk_refit->get_ndf();
			trk_quality[ntrack] = trk_refit->get_quality();

			trk_x[ntrack] = trk_refit->get_x();
			trk_y[ntrack] = trk_refit->get_y();
			trk_z[ntrack] = trk_refit->get_z();
			trk_px[ntrack] = trk_refit->get_px();
			trk_py[ntrack] = trk_refit->get_py();
			trk_pz[ntrack] = trk_refit->get_pz();

			PHG4Particle *best_mc = trackeval->max_truth_particle_by_nclusters(trk_refit);

			if ( best_mc ){
				trk_mc_p_pid[ntrack] = best_mc->get_parent_id();
				trk_mc_pid[ntrack] = best_mc->get_pid();
				trk_mc_px[ntrack] = best_mc->get_px();
				trk_mc_py[ntrack] = best_mc->get_py();
				trk_mc_pz[ntrack] = best_mc->get_pz();
				trk_mc_e[ntrack] = best_mc->get_e();

				//cout << "PID: " << trk_mc_pid[ntrack] << ", P ID: " << trk_mc_p_pid[ntrack] << ", G ID: " << best_mc->get_primary_id() << endl;

				if ( _truth_container ){
					int vtx_id = best_mc->get_vtx_id();
					PHG4VtxPoint *mc_vtx = _truth_container->GetVtx(vtx_id);
					trk_mc_x[ntrack] = mc_vtx->get_x(); 
					trk_mc_y[ntrack] = mc_vtx->get_y(); 
					trk_mc_z[ntrack] = mc_vtx->get_z(); 
				}
			}//best_mc

			ntrack++;

			if ( ntrack>=2000 ) break;
		}
		//cout << "-------------------------------------------------------------" << endl;
	}//trkmap_refit
	*/

	int count = 0;

	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(&_event);

	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(vtx_gen);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(&nvertex);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(vtx_reco);
	//((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(vtx_reco_err);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(vtx_id);

	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(&npart);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(gen_pid);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(gen_embed);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(gen_primary);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(gen_px);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(gen_py);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(gen_pz);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(gen_vx);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(gen_vy);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(gen_vz);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(gen_ngmaps);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(gen_ngintt);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(gen_ngtpc);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(gen_ngmaps_layer);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(gen_ngintt_layer);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(gen_ngtpc_layer);

	//((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(&ncluster);
	//((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(cluster_layer);
	//((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(cluster_gprimary);
	//((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(cluster_gembed);
	//((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(cluster_gflavor);
	//((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(cluster_gpx);
	//((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(cluster_gpy);
	//((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(cluster_gvx);
	//((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(cluster_gvy);
	//((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(cluster_err);

	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(&ntrack);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_pt);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_eta);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_chiq);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_ndf);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_cluster_n_maps);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_cluster_n_intt);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_cluster_n_tpc);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_cluster_n_layer_maps);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_cluster_n_layer_intt);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_cluster_n_layer_tpc);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_cluster_n_maps_truth);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_cluster_n_intt_truth);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_cluster_n_tpc_truth);
	//((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_cluster_pattern_maps);
	//((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_cluster_pattern_intt);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_gprimary);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_gembed);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_gflavor);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_parent_id);
	//((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_matched_cluster_n);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_gvx);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_gvy);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_gvz);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_gpx);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_gpy);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_gpz);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_dca2d);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_dca3d_xy);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_dca3d_z);

	_events->Fill();

  return;

}

//----------------------------------------------------------------------------//
//-- reset_variables():
//--   Reset all the tree variables to their default values.
//--   Needs to be called at the start of every event
//----------------------------------------------------------------------------//
void AnaTrackingEff::reset_variables()
{

	nvertex = 0;
	vtx_gen[0] = vtx_gen[1] = vtx_gen[2] = -9999;
	for (int ii=0; ii<10; ii++){
		vtx_id[ii] = -9999;
		for (int iv=0; iv<3; iv++){
			vtx_reco[ii][iv] = -9999;
			vtx_reco_err[ii][iv] = -9999;
		}
	}

	/*
	for (int ic=0; ic<5000; ic++){
		cluster_layer[ic] = cluster_gembed[ic] = -1;
		cluster_gprimary[ic] = cluster_gflavor[ic] = 0;
		cluster_gpx[ic] = cluster_gpy[ic] = -9999;
		cluster_gvx[ic] = cluster_gvy[ic] = -9999;
		for (int ii=0; ii<3; ii++){
			cluster_err[ic][ii] = 0;
		}
	}
	*/

	for (int it=0; it<2000; it++){

		gen_status[it] = gen_pid[it] = 0;
		gen_primary[it] = gen_embed[it] = 0; 
		gen_px[it] = gen_py[it] = gen_pz[it] = -9999;
		gen_vx[it] = gen_vy[it] = gen_vz[it] = -9999;
		gen_ngmaps[it] = gen_ngintt[it] = gen_ngtpc[it] = 0;
		gen_ngmaps_layer[it] = gen_ngintt_layer[it] = gen_ngtpc_layer[it] = 0;

		track_eta[it] = track_pt[it] = track_chiq[it] = -9999;
		track_ndf[it] = 0;
		track_cluster_n[it] = track_matched_cluster_n[it] = track_cluster_n_layer[it] = 0; 
		track_cluster_n_maps[it] = track_cluster_n_intt[it] = track_cluster_n_tpc[it] = 0;
		track_cluster_n_layer_maps[it] = track_cluster_n_layer_intt[it] = track_cluster_n_layer_tpc[it] = 0;
		track_cluster_n_maps_truth[it] = track_cluster_n_intt_truth[it] = track_cluster_n_tpc_truth[it] = 0;
		track_cluster_pattern_maps[it] = 0;
		track_cluster_pattern_intt[it] = 0;
		track_gprimary[it] = track_gembed[it] = 0;
		track_gflavor[it] = track_parent_id[it] = -999;
		track_gvx[it] = track_gvy[it] = track_gvz[it] = -999;
		track_gpx[it] = track_gpy[it] = track_gpz[it] = -999;
		track_dca2d[it] = track_dca3d_xy[it] = track_dca3d_z[it] = -999;
		/*
		for (int il=0; il<100; il++){
			track_cluster_layer[it][il] = -1;
			track_cluster_gembed[it][il] = -1;
		}
		*/
	}

}

//----------------------------------------------------------------------------//
//-- GetNodes():
//--   Get all the all the required nodes off the node tree
//----------------------------------------------------------------------------//
void AnaTrackingEff::GetNodes(PHCompositeNode * topNode)
{
  //DST objects
	//
  //PHG4TruthInfoContainer
  _truth_container = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
  if (!_truth_container && _event<2)
  {
    cout << PHWHERE << " PHG4TruthInfoContainer node not found on node tree" << endl;
  }

  //HepMCGenEventMap
  _hepmc_eventmap = findNode::getClass<PHHepMCGenEventMap>(topNode,"PHHepMCGenEventMap");
  if (!_hepmc_eventmap && _event<2)
  {
    cout << PHWHERE << " PHHepMCGenEventMap node not found on node tree" << endl;
  }

  //SvtxVertexMap
  _vtxmap = findNode::getClass<SvtxVertexMap>(topNode,"SvtxVertexMap");
  if (!_vtxmap && _event<2)
  {
    cout << PHWHERE << " SvtxVertexMap node not found on node tree" << endl;
  }

  //SvtxTrackMap
  _trkmap = findNode::getClass<SvtxTrackMap>(topNode,"SvtxTrackMap");
  if (!_trkmap && _event<2)
  {
    cout << PHWHERE << " SvtxTrackMap node not found on node tree" << endl;
  }

  //SvtxClusterMap
  _clustermap = findNode::getClass<SvtxClusterMap>(topNode,"SvtxClusterMap");
  if (!_clustermap && _event<2)
  {
    cout << PHWHERE << " SvtxClusterMap node not found on node tree" << endl;
  }

  //SvtxTrackMap
	/*
  _trkmap_refit = findNode::getClass<SvtxTrackMap>(topNode,"SvtxTrackMapRefit");
  if (!_trkmap_refit && _event<2)
  {
    cout << PHWHERE << " SvtxTrackMapRefit node not found on node tree" << endl;
  }
	*/

}



