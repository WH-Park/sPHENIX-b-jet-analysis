////////////////////////////////////////////////////////////////////////////////
//
// This module is desgined to grab svtx tracks and put truth and cluster
// information into a TTree for GenFit testing
//
////////////////////////////////////////////////////////////////////////////////
//
// Darren McGlinchey
// 1 Apr 2016
//
////////////////////////////////////////////////////////////////////////////////


#include "AnaSvtxClus.h"

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

#include <g4detectors/PHG4CellContainer.h>
#include <g4detectors/PHG4Cell.h>

#include <TTree.h>
#include <TVector3.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#include <HepMC/GenParticle.h>

#include <HFJetTruthGeneration/HFJetDefs.h>

#include <iostream>

using namespace std;

//----------------------------------------------------------------------------//
//-- Constructor:
//--  simple initialization
//----------------------------------------------------------------------------//
AnaSvtxClus::AnaSvtxClus(const string &name):
  SubsysReco( name ),
  _events( NULL )
{
  //initialize
  _event = 0;
  _outfile = "AnaSvtxClus.root";

	_truth_container = NULL;
	_g4cells_maps = NULL;
	_hepmc_event = NULL;

	_trkmap = NULL;
	_clustermap = NULL;
	_svtxevalstack = NULL;

}

//----------------------------------------------------------------------------//
//-- Init():
//--   Intialize all histograms, trees, and ntuples
//----------------------------------------------------------------------------//
int AnaSvtxClus::Init(PHCompositeNode *topNode)
{
  cout << PHWHERE << " Openning file " << _outfile << endl;
  PHTFileServer::get().open( _outfile, "RECREATE");

  // create TTree
  //_events = new TTree("events", "Svtx Event");

	//_events->Branch("vtx_reco",0,"vtx_reco[3]/F");
	//_events->Branch("vtx_reco_err",0,"vtx_reco_err[3]/F");

	/*
	_events->Branch("npart",0,"npart/I");
	_events->Branch("gen_pid",0,"gen_pid[npart]/S");
	_events->Branch("gen_status",0,"gen_status[npart]/S");
	_events->Branch("gen_px",0,"gen_px[npart]/F");
	_events->Branch("gen_py",0,"gen_py[npart]/F");
	_events->Branch("gen_pz",0,"gen_pz[npart]/F");
	_events->Branch("gen_x",0,"gen_x[npart]/F");
	_events->Branch("gen_y",0,"gen_y[npart]/F");
	_events->Branch("gen_z",0,"gen_z[npart]/F");
	*/

	/*
	_events->Branch("ng4hit",0,"ng4hit/I");
	_events->Branch("g4hit_layer",0,"g4hit_layer[ng4hit]/S");
	_events->Branch("g4hit_stave",0,"g4hit_stave[ng4hit]/S");
	_events->Branch("g4hit_half_stave",0,"g4hit_half_stave[ng4hit]/S");
	_events->Branch("g4hit_module",0,"g4hit_module[ng4hit]/S");
	_events->Branch("g4hit_chip",0,"g4hit_chip[ng4hit]/S");
	*/


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

	/*
	_events->Branch("ntrack",0,"ntrack/I");
	_events->Branch("track_pt",0,"track_pt[ntrack]/F");
	_events->Branch("track_eta",0,"track_eta[ntrack]/F");
	_events->Branch("track_chiq",0,"track_chiq[ntrack]/F");
	_events->Branch("track_ndf",0,"track_ndf[ntrack]/S");
	_events->Branch("track_cluster_n",0,"track_cluster_n[ntrack]/S");
	_events->Branch("track_cluster_layer",0,"track_cluster_layer[ntrack][100]/S");
	_events->Branch("track_cluster_gembed",0,"track_cluster_gembed[ntrack][100]/S");
	_events->Branch("track_gprimary",0,"track_gprimary[ntrack]/S");
	_events->Branch("track_gembed",0,"track_gembed[ntrack]/S");
	_events->Branch("track_gflavor",0,"track_gflavor[ntrack]/S");
	_events->Branch("track_parent_id",0,"track_parent_id[ntrack]/S");
	_events->Branch("track_matched_cluster_n",0,"track_matched_cluster_n[ntrack]/S");
	_events->Branch("track_gvx",0,"track_gvx[ntrack]/F");
	_events->Branch("track_gvy",0,"track_gvy[ntrack]/F");
	_events->Branch("track_gvz",0,"track_gvz[ntrack]/F");
	_events->Branch("track_gpx",0,"track_gpx[ntrack]/F");
	_events->Branch("track_gpy",0,"track_gpy[ntrack]/F");
	_events->Branch("track_gpz",0,"track_gpz[ntrack]/F");
	_events->Branch("track_dca2d",0,"track_dca2d[ntrack]/F");
	_events->Branch("track_dca3d_xy",0,"track_dca3d_xy[ntrack]/F");
	_events->Branch("track_dca3d_z",0,"track_dca3d_z[ntrack]/F");
	*/

	_t_cell = new TTree("_t_cell","hit info");

  _t_cell->Branch("event",&_event,"event/I");
	_t_cell->Branch("g4cell_layer",&g4cell_layer,"g4cell_layer/S");
	_t_cell->Branch("g4cell_stave",&g4cell_stave,"g4cell_stave/S");
	_t_cell->Branch("g4cell_chip",&g4cell_chip,"g4cell_chip/S");
	_t_cell->Branch("g4cell_pixel",&g4cell_pixel,"g4cell_pixel/I");
	_t_cell->Branch("g4cell_ind_x",&g4cell_ind_x,"g4cell_ind_x/S");
	_t_cell->Branch("g4cell_ind_z",&g4cell_ind_z,"g4cell_ind_z/S");
	_t_cell->Branch("g4cell_primary",&g4cell_primary,"g4cell_primary/S");

  return 0;
}

//----------------------------------------------------------------------------//
//-- process_event():
//--   Call user instructions for every event.
//--   This function contains the analysis structure.
//----------------------------------------------------------------------------//
int AnaSvtxClus::process_event(PHCompositeNode *topNode)
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
int AnaSvtxClus::End(PHCompositeNode *topNode)
{

	cout << "-----AnaSvtxClus::End------" << endl;

  PHTFileServer::get().cd( _outfile );
	if ( _events )
		_events->Write();
	if ( _t_cell )
		_t_cell->Write();
  //PHTFileServer::get().close();

	//if ( trutheval ) trutheval->Delete();
	//if ( trackeval ) trackeval->Delete();

  return 0;
}


//----------------------------------------------------------------------------//
//-- fill_tree():
//--   Fill the trees with truth, track fit, and cluster information
//----------------------------------------------------------------------------//
void AnaSvtxClus::fill_tree(PHCompositeNode *topNode)
{

  // Make sure to reset all the TTree variables before trying to set them.
  reset_variables();

	//cout << "reset_variables" << endl;

	SvtxClusterEval *clustereval = _svtxevalstack->get_cluster_eval();
	SvtxTrackEval *trackeval = _svtxevalstack->get_track_eval();
	SvtxTruthEval *trutheval = _svtxevalstack->get_truth_eval();
	SvtxHitEval *hiteval = _svtxevalstack->get_hit_eval();

	int ng4hit = 0;
	int ng4cell = 0;
	int ncluster = 0;
	int ntrack = 0;
	int npart = 0;

	if ( 0 ){
		std::set<PHG4Hit*> g4hits = trutheval->all_truth_hits();
		for (std::set<PHG4Hit*>::iterator iter = g4hits.begin(); iter != g4hits.end(); ++iter){
			PHG4Hit *g4hit = *iter;

			g4hit_layer[ng4hit] = g4hit->get_layer();

			if ( g4hit_layer[ng4hit]>2 ) continue;

			g4hit_stave[ng4hit] = g4hit->get_property_int(PHG4Hit::prop_stave_index);
			g4hit_half_stave[ng4hit] = g4hit->get_property_int(PHG4Hit::prop_half_stave_index);
			g4hit_module[ng4hit] = g4hit->get_property_int(PHG4Hit::prop_module_index);
			g4hit_chip[ng4hit] = g4hit->get_property_int(PHG4Hit::prop_chip_index);


			ng4hit++;

			if ( ng4hit>=5000 ) break;
		}
	}

	if ( _hitmap ){

		for (SvtxHitMap::ConstIter iter=_hitmap->begin(); iter!=_hitmap->end(); ++iter){
			SvtxHit *svtx_hit = iter->second;

			PHG4Cell* g4cell = nullptr;
			if ( _g4cells_maps ) g4cell = _g4cells_maps->findCell(svtx_hit->get_cellid());

			if ( g4cell ){

				g4cell_layer = g4cell->get_layer();
				g4cell_stave = g4cell->get_stave_index();
				g4cell_chip = g4cell->get_chip_index();
				g4cell_pixel = g4cell->get_pixel_index();
				g4cell_ind_x = g4cell_pixel/1024;
				g4cell_ind_z = g4cell_pixel%1024;
				//g4cell_ind_x = g4cell->get_ladder_phi_index();
				//g4cell_ind_z = g4cell->get_ladder_z_index();

				g4cell_primary = 0;

				PHG4Hit* g4hit = hiteval->max_truth_hit_by_energy(svtx_hit);
				if ( g4hit ){
					PHG4Particle* g4particle = trutheval->get_particle(g4hit);

					if ( g4particle ){
						g4cell_primary = trutheval->is_primary(g4particle);
					}
				}

				_t_cell->Fill();

			}

		}
	}//_hitmap

	//if ( _g4cells_maps ){
	if ( 0 ){
		PHG4CellContainer::ConstRange cell_begin_end = _g4cells_maps->getCells();
		for (PHG4CellContainer::ConstIterator citer = cell_begin_end.first; citer != cell_begin_end.second; ++citer)
		{   

			PHG4Cell *g4cell = citer->second;

			g4cell_layer = g4cell->get_layer();
			g4cell_stave = g4cell->get_stave_index();
			g4cell_chip = g4cell->get_chip_index();
			g4cell_pixel = g4cell->get_pixel_index();
			g4cell_ind_x = g4cell_pixel%512;
			g4cell_ind_z = g4cell_pixel/512;

			_t_cell->Fill();

			ng4cell++;

		}   
	}

	//cout << "ready" << endl;


	//cout << "-EVENT: " << _event << "---------------------------------------------------------------------------------------------------" << endl;

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

	//cout << "NPART: " << npart << endl;

	/*
	//if ( !_hepmc_event && _truth_container ){
	if ( 0 ){

		PHG4TruthInfoContainer::ConstRange range =  _truth_container->GetParticleRange();

		npart = 0;
		for (PHG4TruthInfoContainer::ConstIterator iter=range.first; iter!=range.second; ++iter)
		{

			PHG4Particle *g4particle = iter->second;
			PHG4VtxPoint *vtx = trutheval->get_vertex(g4particle);

			gen_pid[npart]  = g4particle->get_pid();
			gen_prod_x[npart] = vtx->get_x();
			gen_prod_y[npart] = vtx->get_y();
			gen_prod_z[npart] = vtx->get_z();

			npart++;
		}
	}//truth_container
	*/

	if ( 0 ){

		for (SvtxVertexMap::ConstIter iter=_vtxmap->begin(); iter!=_vtxmap->end(); ++iter){
			SvtxVertex *vtx = iter->second;

			vtx_reco[0] = vtx->get_x();
			vtx_reco[1] = vtx->get_y();
			vtx_reco[2] = vtx->get_z();

			vtx_reco_err[0] = vtx->get_error(0,0);
			vtx_reco_err[1] = vtx->get_error(1,1);
			vtx_reco_err[2] = vtx->get_error(2,2);
		}
	}//vtxmap


	//if ( _clustermap ){
	if ( 0 ){
		bool bFIRST = true;
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

	//if ( _trkmap ){
	if ( 0 ){
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
					track_parent_id[ntrack] = g4particle_p->get_pid();
				}

				track_matched_cluster_n[ntrack] = trackeval->get_nclusters_contribution(svtx_trk, g4particle);

				PHG4VtxPoint* vtx = trutheval->get_vertex(g4particle);
				track_gvx[ntrack] = vtx->get_x();
				track_gvy[ntrack] = vtx->get_y();
				track_gvz[ntrack] = vtx->get_z();

				track_gpx[ntrack] = g4particle->get_px();
				track_gpy[ntrack] = g4particle->get_py();
				track_gpz[ntrack] = g4particle->get_pz();
			}

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

			track_cluster_n[ntrack] = trk_cluster_n;
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

	if ( _events ){
		int count = 0;
		((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(&_event);

		//((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(vtx_reco);
		//((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(vtx_reco_err);
		//
		/*
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(&ng4hit);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(g4hit_layer);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(g4hit_stave);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(g4hit_half_stave);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(g4hit_module);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(g4hit_chip);
			 */

		/*

			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(&ng4cell);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(g4cell_layer);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(g4cell_stave);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(g4cell_chip);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(g4cell_pixel);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(g4cell_ind_x);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(g4cell_ind_z);

			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(&npart);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(gen_pid);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(gen_status);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(gen_px);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(gen_py);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(gen_pz);

			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(gen_x);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(gen_y);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(gen_z);
			 */

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

		/*
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(&ntrack);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_pt);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_eta);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_chiq);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_ndf);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_cluster_n);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_cluster_layer);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_cluster_gembed);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_gprimary);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_gembed);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_gflavor);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_parent_id);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_matched_cluster_n);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_gvx);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_gvy);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_gvz);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_gpx);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_gpy);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_gpz);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_dca2d);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_dca3d_xy);
			 ((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(track_dca3d_z);
			 */

		_events->Fill();
	}

  return;

}

//----------------------------------------------------------------------------//
//-- reset_variables():
//--   Reset all the tree variables to their default values.
//--   Needs to be called at the start of every event
//----------------------------------------------------------------------------//
void AnaSvtxClus::reset_variables()
{

	for (int iv=0; iv<3; iv++){
		vtx_reco[iv] = -9999;
		vtx_reco_err[iv] = -9999;
	}

	for (int ic=0; ic<5000; ic++){
		cluster_layer[ic] = cluster_gembed[ic] = -1;
		cluster_gprimary[ic] = cluster_gflavor[ic] = 0;
		cluster_gpx[ic] = cluster_gpy[ic] = -9999;
		cluster_gvx[ic] = cluster_gvy[ic] = -9999;
		for (int ii=0; ii<3; ii++){
			cluster_err[ic][ii] = 0;
		}
	}

	g4cell_layer = g4cell_stave = g4cell_chip = -1;
	g4cell_ind_x = g4cell_ind_z = -1;
	g4cell_pixel = -1;
	g4cell_primary = -1;

	for (int it=0; it<2000; it++){

		gen_status[it] = gen_pid[it] = 0;
		gen_px[it] = gen_py[it] = gen_pz[it] = -9999;
		gen_x[it] = gen_y[it] = gen_z[it] = -9999;

		track_eta[it] = track_pt[it] = track_chiq[it] = -9999;
		track_ndf[it] = 0;
		track_cluster_n[it] = track_matched_cluster_n[it] = 0; 
		track_gprimary[it] = track_gembed[it] = 0;
		track_gflavor[it] = track_parent_id[it] = -999;
		track_gvx[it] = track_gvy[it] = track_gvz[it] = -999;
		track_gpx[it] = track_gpy[it] = track_gpz[it] = -999;
		track_dca2d[it] = track_dca3d_xy[it] = track_dca3d_z[it] = -999;
		for (int il=0; il<100; il++){
			track_cluster_layer[it][il] = -1;
			track_cluster_gembed[it][il] = -1;
		}
	}

}

//----------------------------------------------------------------------------//
//-- GetNodes():
//--   Get all the all the required nodes off the node tree
//----------------------------------------------------------------------------//
void AnaSvtxClus::GetNodes(PHCompositeNode * topNode)
{
  //DST objects
	//
  //PHG4TruthInfoContainer
  _truth_container = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
  if (!_truth_container && _event<2)
  {
    cout << PHWHERE << " PHG4TruthInfoContainer node not found on node tree" << endl;
  }

	//PHG4CELL
	_g4cells_maps = findNode::getClass<PHG4CellContainer>(topNode,"G4CELL_MAPS");
	if (!_g4cells_maps && _event<2)
	{
		cout << PHWHERE << " G4CELL_MAPS node not found on node tree" << endl;
	}

	_hitmap = findNode::getClass<SvtxHitMap>(topNode,"SvtxHitMap");
	if (!_hitmap && _event<2)
	{
		cout << PHWHERE << " SvtxHitMap node not found on node tree" << endl;
	}

  //HepMCGenEvent
  _hepmc_event = findNode::getClass<PHHepMCGenEvent>(topNode,"PHHepMCGenEvent");
  if (!_hepmc_event && _event<2)
  {
    cout << PHWHERE << " PHHepMCGenEvent node not found on node tree" << endl;
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



