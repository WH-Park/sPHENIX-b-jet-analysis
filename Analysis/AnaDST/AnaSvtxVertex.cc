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


#include "AnaSvtxVertex.h"

#include <phool/phool.h>
#include <phool/getClass.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4VtxPoint.h>
#include <fun4all/PHTFileServer.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack.h>
//#include <trackbase_historic/SvtxClusterMap.h>
//#include <trackbase_historic/SvtxCluster.h>
//#include <trackbase_historic/SvtxHitMap.h>
//#include <trackbase_historic/SvtxHit.h>

#include <g4eval/SvtxEvalStack.h>
#include <g4eval/SvtxTrackEval.h>
#include <g4eval/SvtxClusterEval.h>
#include <g4eval/SvtxTruthEval.h>
#include <g4eval/SvtxVertexEval.h>
#include <g4eval/SvtxHitEval.h>

#include <g4jets/JetMap.h>
#include <g4jets/Jet.h>

#include <TTree.h>
#include <TH2D.h>
#include <TVector3.h>
#include <phhepmc/PHHepMCGenEventMap.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#include <HepMC/GenParticle.h>

#include <phgenfit/Fitter.h>
#include <phfield/PHFieldUtility.h>
#include <phgeom/PHGeomUtility.h>

#include <GenFit/FitStatus.h>											// for FitStatus
#include <GenFit/GFRaveTrackParameters.h>					// for GFRaveTrackParameters
#include <GenFit/GFRaveVertex.h>
#include <GenFit/GFRaveVertexFactory.h>
#include <GenFit/KalmanFittedStateOnPlane.h>			// for KalmanFittedSTateOn..
#include <GenFit/KalmanFitterInfo.h>
#include <GenFit/MeasuredStateOnPlane.h>
#include <GenFit/RKTrackRep.h>
#include <GenFit/Track.h>
#include <GenFit/TrackPoint.h>										// for TrackPoint

#include <TMatrixDSymfwd.h>												// for TMatrixDSym
#include <TMatrixTSym.h>													// for TMatrixTSym
#include <TMatrixTUtils.h>												// for TMatrixTRow

//#include <HFJetTruthGeneration/HFJetDefs.h>
#include <hfjettruthgeneration/HFJetDefs.h>

#include <iostream>

using namespace std;

//----------------------------------------------------------------------------//
//-- Constructor:
//--  simple initialization
//----------------------------------------------------------------------------//
AnaSvtxVertex::AnaSvtxVertex(const string &name):
  SubsysReco( name ),
  _events( NULL )
{
  //initialize
  _event = 0;
  _outfile = "AnaSvtxVertex.root";

	_truth_container = NULL;
	_hepmc_event_map = NULL;
	_vtxmap_rave = NULL;
	_vtxmap_rave_refit = NULL;
	_vtxmap_acts = NULL;
	_vtxmap_acts_refit = NULL;
	_svtxevalstack = NULL;

	_trkmap = NULL;

	h2d_eta_pt_p_mvtx0 = NULL;
	h2d_eta_pt_p_mvtx1 = NULL;
	h2d_eta_pt_p_mvtx2 = NULL;
	h2d_eta_pt_p_mvtx3 = NULL;

	h2d_eta_pt_s_mvtx0 = NULL;
	h2d_eta_pt_s_mvtx1 = NULL;
	h2d_eta_pt_s_mvtx2 = NULL;
	h2d_eta_pt_s_mvtx3 = NULL;
}

//----------------------------------------------------------------------------//
//-- Init():
//--   Intialize all histograms, trees, and ntuples
//----------------------------------------------------------------------------//
int AnaSvtxVertex::Init(PHCompositeNode *topNode)
{
  cout << PHWHERE << " Openning file " << _outfile << endl;
  PHTFileServer::get().open( _outfile, "RECREATE");

  // create TTree
  _events = new TTree("events", "Svtx Event");
  _events->Branch("event",0,"event/I");

	_events->Branch("vtx_gen",0,"vtx_gen[3]/F");

	_events->Branch("ntrack_reco",0,"ntrack_reco/I");
	_events->Branch("ntrack_reco_mvtx",0,"ntrack_reco_mvtx/I");

	_events->Branch("nvertex_reco",0,"nvertex_reco/I");
	_events->Branch("vtx_reco",0,"vtx_reco[nvertex_reco][3]/F");
	_events->Branch("vtx_reco_err",0,"vtx_reco_err[nvertex_reco][3]/F");
	_events->Branch("vtx_reco_ntrack",0,"vtx_reco_ntrack[nvertex_reco]/I");

	_events->Branch("nvertex_rave",0,"nvertex_rave/I");
	_events->Branch("vtx_rave",0,"vtx_rave[nvertex_rave][3]/F");
	_events->Branch("vtx_rave_err",0,"vtx_rave_err[nvertex_rave][3]/F");
	_events->Branch("vtx_rave_ntrack",0,"vtx_rave_ntrack[nvertex_rave]/I");

	_events->Branch("nvertex_rave_refit",0,"nvertex_rave_refit/I");
	_events->Branch("vtx_rave_refit",0,"vtx_rave_refit[nvertex_rave_refit][3]/F");
	_events->Branch("vtx_rave_refit_err",0,"vtx_rave_refit_err[nvertex_rave_refit][3]/F");
	_events->Branch("vtx_rave_refit_ntrack",0,"vtx_rave_refit_ntrack[nvertex_rave_refit]/I");

	_events->Branch("nvertex_acts",0,"nvertex_acts/I");
	_events->Branch("vtx_acts",0,"vtx_acts[nvertex_acts][3]/F");
	_events->Branch("vtx_acts_err",0,"vtx_acts_err[nvertex_acts][3]/F");
	_events->Branch("vtx_acts_ntrack",0,"vtx_acts_ntrack[nvertex_acts]/I");

	_events->Branch("nvertex_acts_refit",0,"nvertex_acts_refit/I");
	_events->Branch("vtx_acts_refit",0,"vtx_acts_refit[nvertex_acts_refit][3]/F");
	_events->Branch("vtx_acts_refit_err",0,"vtx_acts_refit_err[nvertex_acts_refit][3]/F");
	_events->Branch("vtx_acts_refit_ntrack",0,"vtx_acts_refit_ntrack[nvertex_acts_refit]/I");

	_events->Branch("njet04_true",0,"njet04_true/I");
	_events->Branch("jet04_prop_parton",0,"jet04_prop_parton[njet04_true]/I");
	_events->Branch("jet04_prop_hadron",0,"jet04_prop_hadron[njet04_true]/I");
	_events->Branch("jet04_pt",0,"jet04_pt[njet04_true]/F");
	_events->Branch("jet04_eta",0,"jet04_eta[njet04_true]/F");
	_events->Branch("jet04_mass",0,"jet04_mass[njet04_true]/F");

	_events->Branch("rave_sv_pT10_nvtx",0,"rave_sv_pT10_nvtx[njet04_true]/I");
	_events->Branch("rave_sv_pT10_vtx_x",0,"rave_sv_pT10_vtx_x[njet04_true][30]/F");
	_events->Branch("rave_sv_pT10_vtx_y",0,"rave_sv_pT10_vtx_y[njet04_true][30]/F");
	_events->Branch("rave_sv_pT10_vtx_z",0,"rave_sv_pT10_vtx_z[njet04_true][30]/F");
	_events->Branch("rave_sv_pT10_vtx_pt",0,"rave_sv_pT10_vtx_pt[njet04_true][30]/F");
	_events->Branch("rave_sv_pT10_vtx_mass",0,"rave_sv_pT10_vtx_mass[njet04_true][30]/F");
	_events->Branch("rave_sv_pT10_vtx_mass_corr",0,"rave_sv_pT10_vtx_mass_corr[njet04_true][30]/F");
	_events->Branch("rave_sv_pT10_vtx_ntrack",0,"rave_sv_pT10_vtx_ntrack[njet04_true][30]/I");

	h2d_eta_pt_p_mvtx0 = new TH2D("h2d_eta_pt_p_mvtx0","h2d_eta_pt_p_mvtx0",30,-1.5,1.5,100,0,10);
	h2d_eta_pt_p_mvtx1 = new TH2D("h2d_eta_pt_p_mvtx1","h2d_eta_pt_p_mvtx1",30,-1.5,1.5,100,0,10);
	h2d_eta_pt_p_mvtx2 = new TH2D("h2d_eta_pt_p_mvtx2","h2d_eta_pt_p_mvtx2",30,-1.5,1.5,100,0,10);
	h2d_eta_pt_p_mvtx3 = new TH2D("h2d_eta_pt_p_mvtx3","h2d_eta_pt_p_mvtx3",30,-1.5,1.5,100,0,10);

	h2d_eta_pt_s_mvtx0 = new TH2D("h2d_eta_pt_s_mvtx0","h2d_eta_pt_s_mvtx0",30,-1.5,1.5,100,0,10);
	h2d_eta_pt_s_mvtx1 = new TH2D("h2d_eta_pt_s_mvtx1","h2d_eta_pt_s_mvtx1",30,-1.5,1.5,100,0,10);
	h2d_eta_pt_s_mvtx2 = new TH2D("h2d_eta_pt_s_mvtx2","h2d_eta_pt_s_mvtx2",30,-1.5,1.5,100,0,10);
	h2d_eta_pt_s_mvtx3 = new TH2D("h2d_eta_pt_s_mvtx3","h2d_eta_pt_s_mvtx3",30,-1.5,1.5,100,0,10);

  return 0;
}

int AnaSvtxVertex::InitRun(PHCompositeNode *topNode)
{

	TGeoManager *tgeo_manager = PHGeomUtility::GetTGeoManager(topNode);
	PHField *field = PHFieldUtility::GetFieldMapNode(nullptr, topNode);

	_fitter = PHGenFit::Fitter::getInstance(tgeo_manager, field, "DafRef", "RKTrackRep", false);
	
	if (!_fitter)
	{
		cerr << PHWHERE << endl;
		return Fun4AllReturnCodes::ABORTRUN;
	}

	_vertex_finder = new genfit::GFRaveVertexFactory();
	_vertex_finder->setMethod("avr-smoothing:1");

	return 0;

}

//----------------------------------------------------------------------------//
//-- process_event():
//--   Call user instructions for every event.
//--   This function contains the analysis structure.
//----------------------------------------------------------------------------//
int AnaSvtxVertex::process_event(PHCompositeNode *topNode)
{
  _event++;
  if (_event % 1000 == 0)
    cout << PHWHERE << "Events processed: " << _event << endl;

  GetNodes(topNode);

  if (!_svtxevalstack) {
    _svtxevalstack = new SvtxEvalStack(topNode);
    _svtxevalstack->set_strict(false);
    //_svtxevalstack->set_verbosity(verbosity + 1);
    _svtxevalstack->next_event(topNode);
  } else {
    _svtxevalstack->next_event(topNode);
  }

	fill_tree(topNode);

  return 0;
}

//----------------------------------------------------------------------------//
//-- End():
//--   End method, wrap everything up
//----------------------------------------------------------------------------//
int AnaSvtxVertex::End(PHCompositeNode *topNode)
{

	cout << "-----AnaSvtxVertex::End------" << endl;

  PHTFileServer::get().cd( _outfile );
  _events->Write();
	h2d_eta_pt_p_mvtx0->Write();
	h2d_eta_pt_p_mvtx1->Write();
	h2d_eta_pt_p_mvtx2->Write();
	h2d_eta_pt_p_mvtx3->Write();

	h2d_eta_pt_s_mvtx0->Write();
	h2d_eta_pt_s_mvtx1->Write();
	h2d_eta_pt_s_mvtx2->Write();
	h2d_eta_pt_s_mvtx3->Write();
  //PHTFileServer::get().close();

	//if ( trutheval ) trutheval->Delete();
	//if ( trackeval ) trackeval->Delete();

  return 0;
}


//----------------------------------------------------------------------------//
//-- fill_tree():
//--   Fill the trees with truth, track fit, and cluster information
//----------------------------------------------------------------------------//
void AnaSvtxVertex::fill_tree(PHCompositeNode *topNode)
{

  // Make sure to reset all the TTree variables before trying to set them.
  reset_variables();

	SvtxTrackEval *trackeval = _svtxevalstack->get_track_eval();
	SvtxTruthEval *trutheval = _svtxevalstack->get_truth_eval();

	int nvertex_reco = 0;
	int nvertex_rave = 0;
	int nvertex_rave_refit = 0;
	int nvertex_acts = 0;
	int nvertex_acts_refit = 0;

	int ntrack_reco = 0;
	int ntrack_reco_mvtx = 0;

	int njet04_true = 0;

	//cout << "-EVENT: " << _event << "---------------------------------------------------------------------------------------------------" << endl;

	if ( _hepmc_event_map ){
	//if ( 0 ){

		for (PHHepMCGenEventMap::ConstIter iter=_hepmc_event_map->begin(); iter!=_hepmc_event_map->end(); ++iter){
			const PHHepMCGenEvent *hepmc_event = iter->second;

			float xx = hepmc_event->get_collision_vertex().x();
			float yy = hepmc_event->get_collision_vertex().y();
			float zz = hepmc_event->get_collision_vertex().z();

			vtx_gen[0] = xx;
			vtx_gen[1] = yy;
			vtx_gen[2] = zz;

		}
	}//hepmc_event_map


	if ( _vtxmap ){

		for (SvtxVertexMap::ConstIter iter=_vtxmap->begin(); iter!=_vtxmap->end(); ++iter){
			SvtxVertex *vtx = iter->second;

			vtx_reco[nvertex_reco][0] = vtx->get_x();
			vtx_reco[nvertex_reco][1] = vtx->get_y();
			vtx_reco[nvertex_reco][2] = vtx->get_z();

			vtx_reco_err[nvertex_reco][0] = vtx->get_error(0,0);
			vtx_reco_err[nvertex_reco][1] = vtx->get_error(1,1);
			vtx_reco_err[nvertex_reco][2] = vtx->get_error(2,2);

			vtx_reco_ntrack[nvertex_reco] = vtx->size_tracks();

			nvertex_reco++;
		}
	}//vtxmap


	if ( _vtxmap_rave ){

		for (SvtxVertexMap::ConstIter iter=_vtxmap_rave->begin(); iter!=_vtxmap_rave->end(); ++iter){

			SvtxVertex *vtx = iter->second;

			vtx_rave[nvertex_rave][0] = vtx->get_x();
			vtx_rave[nvertex_rave][1] = vtx->get_y();
			vtx_rave[nvertex_rave][2] = vtx->get_z();

			vtx_rave_err[nvertex_rave][0] = vtx->get_error(0,0);
			vtx_rave_err[nvertex_rave][1] = vtx->get_error(1,1);
			vtx_rave_err[nvertex_rave][2] = vtx->get_error(2,2);

			vtx_rave_ntrack[nvertex_rave] = vtx->size_tracks();

			nvertex_rave++;
		}//vtx_iter
	}//vtxmap_rave

	/*
	if ( _vtxmap_rave_refit ){

		for (SvtxVertexMap::ConstIter iter=_vtxmap_rave_refit->begin(); iter!=_vtxmap_rave_refit->end(); ++iter){

			SvtxVertex *vtx = iter->second;

			vtx_rave_refit[nvertex_rave_refit][0] = vtx->get_x();
			vtx_rave_refit[nvertex_rave_refit][1] = vtx->get_y();
			vtx_rave_refit[nvertex_rave_refit][2] = vtx->get_z();

			vtx_rave_refit_err[nvertex_rave_refit][0] = vtx->get_error(0,0);
			vtx_rave_refit_err[nvertex_rave_refit][1] = vtx->get_error(1,1);
			vtx_rave_refit_err[nvertex_rave_refit][2] = vtx->get_error(2,2);

			vtx_rave_refit_ntrack[nvertex_rave_refit] = vtx->size_tracks();

			nvertex_rave_refit++;
		}//vtx_iter
	}//vtxmap_rave_refit
	*/

	if ( _vtxmap_acts ){

		for (SvtxVertexMap::ConstIter iter=_vtxmap_acts->begin(); iter!=_vtxmap_acts->end(); ++iter){

			SvtxVertex *vtx = iter->second;

			vtx_acts[nvertex_acts][0] = vtx->get_x();
			vtx_acts[nvertex_acts][1] = vtx->get_y();
			vtx_acts[nvertex_acts][2] = vtx->get_z();

			vtx_acts_err[nvertex_acts][0] = vtx->get_error(0,0);
			vtx_acts_err[nvertex_acts][1] = vtx->get_error(1,1);
			vtx_acts_err[nvertex_acts][2] = vtx->get_error(2,2);

			vtx_acts_ntrack[nvertex_acts] = vtx->size_tracks();

			nvertex_acts++;
		}//vtx_iter
	}//vtxmap_acts

	if ( _vtxmap_acts_refit ){

		for (SvtxVertexMap::ConstIter iter=_vtxmap_acts_refit->begin(); iter!=_vtxmap_acts_refit->end(); ++iter){

			SvtxVertex *vtx = iter->second;

			vtx_acts_refit[nvertex_acts_refit][0] = vtx->get_x();
			vtx_acts_refit[nvertex_acts_refit][1] = vtx->get_y();
			vtx_acts_refit[nvertex_acts_refit][2] = vtx->get_z();

			vtx_acts_refit_err[nvertex_acts_refit][0] = vtx->get_error(0,0);
			vtx_acts_refit_err[nvertex_acts_refit][1] = vtx->get_error(1,1);
			vtx_acts_refit_err[nvertex_acts_refit][2] = vtx->get_error(2,2);

			vtx_acts_refit_ntrack[nvertex_acts_refit] = vtx->size_tracks();

			nvertex_acts_refit++;
		}//vtx_iter
	}//vtxmap_acts_refit

	//Primary vertex reconstruction
	_vertex_finder->setMethod("avf-smoothing:1");
	vector<genfit::Track*> gf_tracks;
	vector<genfit::GFRaveVertex*> rave_vertices;

	if ( _trkmap ){
		for (SvtxTrackMap::ConstIter iter=_trkmap->begin(); iter!=_trkmap->end();++iter){
			SvtxTrack *svtx_track = iter->second;

			//Truth information
			PHG4Particle *g4particle = trackeval->max_truth_particle_by_nclusters(svtx_track);
			if ( !g4particle ) continue;

			PHG4VtxPoint *vtx = trutheval->get_vertex(g4particle);

			float gvx = vtx->get_x();
			float gvy = vtx->get_y();
			float gvz = vtx->get_z();

			float dist = sqrt(pow(gvx-vtx_gen[0],2) + pow(gvy-vtx_gen[1],2) + pow(gvz-vtx_gen[2],2));
			if ( dist > 2 ) continue;

			int nmvtx = 0;

			for (SvtxTrack::ConstClusterKeyIter iter = svtx_track->begin_cluster_keys(); 
					iter != svtx_track->end_cluster_keys();
					++iter)
			{
				TrkrDefs::cluskey cluster_key = *iter;
				unsigned int layer = TrkrDefs::getLayer(cluster_key);

				if ( layer<3 ){
					++nmvtx;
				}
			}

			float px = svtx_track->get_px();
			float py = svtx_track->get_py();
			float pt = sqrt(px*px + py*py);
			float eta = svtx_track->get_eta();

			if ( dist < 1e-3 )
			{
				if ( nmvtx==0 ){
					h2d_eta_pt_p_mvtx0->Fill(eta, pt);
				}else if ( nmvtx==1 ){
					h2d_eta_pt_p_mvtx1->Fill(eta, pt);
				}else if ( nmvtx==2 ){
					h2d_eta_pt_p_mvtx2->Fill(eta, pt);
				}else {
					h2d_eta_pt_p_mvtx3->Fill(eta, pt);
				}
			}//prompt particale(primary)
			else
			{
				if ( nmvtx==0 ){
					h2d_eta_pt_s_mvtx0->Fill(eta, pt);
				}else if ( nmvtx==1 ){
					h2d_eta_pt_s_mvtx1->Fill(eta, pt);
				}else if ( nmvtx==2 ){
					h2d_eta_pt_s_mvtx2->Fill(eta, pt);
				}else{
					h2d_eta_pt_s_mvtx3->Fill(eta, pt);
				}
			}//non-prompt particle(secondary)

			ntrack_reco++;
			if ( nmvtx>2 ) ntrack_reco_mvtx++;

		}//iter

#if 1//Primary vertex reconstruction
		for (SvtxTrackMap::ConstIter iter=_trkmap->begin(); iter!=_trkmap->end(); ++iter)
		{
			SvtxTrack *svtx_track = iter->second;

			/*
			//Truth information
			PHG4Particle *g4particle = trackeval->max_truth_particle_by_nclusters(svtx_track);
			if ( !g4particle ) continue;

			PHG4VtxPoint *vtx = trutheval->get_vertex(g4particle);

			float gvx = vtx->get_x();
			float gvy = vtx->get_y();
			float gvz = vtx->get_z();

			float dist = sqrt(pow(gvx-vtx_gen[0],2) + pow(gvy-vtx_gen[1],2) + pow(gvz-vtx_gen[2],2));
			if ( dist > 2 ) continue;
			*/

			int nmvtx = 0;

			for (SvtxTrack::ConstClusterKeyIter iter = svtx_track->begin_cluster_keys(); 
					iter != svtx_track->end_cluster_keys(); 
					++iter)
			{
				TrkrDefs::cluskey cluster_key = *iter;
				unsigned int layer = TrkrDefs::getLayer(cluster_key);

				if ( layer<3 ){
					++nmvtx;
				}
			}
	
			if ( nmvtx<3 ) continue;

			auto genfit_track = TranslateSvtxToGenFitTrack(svtx_track);
			if (!genfit_track) continue;
			gf_tracks.push_back(const_cast<genfit::Track*>(genfit_track));

		}
		
		if (gf_tracks.size() >= 2)
		{
			try
			{
				_vertex_finder->findVertices(&rave_vertices, gf_tracks);
			}
			catch (...)
			{
				if (Verbosity() >1)
					std::cout << PHWHERE << "GFRaveVertexFactory::findVertices failed!";
			}
		}

		for (genfit::GFRaveVertex* rave_vtx : rave_vertices)
		{
			if (!rave_vtx) continue;
	
			vtx_rave_refit[nvertex_rave_refit][0]=rave_vtx->getPos().X(); 
			vtx_rave_refit[nvertex_rave_refit][1]=rave_vtx->getPos().Y(); 
			vtx_rave_refit[nvertex_rave_refit][2]=rave_vtx->getPos().Z();

			vtx_rave_refit_err[nvertex_rave_refit][0]=rave_vtx->getCov()[0][0]; 
			vtx_rave_refit_err[nvertex_rave_refit][1]=rave_vtx->getCov()[1][1]; 
			vtx_rave_refit_err[nvertex_rave_refit][2]=rave_vtx->getCov()[2][2];

			vtx_rave_refit_ntrack[nvertex_rave_refit] = rave_vtx->getNTracks();

			nvertex_rave_refit++;
		}
#endif//Primary vertex reconstruction

	}//trkmap


	//Secondary vertex reconstruction
	_vertex_finder->setMethod("avr-smoothing:1");
	gf_tracks.clear();
	rave_vertices.clear();

	if ( _jetmap04 ){
		for (JetMap::ConstIter iter = _jetmap04->begin(); iter!=_jetmap04->end();++iter)
		{
			Jet *jet_true = iter->second;

			if ( !jet_true ) continue;
			//else if(jet_true->get_pt() <= 10) continue;
			else if(jet_true->get_pt() <= 15) continue;
			if ( jet_true->has_property(static_cast<Jet::PROPERTY>(prop_JetPartonFlavor)) )
			{
				jet04_prop_parton[njet04_true] = int(jet_true->get_property(static_cast<Jet::PROPERTY>(prop_JetPartonFlavor)));
				jet04_prop_hadron[njet04_true] = int(jet_true->get_property(static_cast<Jet::PROPERTY>(prop_JetHadronFlavor)));
			}//jet property

			jet04_pt[njet04_true] = jet_true->get_pt();
			jet04_eta[njet04_true] = jet_true->get_eta();
			jet04_mass[njet04_true] = jet_true->get_mass();

			float jet_eta = jet_true->get_eta();
			float jet_phi = jet_true->get_phi();

			if (fabs(jet_eta) > 0.6) continue;

			vector<genfit::Track*> gf_tracks;

			if ( _trkmap )
			{
				for (SvtxTrackMap::ConstIter iter=_trkmap->begin(); iter!=_trkmap->end(); ++iter)
				{
					SvtxTrack *svtx_track = iter->second;

					/*
					//Truth information
					PHG4Particle *g4particle = trackeval->max_truth_particle_by_nclusters(svtx_track);
					if ( !g4particle ) continue;

					PHG4VtxPoint *vtx = trutheval->get_vertex(g4particle);

					float gvx = vtx->get_x();
					float gvy = vtx->get_y();
					float gvz = vtx->get_z();

					float dist = sqrt(pow(gvx-vtx_gen[0],2) + pow(gvy-vtx_gen[1],2) + pow(gvz-vtx_gen[2],2));
					if ( dist > 2 ) continue;
					*/

					int nmvtx=0;

					for (SvtxTrack::ConstClusterKeyIter iter = svtx_track->begin_cluster_keys(); 
							iter != svtx_track->end_cluster_keys(); 
							++iter)
					{
						TrkrDefs::cluskey cluster_key = *iter;
						unsigned int layer = TrkrDefs::getLayer(cluster_key);
				
						if ( layer<3 )
						{
							++nmvtx;
						}
					}

					float track_eta = svtx_track->get_eta();
					float track_phi = svtx_track->get_phi();


					float deta = jet_eta-track_eta;
					float dphi = jet_phi-track_phi;
					dphi = atan2(sin(dphi),cos(dphi));
					
					float dR = sqrt(pow(dphi,2)+pow(deta,2));

					if ( nmvtx<3 ) continue;

					auto genfit_track = TranslateSvtxToGenFitTrack(svtx_track);
					if ( !genfit_track ) continue;
					else if ( dR <= 0.4 )
					{
						gf_tracks.push_back(const_cast<genfit::Track*>(genfit_track));
					}

				}//iter_track

				vector<genfit::GFRaveVertex*> rave_vertices;
			
				if(gf_tracks.size() >= 2)
				{
					try
					{
						_vertex_finder->findVertices(&rave_vertices,gf_tracks);
					}
					catch (...)
					{
						if (Verbosity() >1)
							std::cout << PHWHERE << "GFRaveVertexFactory::findVertices failed!";
					}
					
					for (genfit::GFRaveVertex * rave_vtx : rave_vertices)
					{
						int isv = rave_sv_pT10_nvtx[njet04_true];
						rave_sv_pT10_vtx_x[njet04_true][isv]=rave_vtx->getPos().X();
						rave_sv_pT10_vtx_y[njet04_true][isv]=rave_vtx->getPos().Y();
						rave_sv_pT10_vtx_z[njet04_true][isv]=rave_vtx->getPos().Z();

						float sum_E = 0, sum_px = 0, sum_py = 0, sum_pz = 0;
						int sv_ntrk = 0;

						for (unsigned int itrk=0;itrk < rave_vtx->getNTracks();itrk++)
						{
							float weight_trk = rave_vtx->getParameters(itrk)->getWeight();
							if (weight_trk < 0.6) continue;

							TVector3 mom_trk = rave_vtx->getParameters(itrk)->getMom();
							sum_px += mom_trk.X();
							sum_py += mom_trk.Y();
							sum_pz += mom_trk.Z();
							sum_E += sqrt(mom_trk.Mag2() + 0.140*0.140);

							sv_ntrk++;
						}//itrk

						/*
						float vtx_mass, vtx_px, vtx_py, vtx_pz;
						int ntrk_good_pv = 0;
						GetSVMass_mom(rave_vtx,vtx_mass,vtx_px_vtx_py_vtx_pz_ntrk_good_pv);
						*/

						float sv_pt = sqrt(sum_px*sum_px + sum_py*sum_py);
						float sv_mass = sqrt(sum_E*sum_E - sum_px*sum_px - sum_py*sum_py - sum_pz*sum_pz);

						TVector3 vec1(sum_px, sum_py, sum_pz);
						TVector3 vec2(rave_sv_pT10_vtx_x[njet04_true][isv]-vtx_rave_refit[0][0],rave_sv_pT10_vtx_y[njet04_true][isv]-vtx_rave_refit[0][1],rave_sv_pT10_vtx_z[njet04_true][isv]-vtx_rave_refit[0][2]);
						float theta = vec1.Angle(vec2);
						float A = vec1.Mag()*sin(theta);
						float sv_mass_corr = sqrt(sv_mass*sv_mass + A*A) + A;

						rave_sv_pT10_vtx_pt[njet04_true][isv] = sv_pt;
						rave_sv_pT10_vtx_mass[njet04_true][isv] = sv_mass;
						rave_sv_pT10_vtx_mass_corr[njet04_true][isv] = sv_mass_corr;
						rave_sv_pT10_vtx_ntrack[njet04_true][isv] = sv_ntrk;

						rave_sv_pT10_nvtx[njet04_true]++;

						if (rave_sv_pT10_nvtx[njet04_true] >= 30) break;
					}//iter_vertex
				}//gf_tracks.size() >= 2

			}//trkmap

			njet04_true++;

			if ( njet04_true >= 100 ) break;
		}//iter_jet04
	}//jetmap04

	int count = 0;
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(&_event);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(vtx_gen);

	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(&ntrack_reco);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(&ntrack_reco_mvtx);

	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(&nvertex_reco);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(vtx_reco);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(vtx_reco_err);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(vtx_reco_ntrack);

	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(&nvertex_rave);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(vtx_rave);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(vtx_rave_err);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(vtx_rave_ntrack);

	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(&nvertex_rave_refit);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(vtx_rave_refit);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(vtx_rave_refit_err);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(vtx_rave_refit_ntrack);

	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(&nvertex_acts);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(vtx_acts);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(vtx_acts_err);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(vtx_acts_ntrack);

	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(&nvertex_acts_refit);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(vtx_acts_refit);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(vtx_acts_refit_err);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(vtx_acts_refit_ntrack);

	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(&njet04_true);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(jet04_prop_parton);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(jet04_prop_hadron);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(jet04_pt);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(jet04_eta);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(jet04_mass);

	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(rave_sv_pT10_nvtx);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(rave_sv_pT10_vtx_x);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(rave_sv_pT10_vtx_y);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(rave_sv_pT10_vtx_z);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(rave_sv_pT10_vtx_pt);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(rave_sv_pT10_vtx_mass);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(rave_sv_pT10_vtx_mass_corr);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(rave_sv_pT10_vtx_ntrack);

	_events->Fill();

  return;

}

//----------------------------------------------------------------------------//
//-- reset_variables():
//--   Reset all the tree variables to their default values.
//--   Needs to be called at the start of every event
//----------------------------------------------------------------------------//
void AnaSvtxVertex::reset_variables()
{

	vtx_gen[0] = vtx_gen[1] = vtx_gen[2] = -999;

	for (int ii=0; ii<100; ii++){
		vtx_reco[ii][0] = vtx_reco[ii][1] = vtx_reco[ii][2] = -999;
		vtx_reco_err[ii][0] = vtx_reco_err[ii][1] = vtx_reco_err[ii][2] = -999;

		vtx_rave[ii][0] = vtx_rave[ii][1] = vtx_rave[ii][2] = -999;
		vtx_rave_err[ii][0] = vtx_rave_err[ii][1] = vtx_rave_err[ii][2] = -999;

		vtx_acts[ii][0] = vtx_acts[ii][1] = vtx_acts[ii][2] = -999;
		vtx_acts_err[ii][0] = vtx_acts_err[ii][1] = vtx_acts_err[ii][2] = -999;

		vtx_rave_refit[ii][0] = vtx_rave_refit[ii][1] = vtx_rave_refit[ii][2] = -999;
		vtx_rave_refit_err[ii][0] = vtx_rave_refit_err[ii][1] = vtx_rave_refit_err[ii][2] = -999;

		vtx_acts_refit[ii][0] = vtx_acts_refit[ii][1] = vtx_acts_refit[ii][2] = -999;
		vtx_acts_refit_err[ii][0] = vtx_acts_refit_err[ii][1] = vtx_acts_refit_err[ii][2] = -999;

		vtx_reco_ntrack[ii] = vtx_rave_ntrack[ii] = vtx_rave_refit_ntrack[ii] = 0;
		vtx_acts_ntrack[ii] = vtx_acts_refit_ntrack[ii] = 0;

		jet04_pt[ii] = jet04_eta[ii] = jet04_mass[ii] = 0;
		jet04_prop_parton[ii] = jet04_prop_hadron[ii] = 0;

		rave_sv_pT10_nvtx[ii] = 0;
		for (int jj=0; jj<30; jj++)
		{
			rave_sv_pT10_vtx_x[ii][jj] = 0;
			rave_sv_pT10_vtx_y[ii][jj] = 0;
			rave_sv_pT10_vtx_z[ii][jj] = 0;

			rave_sv_pT10_vtx_pt[ii][jj] = 0;
			rave_sv_pT10_vtx_mass[ii][jj] = 0;
			rave_sv_pT10_vtx_mass_corr[ii][jj] = 0;
			rave_sv_pT10_vtx_ntrack[ii][jj] = 0;
		}//jj
	}//ii
}

//----------------------------------------------------------------------------//
//-- GetNodes():
//--   Get all the all the required nodes off the node tree
//----------------------------------------------------------------------------//
void AnaSvtxVertex::GetNodes(PHCompositeNode * topNode)
{
  //DST objects
	//
  //PHG4TruthInfoContainer
	/*
  _truth_container = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
  if (!_truth_container && _event<2)
  {
    cout << PHWHERE << " PHG4TruthInfoContainer node not found on node tree" << endl;
  }
	*/

	//HepMCGenEventMap
	_hepmc_event_map = findNode::getClass<PHHepMCGenEventMap>(topNode,"PHHepMCGenEventMap");
	if (!_hepmc_event_map && _event<2)
	{
		cout << PHWHERE << " PHHepMCGenEventMap node not found on node tree" << endl;
	}

  //SvtxVertexMap
  _vtxmap = findNode::getClass<SvtxVertexMap>(topNode,"SvtxVertexMap");
  if (!_vtxmap && _event<2)
  {
    cout << PHWHERE << " SvtxVertexMap node not found on node tree" << endl;
  }

  //SvtxVertexMap
  _vtxmap_rave = findNode::getClass<SvtxVertexMap>(topNode,"SvtxVertexMapRave");
  if (!_vtxmap_rave && _event<2)
  {
    cout << PHWHERE << " SvtxVertexMapRave node not found on node tree" << endl;
  }

  //SvtxVertexMap
  _vtxmap_rave_refit = findNode::getClass<SvtxVertexMap>(topNode,"SvtxVertexMapRaveRefitMvtx");
  if (!_vtxmap_rave_refit && _event<2)
  {
    cout << PHWHERE << " SvtxVertexMapRaveRefitMvtx node not found on node tree" << endl;
  }

  //SvtxVertexMap
  _vtxmap_acts = findNode::getClass<SvtxVertexMap>(topNode,"SvtxVertexMapActsRefit");
  if (!_vtxmap_acts && _event<2)
  {
    cout << PHWHERE << " SvtxVertexMapActs node not found on node tree" << endl;
  }

  //SvtxVertexMap
  _vtxmap_acts_refit = findNode::getClass<SvtxVertexMap>(topNode,"SvtxVertexMapActsRefitMvtx");
  if (!_vtxmap_acts_refit && _event<2)
  {
    cout << PHWHERE << " SvtxVertexMapActsRefitMvtx node not found on node tree" << endl;
  }

  //SvtxTrackMap
  _trkmap = findNode::getClass<SvtxTrackMap>(topNode,"SvtxTrackMap");
  if (!_trkmap && _event<2)
  {
    cout << PHWHERE << " SvtxTrackMap node not found on node tree" << endl;
  }

  //JetMap04
  _jetmap04 = findNode::getClass<JetMap>(topNode,"AntiKt_Truth_r04");
  if (!_jetmap04 && _event<2)
  {
    cout << PHWHERE << " AntiKt_Truth_r04 node not found on node tree" << endl;
  }

}

genfit::Track *AnaSvtxVertex::TranslateSvtxToGenFitTrack(SvtxTrack* svtx_track)
{
	try
	{
		// The first state is extracted to PCA, second one is the one with measurement
		SvtxTrackState *svtx_state(nullptr);
		//SvtxTrackState *svtx_state = (svtx_track->begin_states())->second;

		if (svtx_track->begin_states() == svtx_track->end_states())
		{
			cout << PHWHERE << "TranslateSvtxToGenFitTrack no state in track!" << endl;
			return nullptr;
		}
		else if (++(svtx_track->begin_states()) == svtx_track->end_states())
		{
			// only one state in track
			svtx_state = (svtx_track->begin_states())->second;
		}
		else
		{
			// multiple state in track
			// The first state is extracted to PCA, second one is the one with measurement
			svtx_state = (++(svtx_track->begin_states()))->second;
		}

		if (!svtx_state)
		{
			cout << PHWHERE << "TranslateSvtxToGenFitTrack invalid state found on track!" << endl;
			return nullptr;
		}

		TVector3 pos(svtx_state->get_x(), svtx_state->get_y(), svtx_state->get_z());
		TVector3 mom(svtx_state->get_px(), svtx_state->get_py(), svtx_state->get_pz());
		TMatrixDSym cov(6);
		for (int i = 0; i < 6; ++i)
		{
			for (int j = 0; j < 6; ++j)
			{
				cov[i][j] = svtx_state->get_error(i, j);
			}
		}

		genfit::AbsTrackRep *rep = new genfit::RKTrackRep(211);
		genfit::Track *genfit_track = new genfit::Track(rep, TVector3(0, 0, 0), TVector3(0, 0, 0));

		genfit::FitStatus *fs = new genfit::FitStatus();
		fs->setCharge(svtx_track->get_charge());
		fs->setChi2(svtx_track->get_chisq());
		fs->setNdf(svtx_track->get_ndf());
		fs->setIsFitted(true);
		fs->setIsFitConvergedFully(true);

		genfit_track->setFitStatus(fs, rep);

		genfit::TrackPoint *tp = new genfit::TrackPoint(genfit_track);

		genfit::KalmanFitterInfo *fi = new genfit::KalmanFitterInfo(tp, rep);
		tp->setFitterInfo(fi);

		genfit::MeasuredStateOnPlane *ms = new genfit::MeasuredStateOnPlane(rep);
		ms->setPosMomCov(pos, mom, cov);
		genfit::KalmanFittedStateOnPlane *kfs = new genfit::KalmanFittedStateOnPlane(*ms, 1., 1.);

		//< According to the special order of using the stored states
		fi->setForwardUpdate(kfs);

		genfit_track->insertPoint(tp);

		return genfit_track;
	}
	catch (...)
	{
		cout << PHWHERE << "TranslateSvtxToGenFitTrack failed!" << endl;
	}

	return nullptr;
}


/*
int SVReco::GetSVMass_mom(
		const genfit::GFRaveVertex* rave_vtx,
		float & vtx_mass,
		float & vtx_px,
		float & vtx_py,
		float & vtx_pz,
		int & ntrk_good_pv
		){

	float sum_E = 0, sum_px = 0, sum_py = 0, sum_pz = 0;

	int N_good = 0, N_good_pv = 0;

	for (unsigned int itrk=0; itrk<rave_vtx->getNTracks(); itrk++){
		TVector3 mom_trk = rave_vtx->getParameters(itrk)->getMom(); 

		double w_trk = rave_vtx->getParameters(itrk)->getWeight();

		sum_px += mom_trk.X();
		sum_py += mom_trk.Y();
		sum_pz += mom_trk.Z();
		sum_E += sqrt(mom_trk.Mag2() + 0.140*0.140);

		//cout << "W: " << w_trk << endl;
		if ( w_trk>0.7 ){
			N_good++;

			unsigned int rvtrk_mc_id = rave_vtx->getParameters(itrk)->getTrack()->getMcTrackId();
			//cout << "mc_id: " << rvtrk_mc_id << ", wt: " << svtxtrk_wt_map[rvtrk_mc_id] << endl;
			if ( svtxtrk_wt_map[rvtrk_mc_id]>0.7 ){
				N_good_pv++;
			}
		}//

	}//itrk

	vtx_mass =  sqrt(sum_E*sum_E - sum_px*sum_px - sum_py*sum_py - sum_pz*sum_pz);
	vtx_px = sum_px;
	vtx_py = sum_py;
	vtx_pz = sum_pz;

	ntrk_good_pv = N_good_pv;

	//cout << "Mass: " << vtx_mass << ", pT: " << vtx_pT << endl;
	return N_good;
}
*/
