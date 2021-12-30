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


#include "AnaSvtxTracksForGenFit.h"

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

#include <TTree.h>
#include <TVector3.h>

#include <iostream>

using namespace std;

//----------------------------------------------------------------------------//
//-- Constructor:
//--  simple initialization
//----------------------------------------------------------------------------//
AnaSvtxTracksForGenFit::AnaSvtxTracksForGenFit(const string &name):
  SubsysReco( name ),
  _flags( NONE ),
  _tracks( NULL ),
  _svtxevalstack( NULL )
{
  //initialize
  _event = 0;
  _outfile = "AnaSvtxTracksForGenFit.root";
}

//----------------------------------------------------------------------------//
//-- Init():
//--   Intialize all histograms, trees, and ntuples
//----------------------------------------------------------------------------//
int AnaSvtxTracksForGenFit::Init(PHCompositeNode *topNode)
{
  cout << PHWHERE << " Openning file " << _outfile << endl;
  PHTFileServer::get().open( _outfile, "RECREATE");


  // create TTree
  _tracks = new TTree("tracks", "Svtx Tracks");
  _tracks->Branch("event", &event, "event/I");
  _tracks->Branch("ntracks", &ntracks, "ntracks/I");
  _tracks->Branch("gtrackID", gtrackID, "gtrackID[ntracks]/I");
  _tracks->Branch("gflavor", gflavor, "gflavor[ntracks]/I");
  _tracks->Branch("gpx", gpx, "gpx[ntracks]/F");
  _tracks->Branch("gpy", gpy, "gpy[ntracks]/F");
  _tracks->Branch("gpz", gpz, "gpz[ntracks]/F");
  _tracks->Branch("gvx", gvx, "gvx[ntracks]/F");
  _tracks->Branch("gvy", gvy, "gvy[ntracks]/F");
  _tracks->Branch("gvz", gvz, "gvz[ntracks]/F");
  _tracks->Branch("trackID", trackID, "trackID[ntracks]/I");
  _tracks->Branch("charge", charge, "charge[ntracks]/I");
  _tracks->Branch("nhits", nhits, "nhits[ntracks]/I");
  _tracks->Branch("px", px, "px[ntracks]/F");
  _tracks->Branch("py", py, "py[ntracks]/F");
  _tracks->Branch("pz", pz, "pz[ntracks]/F");
  _tracks->Branch("dca2d", dca2d, "dca2d[ntracks]/F");
  _tracks->Branch("clusterID", clusterID, "clusterID[ntracks][7]/I");
  _tracks->Branch("layer", layer, "layer[ntracks][7]/I");
  _tracks->Branch("x", x, "x[ntracks][7]/F");
  _tracks->Branch("y", y, "y[ntracks][7]/F");
  _tracks->Branch("z", z, "z[ntracks][7]/F");
  _tracks->Branch("size_dphi", size_dphi, "size_dphi[ntracks][7]/F");
  _tracks->Branch("size_dz", size_dz, "size_dz[ntracks][7]/F");


  return 0;
}

//----------------------------------------------------------------------------//
//-- process_event():
//--   Call user instructions for every event.
//--   This function contains the analysis structure.
//----------------------------------------------------------------------------//
int AnaSvtxTracksForGenFit::process_event(PHCompositeNode *topNode)
{
  _event++;
  if (_event % 1000 == 0)
    cout << PHWHERE << "Events processed: " << _event << endl;

  GetNodes(topNode);

  if (!_svtxevalstack) {
    _svtxevalstack = new SvtxEvalStack(topNode);
    _svtxevalstack->set_strict(false);
    //_svtxevalstack->set_verbosity(verbosity + 1);
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
int AnaSvtxTracksForGenFit::End(PHCompositeNode *topNode)
{

  PHTFileServer::get().cd( _outfile );

  _tracks->Write();

  if (_svtxevalstack) delete _svtxevalstack;

  return 0;
}


//----------------------------------------------------------------------------//
//-- fill_tree():
//--   Fill the trees with truth, track fit, and cluster information
//----------------------------------------------------------------------------//
void AnaSvtxTracksForGenFit::fill_tree(PHCompositeNode *topNode)
{
  // Make sure to reset all the TTree variables before trying to set them.
  reset_variables();

  // get evaluators
  SvtxTrackEval*     trackeval = _svtxevalstack->get_track_eval();
  //SvtxClusterEval*   cluseval = _svtxevalstack->get_cluster_eval();
  SvtxTruthEval*     trutheval = _svtxevalstack->get_truth_eval();

	//set<SvtxCluster*> clus_set;

  if (_truth_container)
  {

    PHG4TruthInfoContainer::ConstRange range =  _truth_container->GetPrimaryParticleRange();

		ntracks = 0;
		for (PHG4TruthInfoContainer::ConstIterator iter=range.first; iter!=range.second; ++iter)
		{

			PHG4Particle* g4particle = iter->second;

			gtrackID[ntracks] = g4particle->get_track_id();
			gflavor[ntracks]  = g4particle->get_pid();

			if ( fabs(gflavor[ntracks])!=211 && fabs(gflavor[ntracks])!=321 && fabs(gflavor[ntracks])!=2212 ) continue;

			PHG4VtxPoint* vtx = trutheval->get_vertex(g4particle);
			gvx[ntracks]      = vtx->get_x();
			gvy[ntracks]      = vtx->get_y();
			gvz[ntracks]      = vtx->get_z();

			if ( gvz[ntracks]!=0 ) continue;

			gpx[ntracks]			 = g4particle->get_px();
			gpy[ntracks]			 = g4particle->get_py();
			gpz[ntracks]			 = g4particle->get_pz();

			TVector3 mom(gpx[ntracks], gpy[ntracks], gpz[ntracks]);
			if ( fabs(mom.Eta())>0.5 ) continue;

			/*
			clus_set.clear();
			clus_set = cluseval->all_clusters_from(g4particle);

			cout << "--------------------------------------------" << endl;
			cout << "SIZE : " << clus_set.size() << endl;

			for (set<SvtxCluster*>::iterator iter_clus=clus_set.begin(); iter_clus!=clus_set.end(); iter_clus++){
				unsigned int l = (*iter_clus)->get_layer();
				float xx = (*iter_clus)->get_x();
				float yy = (*iter_clus)->get_y();
				float zz = (*iter_clus)->get_z();
				//float ee = (*iter_clus)->get_e();
				int tmp_adc = (*iter_clus)->get_adc();

				if ( tmp_adc>adc[ntracks][l] ){ 
					layer[ntracks][l] = (int)l;
					x[ntracks][l] = xx;
					y[ntracks][l] = yy;
					z[ntracks][l] = zz;
					size_dphi[ntracks][l] = (*iter_clus)->get_phi_size();
					size_dz[ntracks][l] = (*iter_clus)->get_z_size();
				}//if
			}//for
			*/

			SvtxTrack* track = trackeval->best_track_from(g4particle);
			if (track)
			{
				trackID[ntracks]   = track->get_id();
				charge[ntracks]    = track->get_charge();
				nhits[ntracks]     = track->size_clusters();
				px[ntracks]        = track->get_px();
				py[ntracks]        = track->get_py();
				pz[ntracks]        = track->get_pz();
				dca2d[ntracks]     = track->get_dca2d();


				int iclus = 0;
				for (SvtxTrack::ConstClusterIter iter = track->begin_clusters();
						iter != track->end_clusters();
						++iter)
				{
					unsigned int cluster_id = *iter;
					SvtxCluster* cluster = _clustermap->get(cluster_id);
					unsigned int l = cluster->get_layer();

					//cout << l << endl;

					clusterID[ntracks][iclus] = (int)cluster_id;
					layer[ntracks][iclus] = (int)l;
					x[ntracks][iclus] = cluster->get_x();
					y[ntracks][iclus] = cluster->get_y();
					z[ntracks][iclus] = cluster->get_z();
					size_dphi[ntracks][iclus] = cluster->get_phi_size();
					size_dz[ntracks][iclus] = cluster->get_z_size();

					++iclus;
				}
				ntracks++;
			} // if(track)


			if ( ntracks>=100 ) break;
		} // for( iter)
	} //if (_truth_container)

  _tracks->Fill();
  return;

}

//----------------------------------------------------------------------------//
//-- reset_variables():
//--   Reset all the tree variables to their default values.
//--   Needs to be called at the start of every event
//----------------------------------------------------------------------------//
void AnaSvtxTracksForGenFit::reset_variables()
{
  event = -9999;

	for (int itrk=0; itrk<100; itrk++){
		//-- truth
		gtrackID[itrk] = -9999;
		gflavor[itrk] = -9999;
		gpx[itrk] = -9999;
		gpy[itrk] = -9999;
		gpz[itrk] = -9999;
		gvx[itrk] = -9999;
		gvy[itrk] = -9999;
		gvz[itrk] = -9999;

		//-- reco
		trackID[itrk] = -9999;
		charge[itrk] = -9999;
		nhits[itrk] = -9999;
		px[itrk] = -9999;
		py[itrk] = -9999;
		pz[itrk] = -9999;
		dca2d[itrk] = -9999;

		//-- clusters
		for (int i = 0; i < 7; i++)
		{
			clusterID[itrk][i] = -9999;
			layer[itrk][i] = -9999;
			x[itrk][i] = -9999;
			y[itrk][i] = -9999;
			z[itrk][i] = -9999;
			size_dphi[itrk][i] = -9999;
			size_dz[itrk][i] = -9999;
			adc[itrk][i] = -9999;
		}
	}

}

//----------------------------------------------------------------------------//
//-- GetNodes():
//--   Get all the all the required nodes off the node tree
//----------------------------------------------------------------------------//
void AnaSvtxTracksForGenFit::GetNodes(PHCompositeNode * topNode)
{
  //DST objects
  //Truth container
  _truth_container = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!_truth_container && _event < 2)
  {
    cout << PHWHERE
         << " PHG4TruthInfoContainer node not found on node tree"
         << endl;
  }

  //Svtx Clusters
  _clustermap = findNode::getClass<SvtxClusterMap>(topNode, "SvtxClusterMap");
  if (!_clustermap && _event < 2)
  {
    cout << PHWHERE
         << " SvtxClusterMap node not found on node tree"
         << endl;
  }


}




