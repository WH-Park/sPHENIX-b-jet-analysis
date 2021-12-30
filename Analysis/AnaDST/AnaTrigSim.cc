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


#include "AnaTrigSim.h"

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

#include <g4eval/CaloEvalStack.h>
#include <g4eval/CaloRawClusterEval.h>

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeomContainer.h>

#include <TTree.h>
#include <TVector3.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#include <HepMC/GenParticle.h>

#include <HFJetTruthGeneration/HFJetDefs.h>
#include <calotrigger/CaloTriggerInfo.h>

#include <iostream>

using namespace std;

//----------------------------------------------------------------------------//
//-- Constructor:
//--  simple initialization
//----------------------------------------------------------------------------//
AnaTrigSim::AnaTrigSim(const string &name):
  SubsysReco( name ),
  _T( NULL )
{
  //initialize
  _event = 0;
  _outfile = "AnaTrigSim.root";
	_triggerinfo = 0;
	_truthcontainer = 0;
	_cemc_clusters = 0;
	_caloevalstack = 0;

}

void AnaTrigSim::set_truncation( int emulate_truncation ) { 
	_emulate_truncation = emulate_truncation;
}

//----------------------------------------------------------------------------//
//-- Init():
//--   Intialize all histograms, trees, and ntuples
//----------------------------------------------------------------------------//
int AnaTrigSim::Init(PHCompositeNode *topNode)
{
  cout << PHWHERE << " Openning file " << _outfile << endl;
  PHTFileServer::get().open( _outfile, "RECREATE");

  // create TTree
	_T = new TTree("T", "Trig Sim");
  _T->Branch("event",0,"event/I");
  _T->Branch("npart",0,"npart/I");
  _T->Branch("part_pid",0,"part_pid[npart]/I");
  _T->Branch("part_px",0,"part_px[npart]/F");
  _T->Branch("part_py",0,"part_py[npart]/F");
  _T->Branch("part_pz",0,"part_pz[npart]/F");
  _T->Branch("part_eta",0,"part_eta[npart]/F");
  _T->Branch("trig_best_EMC_4x4",0,"trig_best_EMC_4x4[3]/F");
  _T->Branch("trig_best2_EMC_4x4",0,"trig_best2_EMC_4x4[3]/F");
  _T->Branch("nclus",0,"nclus/I");
  _T->Branch("clus_pid",0,"clus_pid[nclus]/I");
  _T->Branch("clus_ntower",0,"clus_ntower[nclus]/I");
  _T->Branch("clus_e",0,"clus_e[nclus]/F");


  return 0;
}

//----------------------------------------------------------------------------//
//-- process_event():
//--   Call user instructions for every event.
//--   This function contains the analysis structure.
//----------------------------------------------------------------------------//
int AnaTrigSim::process_event(PHCompositeNode *topNode)
{
  _event++;
  //if (1)
	if (_event % 10 == 0)
    cout << PHWHERE << "Events processed: " << _event << endl;

  GetNodes(topNode);

  fill_tree(topNode);

  return 0;
}

//----------------------------------------------------------------------------//
//-- End():
//--   End method, wrap everything up
//----------------------------------------------------------------------------//
int AnaTrigSim::End(PHCompositeNode *topNode)
{

	cout << "-----AnaTrigSim::End------" << endl;

  PHTFileServer::get().cd( _outfile );
  _T->Write();
  //PHTFileServer::get().close();

  return 0;
}


//----------------------------------------------------------------------------//
//-- fill_tree():
//--   Fill the trees with truth, track fit, and cluster information
//----------------------------------------------------------------------------//
void AnaTrigSim::fill_tree(PHCompositeNode *topNode)
{

  // Make sure to reset all the TTree variables before trying to set them.
  reset_variables();

	//cout << "reset_variables" << endl;

	//cout << "ready" << endl;

	cout << "-EVENT: " << _event << "---------------------------------------------------------------------------------------------------" << endl;
	
	if ( _triggerinfo ){

	}

	float pv_z = -9999;

	if ( _truthcontainer ){

		PHG4TruthInfoContainer::ConstVtxRange vtx_range =  _truthcontainer->GetPrimaryVtxRange();
		for (PHG4TruthInfoContainer::ConstVtxIterator iter=vtx_range.first; iter!=vtx_range.second; ++iter)
		{    
			PHG4VtxPoint *g4vertex = iter->second;
			pv_z = g4vertex->get_z();
			/*
			cout 
				<< "vertex x: " << g4vertex->get_x()
				<< "vertex y: " << g4vertex->get_y()
				<< "vertex z: " << g4vertex->get_z()
				<< endl;
			*/
		}
	}//_truthcontainer

	if ( fabs(pv_z)>10.0 ) return;

	int npart = 0;
	int part_pid[100] = {0};
	float part_px[100] = {0.}, part_py[100] = {0.}, part_pz[100] = {0.}, part_eta[100] = {0.};

	if ( _truthcontainer ){
		PHG4TruthInfoContainer::ConstRange range =  _truthcontainer->GetPrimaryParticleRange();
		for (PHG4TruthInfoContainer::ConstIterator iter=range.first; iter!=range.second; ++iter)
		{    

			PHG4Particle *g4particle = iter->second;
			if ( abs(g4particle->get_pid())!=11 ) continue;

			part_pid[npart] = g4particle->get_pid();
			part_px[npart] = g4particle->get_px();
			part_py[npart] = g4particle->get_py();
			part_pz[npart] = g4particle->get_pz();

			TVector3 vec(part_px[npart], part_py[npart], part_pz[npart]);
			part_eta[npart] = vec.Eta();

			npart++;
			/*
			cout 
				<< "pid: " << g4particle->get_pid()
				<< ", px: " << g4particle->get_px()
				<< ", py: " << g4particle->get_py()
				<< ", pz: " << g4particle->get_pz()
				<< endl;
			*/
		}
	}

	float trig_best_EMC_4x4[3] = {0};
	float trig_best2_EMC_4x4[3] = {0};

	if ( _triggerinfo ){
		trig_best_EMC_4x4[0] = _triggerinfo->get_best_EMCal_4x4_E();
		trig_best_EMC_4x4[1] = _triggerinfo->get_best_EMCal_4x4_eta();
		trig_best_EMC_4x4[2] = _triggerinfo->get_best_EMCal_4x4_phi();

		trig_best2_EMC_4x4[0] = _triggerinfo->get_best2_EMCal_4x4_E();
		trig_best2_EMC_4x4[1] = _triggerinfo->get_best2_EMCal_4x4_eta();
		trig_best2_EMC_4x4[2] = _triggerinfo->get_best2_EMCal_4x4_phi();
	}//_triggerinfo

	int nclus = 0;
	int clus_pid[100] = {0}, clus_ntower[100] = {0};
	float clus_e[100] = {0.};

	if ( _cemc_clusters && _caloevalstack ){

		CaloRawClusterEval* clustereval = _caloevalstack->get_rawcluster_eval();
		for (const auto & iterator : _cemc_clusters->getClustersMap()){

			RawCluster *cluster = iterator.second;
			PHG4Particle* primary = clustereval->max_truth_primary_particle_by_energy(cluster);

			clus_ntower[nclus] = cluster->getNTowers();
			clus_e[nclus] = cluster->get_energy(); 

			if ( clus_e[nclus]<0.2 ) continue;

			if ( primary ){
				clus_pid[nclus] = primary->get_pid();
			}else{
				clus_pid[nclus] = -999;
			}

			//cout << clus_pid[nclus] << ", " << clus_ntower[nclus] << ", " << clus_e[nclus] << endl;

			nclus++;
		}//iter
	}//_cemc_clusters

	int count = 0;

	((TBranch*) _T->GetListOfBranches()->At(count++))->SetAddress(&_event);
	((TBranch*) _T->GetListOfBranches()->At(count++))->SetAddress(&npart);
	((TBranch*) _T->GetListOfBranches()->At(count++))->SetAddress(part_pid);
	((TBranch*) _T->GetListOfBranches()->At(count++))->SetAddress(part_px);
	((TBranch*) _T->GetListOfBranches()->At(count++))->SetAddress(part_py);
	((TBranch*) _T->GetListOfBranches()->At(count++))->SetAddress(part_pz);
	((TBranch*) _T->GetListOfBranches()->At(count++))->SetAddress(part_eta);
	((TBranch*) _T->GetListOfBranches()->At(count++))->SetAddress(trig_best_EMC_4x4);
	((TBranch*) _T->GetListOfBranches()->At(count++))->SetAddress(trig_best2_EMC_4x4);
	((TBranch*) _T->GetListOfBranches()->At(count++))->SetAddress(&nclus);
	((TBranch*) _T->GetListOfBranches()->At(count++))->SetAddress(clus_pid);
	((TBranch*) _T->GetListOfBranches()->At(count++))->SetAddress(clus_ntower);
	((TBranch*) _T->GetListOfBranches()->At(count++))->SetAddress(clus_e);

	_T->Fill();

  return;

}

//----------------------------------------------------------------------------//
//-- reset_variables():
//--   Reset all the tree variables to their default values.
//--   Needs to be called at the start of every event
//----------------------------------------------------------------------------//
void AnaTrigSim::reset_variables()
{


}

//----------------------------------------------------------------------------//
//-- GetNodes():
//--   Get all the all the required nodes off the node tree
//----------------------------------------------------------------------------//
void AnaTrigSim::GetNodes(PHCompositeNode * topNode)
{
  //DST objects
	//
  //CaloTriggerInfo
  _triggerinfo = findNode::getClass<CaloTriggerInfo>(topNode, !_emulate_truncation ? "CaloTriggerInfo" : "CaloTriggerInfo_Truncate");
  if (!_triggerinfo && _event<2)
  {
    cout << PHWHERE << " CaloTriggerInfo node not found on node tree" << endl;
  }

  //TruthInfoTrigger
	_truthcontainer = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!_truthcontainer && _event<2)
  {
		cout << PHWHERE << " G4TruthInfo nodeContainer not found on node tree" << endl;
  }

  //RawClusterContainer
	_cemc_clusters = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_CEMC");
  if (!_cemc_clusters && _event<2)
  {
		cout << PHWHERE << " RawClusterContainer nodeContainer not found on node tree" << endl;
  }

	if (!_caloevalstack)
	{
		_caloevalstack = new CaloEvalStack(topNode, "CEMC");
		_caloevalstack->set_strict(false);
		_caloevalstack->set_verbosity(0);
	}
	else
	{
		_caloevalstack->next_event(topNode);
	}

}



