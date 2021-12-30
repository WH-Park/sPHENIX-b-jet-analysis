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


#include "AnaEventGen.h"

#include <phool/phool.h>
#include <phool/getClass.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4VtxPoint.h>
#include <fun4all/PHTFileServer.h>
#include <fun4all/Fun4AllServer.h>

#include <g4jets/JetMap.h>
#include <g4jets/Jet.h>

#include <TTree.h>
#include <TVector3.h>
#include <phhepmc/PHHepMCGenEventMap.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#include <HepMC/GenParticle.h>
#include <HepMC/HeavyIon.h>

//#include <HFJetTruthGeneration/HFJetDefs.h>

#include <iostream>

using namespace std;

//----------------------------------------------------------------------------//
//-- Constructor:
//--  simple initialization
//----------------------------------------------------------------------------//
AnaEventGen::AnaEventGen(const string &name):
  SubsysReco( name ),
  _events( NULL )
{
  //initialize
  _event = 0;
  _outfile = "AnaEventGen.root";

	_jet_pt_min = 10.0;

	_truth_container = NULL;
	_hepmc_event_map = NULL;

	_jetmap04 = NULL;
	_jetmap08 = NULL;

}

//----------------------------------------------------------------------------//
//-- Init():
//--   Intialize all histograms, trees, and ntuples
//----------------------------------------------------------------------------//
int AnaEventGen::Init(PHCompositeNode *topNode)
{
  cout << PHWHERE << " Openning file " << _outfile << endl;
  PHTFileServer::get().open( _outfile, "RECREATE");

  // create TTree
	_events = new TTree("T", "T");
	_events->Branch("event",0,"event/I");

	_events->Branch("ncoll",0,"ncoll/S");
	_events->Branch("ip",0,"ip/F");

	_events->Branch("npart",0,"npart/I");
	_events->Branch("part_pid",0,"part_pid[npart]/I");
	_events->Branch("part_px",0,"part_px[npart]/F");
	_events->Branch("part_py",0,"part_py[npart]/F");
	_events->Branch("part_pz",0,"part_pz[npart]/F");

	_events->Branch("njet04",0,"njet04/I");
	_events->Branch("jet04_pT",0,"jet04_pT[njet04]/F");
	_events->Branch("jet04_phi",0,"jet04_phi[njet04]/F");
	_events->Branch("jet04_eta",0,"jet04_eta[njet04]/F");

	/*
	_events->Branch("njet08",0,"njet08/I");
	_events->Branch("jet08_pT",0,"jet08_pT[njet08]/F");
	_events->Branch("jet08_phi",0,"jet08_phi[njet08]/F");
	_events->Branch("jet08_eta",0,"jet08_eta[njet08]/F");
	*/

  return 0;
}

//----------------------------------------------------------------------------//
//-- process_event():
//--   Call user instructions for every event.
//--   This function contains the analysis structure.
//----------------------------------------------------------------------------//
int AnaEventGen::process_event(PHCompositeNode *topNode)
{
  _event++;
  if (_event % 1000 == 0)
    cout << PHWHERE << "Events processed: " << _event << endl;

  GetNodes(topNode);

  fill_tree(topNode);

  return 0;
}

//----------------------------------------------------------------------------//
//-- End():
//--   End method, wrap everything up
//----------------------------------------------------------------------------//
int AnaEventGen::End(PHCompositeNode *topNode)
{

	cout << "-----AnaEventGen::End------" << endl;

  PHTFileServer::get().cd( _outfile );
  _events->Write();
  //PHTFileServer::get().close();

	//if ( trutheval ) trutheval->Delete();
	//if ( trackeval ) trackeval->Delete();

  return 0;
}


//----------------------------------------------------------------------------//
//-- fill_tree():
//--   Fill the trees with truth, track fit, and cluster information
//----------------------------------------------------------------------------//
void AnaEventGen::fill_tree(PHCompositeNode *topNode)
{

  // Make sure to reset all the TTree variables before trying to set them.
  reset_variables();

	//! read event variable

	int njet04 = 0;

	//! check tracks in jet
	if ( _jetmap04 ){
		for(JetMap::ConstIter iter = _jetmap04->begin(); iter != _jetmap04->end(); ++iter)
		{
			Jet *jet_true = iter->second;
			if ( jet_true->get_pt()<10.0 ) continue; 
			if ( jet_true->get_pt()<_jet_pt_min ) continue; 
			if ( fabs(jet_true->get_eta())>2.5 ) continue; 

			//cout << "JET04 ID: " << jet_true->get_id() << ", # of comp: " << jet_true->size_comp() << ", pT: " << jet_true->get_pt() << endl;
			jet04_pT[njet04] = jet_true->get_pt();
			jet04_eta[njet04] = jet_true->get_eta();
			jet04_phi[njet04] = jet_true->get_phi();

			njet04++;

			if ( njet04>=100 ) break;
		}//JepMap
	}//_jetmap

	int njet08 = 0;

	//! check tracks in jet
	if ( _jetmap08 ){
		for(JetMap::ConstIter iter = _jetmap08->begin(); iter != _jetmap08->end(); ++iter)
		{
			Jet *jet_true = iter->second;
			if ( jet_true->get_pt()<10.0 ) continue; 
			if ( jet_true->get_pt()<_jet_pt_min ) continue; 
			if ( fabs(jet_true->get_eta())>2.5 ) continue; 

			//cout << "JET08 ID: " << jet_true->get_id() << ", # of comp: " << jet_true->size_comp() << ", pT: " << jet_true->get_pt() << endl;
			jet08_pT[njet08] = jet_true->get_pt();
			jet08_eta[njet08] = jet_true->get_eta();
			jet08_phi[njet08] = jet_true->get_phi();

			njet08++;

			if ( njet08>=100 ) break;
		}//JepMap
	}//_jetmap

	if ( _jet_pt_min>0 && njet04<1 ) return;

	int npart = 0;

	if ( _hepmc_event_map ){
	//
		for (PHHepMCGenEventMap::ConstIter iter=_hepmc_event_map->begin(); iter!=_hepmc_event_map->end(); ++iter)
		{
			PHHepMCGenEvent *phhepmc_event = iter->second;
			HepMC::GenEvent *hepmc_event = phhepmc_event->getEvent();

			HepMC::HeavyIon *hi = hepmc_event->heavy_ion();

			ncoll = hi->Ncoll();
			ip = hi->impact_parameter();

			/*
			cout 
				<< hi->Ncoll() << " "
				<< hi->Npart_proj() << " "
				<< hi->Npart_targ() << " "
				<< hi->impact_parameter() << " "
				<< endl;
			*/

			for ( HepMC::GenEvent::particle_iterator p=hepmc_event->particles_begin(); p!=hepmc_event->particles_end(); ++p ){   

				HepMC::GenParticle *part = hepmc_event->barcode_to_particle((*p)->barcode());

				int tmp_pid = part->pdg_id();
				int tmp_status = part->status();

				if ( tmp_status!=1 ) continue;
				if ( !(abs(tmp_pid)==211 || abs(tmp_pid)==321 || abs(tmp_pid)==2212) ) continue;

				TVector3 vec(part->momentum().x(), part->momentum().y(), part->momentum().z());
				if ( fabs(vec.Eta())>6.0 ) continue;
				if ( vec.Pt()<_pt_min && fabs(vec.Eta())<3.0 ) continue;

				part_px[npart] = part->momentum().x();
				part_py[npart] = part->momentum().y();
				part_pz[npart] = part->momentum().z();

				part_pid[npart] = tmp_pid;

				npart++;

				if ( npart>=2000 ) break;

			}//particle 
		}//

	}//hepmc_event_map

	if ( 0 ){

		PHG4TruthInfoContainer::ConstRange range =  _truth_container->GetPrimaryParticleRange();

		for (PHG4TruthInfoContainer::ConstIterator iter=range.first; iter!=range.second; ++iter)
		{

			PHG4Particle *g4particle = iter->second;
			//PHG4VtxPoint *vtx = trutheval->get_vertex(g4particle);

			int pid = g4particle->get_pid();
			if ( !(abs(pid)==211 || abs(pid)==321 || abs(pid)==2212) ) continue;

			TVector3 vec(g4particle->get_px(), g4particle->get_py(), g4particle->get_pz());

			if ( fabs(vec.Eta())>5.0 ) continue;


			part_pid[npart]  = pid;
			//part_prod_x[npart] = vtx->get_x();
			//part_prod_y[npart] = vtx->get_y();
			//part_prod_z[npart] = vtx->get_z();

			part_px[npart] = g4particle->get_px();
			part_py[npart] = g4particle->get_py();
			part_pz[npart] = g4particle->get_pz();
			part_status[npart] = g4particle->get_primary_id();

			npart++;
		}
	}//truth_container

	
	int count = 0;
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(&_event);

	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(&ncoll);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(&ip);

	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(&npart);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(part_pid);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(part_px);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(part_py);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(part_pz);

	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(&njet04);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(jet04_pT);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(jet04_phi);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(jet04_eta);

	/*
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(&njet08);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(jet08_pT);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(jet08_phi);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(jet08_eta);
	*/

	_events->Fill();

  return;

}

//----------------------------------------------------------------------------//
//-- reset_variables():
//--   Reset all the tree variables to their default values.
//--   Needs to be called at the start of every event
//----------------------------------------------------------------------------//
void AnaEventGen::reset_variables()
{

	ncoll = 0;
	ip = -1;

	for (int ii=0; ii<2000; ii++){
		part_pid[ii] = part_status[ii] = -999;
		part_px[ii] = part_py[ii] = part_pz[ii] = -999;
		part_prod_x[ii] = part_prod_y[ii] = part_prod_z[ii] = -999;
		part_end_x[ii] = part_end_y[ii] = part_end_z[ii] = -999;
		part_prod_id[ii] = part_end_id[ii] = -999;
	}//ii

	for (int ii=0; ii<100; ii++){
		jet04_pT[ii] = jet04_eta[ii] = jet04_phi[ii] = -999;
		jet08_pT[ii] = jet08_eta[ii] = jet08_phi[ii] = -999;
	}//ii


}

//----------------------------------------------------------------------------//
//-- GetNodes():
//--   Get all the all the required nodes off the node tree
//----------------------------------------------------------------------------//
void AnaEventGen::GetNodes(PHCompositeNode * topNode)
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
	_hepmc_event_map = findNode::getClass<PHHepMCGenEventMap>(topNode,"PHHepMCGenEventMap");
	if (!_hepmc_event_map && _event<2)
	{
		cout << PHWHERE << " PHHepMCGenEventMap node not found on node tree" << endl;
	}

	// JetMap
	_jetmap04 = findNode::getClass<JetMap>(topNode,"AntiKt_Truth_r04");
	if (!_jetmap04 && _event < 2)
	{
		cout << PHWHERE << " AntiKt_Truth_r04 node not found on node tree" << endl;
	}

	// JetMap
	_jetmap08 = findNode::getClass<JetMap>(topNode,"AntiKt_Truth_r08");
	if (!_jetmap08 && _event < 2)
	{
		cout << PHWHERE << " AntiKt_Truth_r08 node not found on node tree" << endl;
	}

}



