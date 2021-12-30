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


#include "AnaPythiaEventGen.h"

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
AnaPythiaEventGen::AnaPythiaEventGen(const string &name):
  SubsysReco( name ),
  _events( NULL )
{
  //initialize
  _event = 0;
  _outfile = "AnaPythiaEventGen.root";

	_hepmc_event_map = NULL;

}

//----------------------------------------------------------------------------//
//-- Init():
//--   Intialize all histograms, trees, and ntuples
//----------------------------------------------------------------------------//
int AnaPythiaEventGen::Init(PHCompositeNode *topNode)
{
  cout << PHWHERE << " Openning file " << _outfile << endl;
  PHTFileServer::get().open( _outfile, "RECREATE");

  // create TTree
	_events = new TTree("T", "T");
	_events->Branch("event",0,"event/I");

	_events->Branch("ncoll",0,"ncoll/S");
	_events->Branch("ip",0,"ip/F");

	_events->Branch("np",0,"np/I");
	_events->Branch("part_pid",0,"part_pid[np]/I");
	//_events->Branch("part_px",0,"part_px[np]/F");
	//_events->Branch("part_py",0,"part_py[np]/F");
	//_events->Branch("part_pz",0,"part_pz[np]/F");
	_events->Branch("part_eta",0,"part_eta[np]/F");
	_events->Branch("part_phi",0,"part_phi[np]/F");
	_events->Branch("part_pt",0,"part_pt[np]/F");

	_events->Branch("njet04",0,"njet04/I");
	//_events->Branch("jet04_px",0,"jet04_px[njet04]/F");
	//_events->Branch("jet04_py",0,"jet04_py[njet04]/F");
	//_events->Branch("jet04_pz",0,"jet04_pz[njet04]/F");
	_events->Branch("jet04_eta",0,"jet04_eta[njet04]/F");
	_events->Branch("jet04_phi",0,"jet04_phi[njet04]/F");
	_events->Branch("jet04_pt",0,"jet04_pt[njet04]/F");

	/*
	_events->Branch("njet07",0,"njet07/I");
	_events->Branch("jet07_px",0,"jet07_px[njet07]/F");
	_events->Branch("jet07_py",0,"jet07_py[njet07]/F");
	_events->Branch("jet07_pz",0,"jet07_pz[njet07]/F");
	*/

  return 0;
}

//----------------------------------------------------------------------------//
//-- process_event():
//--   Call user instructions for every event.
//--   This function contains the analysis structure.
//----------------------------------------------------------------------------//
int AnaPythiaEventGen::process_event(PHCompositeNode *topNode)
{
  _event++;
  if (_event % 100 == 0)
    cout << PHWHERE << "Events processed: " << _event << endl;

  GetNodes(topNode);

  fill_tree(topNode);

  return 0;
}

//----------------------------------------------------------------------------//
//-- End():
//--   End method, wrap everything up
//----------------------------------------------------------------------------//
int AnaPythiaEventGen::End(PHCompositeNode *topNode)
{

	cout << "-----AnaPythiaEventGen::End------" << endl;

  PHTFileServer::get().cd( _outfile );
  //_events->Write();
	PHTFileServer::get().write( _outfile );
	PHTFileServer::get().close();

	//if ( trutheval ) trutheval->Delete();
	//if ( trackeval ) trackeval->Delete();

  return 0;
}


//----------------------------------------------------------------------------//
//-- fill_tree():
//--   Fill the trees with truth, track fit, and cluster information
//----------------------------------------------------------------------------//
void AnaPythiaEventGen::fill_tree(PHCompositeNode *topNode)
{

  // Make sure to reset all the TTree variables before trying to set them.
  reset_variables();

	//! read event variable

	int npart = 0;

	if ( _hepmc_event_map ){
	//
		for (PHHepMCGenEventMap::ConstIter iter=_hepmc_event_map->begin(); iter!=_hepmc_event_map->end(); ++iter)
		{
			PHHepMCGenEvent *phhepmc_event = iter->second;
			HepMC::GenEvent *hepmc_event = phhepmc_event->getEvent();

			HepMC::HeavyIon *hi = hepmc_event->heavy_ion();

			if ( hi!=NULL ){
				ncoll = hi->Ncoll();
				ip = hi->impact_parameter();
			}

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

				part_eta[npart] = vec.Eta();
				part_phi[npart] = vec.Phi();
				part_pt[npart] = vec.Pt();

				part_pid[npart] = tmp_pid;

				npart++;

				if ( npart>=2000 ) break;

			}//particle 
		}//

	}//hepmc_event_map

	int njets_04 = 0;
	for (JetMap::ConstIter iter=_jetmap_04->begin(); iter!=_jetmap_04->end(); ++iter)
	{
		Jet *jet = iter->second;
		if ( jet->get_pt()<15.0 ) continue;
		//if ( fabs(jet->get_eta())>_cut_jet_eta ) continue;

		jet04_px[njets_04] = jet->get_px();
		jet04_py[njets_04] = jet->get_py();
		jet04_pz[njets_04] = jet->get_pz();

		jet04_eta[njets_04] = jet->get_eta();
		jet04_phi[njets_04] = jet->get_phi();
		jet04_pt[njets_04] = jet->get_pt();

		njets_04++;

		if ( njets_04>=100 ) break;
	}

	int njets_07 = 0;
	for (JetMap::ConstIter iter=_jetmap_07->begin(); iter!=_jetmap_07->end(); ++iter)
	{
		Jet *jet = iter->second;
		if ( jet->get_pt()<15.0 ) continue;
		//if ( fabs(jet->get_eta())>_cut_jet_eta ) continue;

		jet07_px[njets_07] = jet->get_px();
		jet07_py[njets_07] = jet->get_py();
		jet07_pz[njets_07] = jet->get_pz();

		njets_07++;

		if ( njets_07>=100 ) break;
	}

	
	int count = 0;
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(&_event);

	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(&ncoll);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(&ip);

	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(&npart);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(part_pid);
	//((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(part_px);
	//((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(part_py);
	//((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(part_pz);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(part_eta);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(part_phi);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(part_pt);

	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(&njets_04);
	//((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(jet04_px);
	//((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(jet04_py);
	//((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(jet04_pz);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(jet04_eta);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(jet04_phi);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(jet04_pt);

	/*
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(&njets_07);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(jet07_px);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(jet07_py);
	((TBranch*) _events->GetListOfBranches()->At(count++))->SetAddress(jet07_pz);
	*/

	_events->Fill();

  return;

}

//----------------------------------------------------------------------------//
//-- reset_variables():
//--   Reset all the tree variables to their default values.
//--   Needs to be called at the start of every event
//----------------------------------------------------------------------------//
void AnaPythiaEventGen::reset_variables()
{

	ncoll = 0;
	ip = -1;

	for (int ii=0; ii<2000; ii++){
		part_pid[ii] = part_status[ii] = -999;
		part_eta[ii] = part_phi[ii] = part_pt[ii] = -999;
		part_px[ii] = part_py[ii] = part_pz[ii] = -999;
		part_prod_x[ii] = part_prod_y[ii] = part_prod_z[ii] = -999;
		part_end_x[ii] = part_end_y[ii] = part_end_z[ii] = -999;
		part_prod_id[ii] = part_end_id[ii] = -999;
	}//ii

	for (int ii=0; ii<100; ii++){
		jet04_px[ii] = jet04_py[ii] = jet04_pz[ii] = -999;
		jet04_eta[ii] = jet04_phi[ii] = jet04_pt[ii] = -999;
		jet07_px[ii] = jet07_py[ii] = jet07_pz[ii] = -999;
	}//ii


}

//----------------------------------------------------------------------------//
//-- GetNodes():
//--   Get all the all the required nodes off the node tree
//----------------------------------------------------------------------------//
void AnaPythiaEventGen::GetNodes(PHCompositeNode * topNode)
{
  //DST objects
	//HepMCGenEventMap
	_hepmc_event_map = findNode::getClass<PHHepMCGenEventMap>(topNode,"PHHepMCGenEventMap");
	if (!_hepmc_event_map && _event<2)
	{
		cout << PHWHERE << " PHHepMCGenEventMap node not found on node tree" << endl;
	}

	_jetmap_04 = findNode::getClass<JetMap>(topNode,"AntiKt_Truth_r04");
	if (!_jetmap_04 && _event<2){
		cout << PHWHERE << " JetMap node not found on node tree" << endl;
	}

	_jetmap_07 = findNode::getClass<JetMap>(topNode,"AntiKt_Truth_r07");
	if (!_jetmap_07 && _event<2){
		cout << PHWHERE << " JetMap node not found on node tree" << endl;
	}

}



