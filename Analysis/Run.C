#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <G4_Magnet.C>
#include <trackreco/MakeActsGeometry.h>
#include "FROG.h"

//#include "anadst/AnaSvtxVertex.h"
#include "/gpfs/mnt/gpfs02/sphenix/user/whpark/Analysis/install/include/anadst/AnaSvtxVertex.h"
#include "/gpfs/mnt/gpfs02/sphenix/user/whpark/temp_truth_tagger/analysis/HF-Jet/TruthGeneration/HFJetTruthTrigger.h"

//#include "/gpfs/mnt/gpfs02/sphenix/user/whpark/test/install/include/trackreco/PHRaveVertexing.h"

//#include "/gpfs/mnt/gpfs02/sphenix/user/whpark/test/install/include/trackreco/PHActsInitialVertexFinder.h"
	
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libkfparticle_sphenix.so)
R__LOAD_LIBRARY(libHFJetTruthGeneration.so)
R__LOAD_LIBRARY(libAnaDST.so)

int Run(
		const int nEvents = 0,
		const char * inputFile = "G4sPHENIX.root"
		)
//G4sPHENIX_0301.root
{
  //---------------
  // Load libraries
  //---------------

	gSystem->Load("libfun4all.so");
  gSystem->Load("libg4dst.so");

  // Enabling file finding in dCache
  FROG*fr = new FROG();

  //---------------
  // Fun4All server
  //---------------

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);
  // just if we set some flags somewhere in this macro
  //recoConsts *rc = recoConsts::instance();
  // By default every random number generator uses
  // PHRandomSeed() which reads /dev/urandom to get its seed
  // if the RANDOMSEED flag is set its value is taken as seed
  // You ca neither set this to a random value using PHRandomSeed()
  // which will make all seeds identical (not sure what the point of
  // this would be:
  //  rc->set_IntFlag("RANDOMSEED",PHRandomSeed());
  // or set it to a fixed value so you can debug your code
  // rc->set_IntFlag("RANDOMSEED", 12345);

	MagnetInit();
	MagnetFieldInit();

	MakeActsGeometry* geom = new MakeActsGeometry();
	geom->Verbosity(0);
	geom->setMagField(G4MAGNET::magfield);
	geom->setMagFieldRescale(G4MAGNET::magfield_rescale);
	se->registerSubsystem(geom);

	/*
	PHActsInitialVertexFinder* mod3 = new PHActsInitialVertexFinder("MyVertexFinder");
	mod3->setSvtxTrackMapName("SvtxTrackMap");
	mod3->setSvtxVertexMapName("SvtxVertexMapActsRefit");
	mod3->setMaxVertices(1);
	mod3->Verbosity(0);
	se->registerSubsystem(mod3);

	PHActsInitialVertexFinder* mod4 = new PHActsInitialVertexFinder("MyVertexFinder");
	mod4->setSvtxTrackMapName("SvtxTrackMap");
	mod4->setSvtxVertexMapName("SvtxVertexMapActsRefitMvtx");
	mod4->setMaxVertices(1);
	mod4->requireMvtxCluster(true);
	mod4->Verbosity(0);
	se->registerSubsystem(mod4);

	PHRaveVertexing *mod2 = new PHRaveVertexing();
	mod2->set_svtxvertexmaprefit_node_name("SvtxVertexMapRaveRefitMvtx");
	se->registerSubsystem(mod2);
	*/

	HFJetTruthTrigger *jt_04 = new HFJetTruthTrigger("HFJetTruthTrigger.root",5,"AntiKt_Truth_r04");
	jt_04->set_eta_min(-0.6);
	jt_04->set_eta_max(+0.6);
	jt_04->set_pt_min(10.0);
	jt_04->set_pt_max(100.0);
	se->registerSubsystem(jt_04);

	AnaSvtxVertex *mod1 = new AnaSvtxVertex();
	se->registerSubsystem(mod1);

	Fun4AllInputManager *hitsin = new Fun4AllDstInputManager("DSTin");
	if (strstr(inputFile,"root")){
		hitsin->fileopen(inputFile);
	}else{
		hitsin->AddListFile(inputFile);
	}
	se->registerInputManager(hitsin);

	gSystem->ListLibraries();
  se->run(nEvents);

  //-----
  // Exit
  //-----

  se->End();
  std::cout << "All done" << std::endl;
  delete se;
  gSystem->Exit(0);
	return 0;
}
