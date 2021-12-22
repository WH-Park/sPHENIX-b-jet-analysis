#ifndef __AnaSvtxVertex_H__
#define __AnaSvtxVertex_H__

#include <fun4all/SubsysReco.h>
#include <string>

namespace genfit
{
	class GFRaveVErtex;
	class GFRaveVertexFactory;
	class Track;
}//namespace PHGenFit

namespace PHGenFit
{
	class Fitter;
}//namespace PHGenFit

//Forward declerations
class PHCompositeNode;
class PHG4TruthInfoContainer;
class SvtxClusterMap;
class SvtxEvalStack;
class PHHepMCGenEventMap;
class SvtxVertexMap;
class SvtxTrack;
class SvtxTrackMap;
class JetMap;
class TFile;
class TH2D;
class TTree;


//Brief: basic ntuple and histogram creation for sim evaluation
class AnaSvtxVertex: public SubsysReco
{
 public: 
  //Default constructor
  AnaSvtxVertex(const std::string &name="AnaSvtxVertex");

  //Initialization, called for initialization
  int Init(PHCompositeNode *);

	int InitRun(PHCompositeNode *);

  //Process Event, called for each event
  int process_event(PHCompositeNode *);

  //End, write and close files
  int End(PHCompositeNode *);

  //Change output filename
  void set_filename(const char* file)
  { if(file) _outfile = file; }

  //User modules
  void fill_tree(PHCompositeNode*);
  void reset_variables();

 private:
  //output filename
  std::string _outfile;
   
  //Event counter
  int _event;

  //Get all the nodes
  void GetNodes(PHCompositeNode *);
  
  //flags
  unsigned int _flags;

  //TTrees
  TTree* _events;
  int event;

	//generated vertex
	float vtx_gen[3];

	//reconstructed vertex
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

	float tmp_prod_x, tmp_prod_y, tmp_prod_z;
	int tmp_prod_pid;

	//generated jet
	int jet04_prop_parton[100];
	int jet04_prop_hadron[100];
	float jet04_pt[100];
	float jet04_eta[100];
	float jet04_mass[100];

	float jet08_pt[100];
	float jet08_eta[100];
	float jet08_mass[100];

	int rave_sv_pT10_nvtx[100];
	float rave_sv_pT10_vtx_x[100][30];
	float rave_sv_pT10_vtx_y[100][30];
	float rave_sv_pT10_vtx_z[100][30];
	float rave_sv_pT10_vtx_pt[100][30];
	float rave_sv_pT10_vtx_mass[100][30];
	float rave_sv_pT10_vtx_mass_corr[100][30];
	int rave_sv_pT10_vtx_ntrack[100][30];


  //Node pointers
  PHG4TruthInfoContainer *_truth_container;
	PHHepMCGenEventMap *_hepmc_event_map;
	SvtxVertexMap *_vtxmap;
	SvtxVertexMap *_vtxmap_rave;
	SvtxVertexMap *_vtxmap_rave_refit;
	SvtxVertexMap *_vtxmap_acts;
	SvtxVertexMap *_vtxmap_acts_refit;

	SvtxEvalStack *_svtxevalstack;

	SvtxTrackMap *_trkmap;

	JetMap *_jetmap04;
	JetMap *_jetmap08;

	TH2D *h2d_eta_pt_p_mvtx0;
	TH2D *h2d_eta_pt_p_mvtx1;
	TH2D *h2d_eta_pt_p_mvtx2;
	TH2D *h2d_eta_pt_p_mvtx3;

	TH2D *h2d_eta_pt_s_mvtx0;
	TH2D *h2d_eta_pt_s_mvtx1;
	TH2D *h2d_eta_pt_s_mvtx2;
	TH2D *h2d_eta_pt_s_mvtx3;

	genfit::Track *TranslateSvtxToGenFitTrack(SvtxTrack *svtx);
	PHGenFit::Fitter *_fitter;
	genfit::GFRaveVertexFactory *_vertex_finder;

};

#endif //* __AnaSvtxVertex_H__ *//
