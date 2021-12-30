#ifndef __AnaPileupPP_H__
#define __AnaPileupPP_H__

#include <fun4all/SubsysReco.h>
#include <string>

//Forward declerations
class PHCompositeNode;
class PHG4TruthInfoContainer;
class SvtxClusterMap;
class SvtxEvalStack;
class PHHepMCGenEventMap;
class SvtxVertexMap;
class SvtxTrackMap;
class PHG4HitContainer;
class JetMap;
class TFile;
class TTree;
class TH1D;
class TH2D;
class TH3D;


//Brief: basic ntuple and histogram creation for sim evaluation
class AnaPileupPP: public SubsysReco
{
 public: 
  //Default constructor
  AnaPileupPP(const std::string &name="AnaPileupPP");

  //Initialization, called for initialization
  int Init(PHCompositeNode *);

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

	void set_n_maps_layers(int n){ _n_maps_layers = n; }
	void set_n_intt_layers(int n){ _n_intt_layers = n; }
	void set_n_tpc_layers(int n){ _n_tpc_layers = n; }

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

	//Histogram
	TH1D *hz_in;
	TH1D *hz_out;
	TH1D *hz_diff;
	TH1D *hn_out;
	TH2D *hpT_eta_in;
	TH2D *hpT_eta_in_tpc;
	TH2D *hpT_eta_in_tpc_mvtx;
	TH2D *hpT_eta_in_tpc_mvtx_intt1;
	TH2D *hpT_eta_in_tpc_mvtx_intt2;
	TH2D *hpT_eta_in_tpc_mvtx_intt3;
	TH2D *hpT_eta_in_tpc_mvtx_intt4;
	TH2D *hpT_eta_out;
	TH2D *hpT_eta_out_tpc;
	TH2D *hpT_eta_out_tpc_mvtx;
	TH2D *hpT_eta_out_tpc_mvtx_intt1;
	TH2D *hpT_eta_out_tpc_mvtx_intt2;
	TH2D *hpT_eta_out_tpc_mvtx_intt3;
	TH2D *hpT_eta_out_tpc_mvtx_intt4;

	TH1D *hz_out_trk_tpc;
	TH1D *hz_out_trk_tpc_mvtx;

	TH2D *hclus_intt_zphi_gen[8];
	TH2D *hclus_intt_zphi_reco[8];

	TH2D *hnclus_pT_tpc;
	TH2D *hnclus_pT_intt;
	TH2D *hnclus_pT_mvtx;

	TH2D *hnclus_pT_mvtx_tpc35;
	TH2D *hnclus_pT_intt_tpc35;

	//SvtxCluster
	/*
	short cluster_layer[5000];
	short cluster_gprimary[5000];
	short cluster_gembed[5000];
	short cluster_gflavor[5000];
	float cluster_gvx[5000];
	float cluster_gvy[5000];
	float cluster_gpx[5000];
	float cluster_gpy[5000];
	float cluster_err[5000][3];
	*/

	short gen_pid[20000];
	short gen_status[20000];

	float gen_px[20000];
	float gen_py[20000];
	float gen_pz[20000];

	float gen_vx[20000];
	float gen_vy[20000];
	float gen_vz[20000];
	float gen_vt[20000];

	short gen_embed[20000];
	short gen_primary[20000];
	short gen_ngmaps[20000];
	short gen_ngintt[20000];
	short gen_ngtpc[20000];
	short gen_ngmaps_layer[20000];
	short gen_ngintt_layer[20000];
	short gen_ngtpc_layer[20000];

	//G4Hit
	short g4hit_layer[2000];
	float g4hit_x[2000];
	float g4hit_y[2000];
	float g4hit_z[2000];

	//SvtxCluster
	short cluster_layer[2000];
	float cluster_x[2000];
	float cluster_y[2000];
	float cluster_z[2000];

	//SvtxTrack
	float track_pt[2000];
	float track_eta[2000];
	float track_chiq[2000];
	short track_ndf[2000];
	short track_gprimary[2000];
	short track_gembed[2000];
	short track_gflavor[2000];
	short track_matched_cluster_n[2000];
	short track_parent_id[2000];
	float track_gvx[2000];
	float track_gvy[2000];
	float track_gvz[2000];
	float track_gvt[2000];
	float track_gpx[2000];
	float track_gpy[2000];
	float track_gpz[2000];
	float track_dca2d[2000];
	float track_dca3d_xy[2000];
	float track_dca3d_z[2000];

	short track_cluster_n[2000];
	short track_cluster_n_maps[2000];
	short track_cluster_n_intt[2000];
	short track_cluster_n_tpc[2000];

	short track_cluster_n_layer[2000];
	short track_cluster_n_layer_maps[2000];
	short track_cluster_n_layer_intt[2000];
	short track_cluster_n_layer_tpc[2000];

	short track_cluster_n_maps_truth[2000];
	short track_cluster_n_intt_truth[2000];
	short track_cluster_n_tpc_truth[2000];

	short track_cluster_n_maps_embed[2000];
	short track_cluster_n_intt_embed[2000];
	short track_cluster_n_tpc_embed[2000];

	short track_cluster_pattern_maps[2000];
	short track_cluster_pattern_intt[2000];
	//short track_cluster_layer[2000][100];
	//short track_cluster_gembed[2000][100];
	short track_cluster_gflavor[2000];

	int nvertex;
	short vtx_id[10];
	float vtx_gen[3];
	float vtx_reco[10][3];
	float vtx_reco_err[10][3];

	int _n_maps_layers;
	int _n_intt_layers;
	int _n_tpc_layers;

  //Node pointers
  PHG4TruthInfoContainer *_truth_container;
	PHHepMCGenEventMap *_hepmc_eventmap;
	SvtxVertexMap *_vtxmap;
	SvtxTrackMap *_trkmap;
	SvtxClusterMap *_clustermap;
	SvtxEvalStack *_svtxevalstack;
	PHG4HitContainer *_g4hit_container;

};

#endif //* __AnaPileupPP_H__ *//
