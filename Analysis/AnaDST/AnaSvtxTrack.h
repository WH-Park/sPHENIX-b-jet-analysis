#ifndef __AnaSvtxTrack_H__
#define __AnaSvtxTrack_H__

#include <fun4all/SubsysReco.h>
#include <string>

//Forward declerations
class PHCompositeNode;
class PHG4TruthInfoContainer;
class SvtxClusterMap;
class SvtxEvalStack;
class PHHepMCGenEvent;
class SvtxVertexMap;
class SvtxTrackMap;
class JetMap;
class TFile;
class TTree;


//Brief: basic ntuple and histogram creation for sim evaluation
class AnaSvtxTrack: public SubsysReco
{
 public: 
  //Default constructor
  AnaSvtxTrack(const std::string &name="AnaSvtxTrack");

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

	short gen_pid[2000];
	short gen_status[2000];

	float gen_px[2000];
	float gen_py[2000];
	float gen_pz[2000];

	float gen_x[2000];
	float gen_y[2000];
	float gen_z[2000];

	short gen_embed[2000];
	short gen_primary[2000];
	short gen_ngmaps[2000];
	short gen_ngintt[2000];
	short gen_ngtpc[2000];

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
	short track_cluster_n_maps_truth[2000];
	short track_cluster_n_intt_truth[2000];
	short track_cluster_n_tpc_truth[2000];
	short track_cluster_pattern_maps[2000];
	short track_cluster_pattern_intt[2000];
	//short track_cluster_layer[2000][100];
	//short track_cluster_gembed[2000][100];
	short track_cluster_gflavor[2000];

	int nvertex;
	short vtx_id[10];
	float vtx_gen[10][3];
	float vtx_reco[10][3];
	float vtx_reco_err[10][3];

	int _n_maps_layers;
	int _n_intt_layers;

  //Node pointers
  PHG4TruthInfoContainer *_truth_container;
	PHHepMCGenEvent *_hepmc_event;
	SvtxVertexMap *_vtxmap;
	SvtxTrackMap *_trkmap;
	SvtxClusterMap *_clustermap;
	SvtxEvalStack *_svtxevalstack;

};

#endif //* __AnaSvtxTrack_H__ *//
