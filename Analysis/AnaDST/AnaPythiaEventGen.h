#ifndef __AnaPythiaEventGen_H__
#define __AnaPythiaEventGen_H__

#include <fun4all/SubsysReco.h>
#include <string>
#include <fstream>

//Forward declerations
class PHCompositeNode;
class PHG4TruthInfoContainer;
class SvtxClusterMap;
class SvtxEvalStack;
class PHHepMCGenEventMap;
class SvtxVertexMap;
class SvtxTrackMap;
class JetMap;
class TFile;
class TTree;


//Brief: basic ntuple and histogram creation for sim evaluation
class AnaPythiaEventGen: public SubsysReco
{
 public: 
  //Default constructor
  AnaPythiaEventGen(const std::string &name="AnaPythiaEventGen");

  //Initialization, called for initialization
  int Init(PHCompositeNode *);

  //Process Event, called for each event
  int process_event(PHCompositeNode *);

  //End, write and close files
  int End(PHCompositeNode *);

  //Change output filename
  void set_filename(const char* file)
  { if(file) _outfile = file; }

	void set_pt_min(float min)
	{ _pt_min = min; }

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

	float _pt_min;

  //TTrees
  TTree* _events;
  int event;

	//primary particles
	int part_pid[2000], part_status[2000];
	float part_px[2000], part_py[2000], part_pz[2000];
	float part_eta[2000], part_phi[2000], part_pt[2000];
	float part_prod_x[2000], part_prod_y[2000], part_prod_z[2000];
	float part_end_x[2000], part_end_y[2000], part_end_z[2000];
	int part_prod_id[2000], part_end_id[2000];

	//truth jet
	float jet04_px[100], jet04_py[100], jet04_pz[100];
	float jet04_eta[100], jet04_phi[100], jet04_pt[100];
	float jet07_px[100], jet07_py[100], jet07_pz[100];

	//event info
	short ncoll;
	float ip;

  //Node pointers
	PHHepMCGenEventMap *_hepmc_event_map;
	JetMap *_jetmap_04;
	JetMap *_jetmap_07;

};

#endif //* __AnaPythiaEventGen_H__ *//
