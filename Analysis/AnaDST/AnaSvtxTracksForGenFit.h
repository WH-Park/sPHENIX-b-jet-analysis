#ifndef __AnaSvtxTracksForGenFit_H__
#define __AnaSvtxTracksForGenFit_H__

#include <fun4all/SubsysReco.h>
#include <string>

//Forward declerations
class PHCompositeNode;
class PHG4TruthInfoContainer;
class SvtxClusterMap;
class SvtxEvalStack;
class TFile;
class TTree;


//Brief: basic ntuple and histogram creation for sim evaluation
class AnaSvtxTracksForGenFit: public SubsysReco
{
 public: 
  //Default constructor
  AnaSvtxTracksForGenFit(const std::string &name="AnaSvtxTracksForGenFit");

  //Initialization, called for initialization
  int Init(PHCompositeNode *);

  //Process Event, called for each event
  int process_event(PHCompositeNode *);

  //End, write and close files
  int End(PHCompositeNode *);

  //Change output filename
  void set_filename(const char* file)
  { if(file) _outfile = file; }

  //Flags of different kinds of outputs
  enum Flag
  {
    //all disabled
    NONE = 0,
  };

  //Set the flag
  //Flags should be set like set_flag(AnaSvtxTracksForGenFit::TRUTH, true) from macro
  void set_flag(const Flag& flag, const bool& value)
  {
   if(value) _flags |= flag;
   else _flags &= (~flag);
  }

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
  TTree* _tracks;
  int event;
	int ntracks;
  //-- truth
  int gtrackID[100];
  int gflavor[100];
  float gpx[100];
  float gpy[100];
  float gpz[100];
  float gvx[100];
  float gvy[100];
  float gvz[100];
  //-- reco
  int trackID[100];
  int charge[100];
  int nhits[100];
  float px[100];
  float py[100];
  float pz[100];
  float dca2d[100];
  //-- clusters
  int clusterID[100][7];
  int layer[100][7];
	int adc[100][7];
  float x[100][7];
  float y[100][7];
  float z[100][7];
  float size_dphi[100][7];
  float size_dz[100][7];

  //Node pointers
  PHG4TruthInfoContainer* _truth_container;
  SvtxClusterMap* _clustermap;

  // eval stack
  SvtxEvalStack* _svtxevalstack;

};

#endif //* __AnaSvtxTracksForGenFit_H__ *//
