#ifndef __AnaTrigSim_H__
#define __AnaTrigSim_H__

#include <fun4all/SubsysReco.h>
#include <string>

//Forward declerations
class PHCompositeNode;
class CaloTriggerInfo;
class RawClusterContainer;
class PHG4TruthInfoContainer;
class CaloEvalStack;
class TFile;
class TTree;


//Brief: basic ntuple and histogram creation for sim evaluation
class AnaTrigSim: public SubsysReco
{
 public: 
  //Default constructor
  AnaTrigSim(const std::string &name="AnaTrigSim");

  //Initialization, called for initialization
  int Init(PHCompositeNode *);

  //Process Event, called for each event
  int process_event(PHCompositeNode *);

  //End, write and close files
  int End(PHCompositeNode *);

  //Change output filename
  void set_filename(const char* file)
  { if(file) _outfile = file; }

	void set_truncation( int emulate_truncation );

  //User modules
  void fill_tree(PHCompositeNode*);
  void reset_variables();

 private:
  //output filename
  std::string _outfile;

	int _emulate_truncation;
   
  //Event counter
  int _event;

  //Get all the nodes
  void GetNodes(PHCompositeNode *);
  
  //flags
  unsigned int _flags;

  //TTrees
  TTree* _T;

  //Node pointers
	CaloTriggerInfo *_triggerinfo;
	PHG4TruthInfoContainer *_truthcontainer;
	RawClusterContainer *_cemc_clusters;
	CaloEvalStack *_caloevalstack;

};

#endif //* __AnaTrigSim_H__ *//
