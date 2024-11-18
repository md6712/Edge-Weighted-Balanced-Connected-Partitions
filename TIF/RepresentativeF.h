#pragma once
#include "Cplex.h"
#include "_g.h"

class RepresentativeF :    public Cplex
{
public:
	NumVarMatrix x;	
	IloNumVar z;


	RepresentativeF(_g*, bool);
	~RepresentativeF();

	void DefVar();
	void DefVarX();	

	void AddObj();

	void AddCons();
	

	RepresentativeF* Run();
	RepresentativeF* PrintSol();
	RepresentativeF* PrintModel();
	void SaveOpt();

private :

	void AddConsComputeZ ();
	void AddConsKeyEdges ();
	void AddConsIFeeTHENef ();
	void AddConsBelongToOneTree ();
	void AddConsTotalNumberOfEdgesSelected ();
	void AddConsAllCycles ();
};

