#pragma once
#include "Cplex.h"
#include "_g.h"

class CutF : public Cplex
{

public:
	NumVarMatrix x;	

	int usercallback_count = 0;
	int lazycallback_count = 0;

	// callback arrays
	double** opt_x ;
	int* S;

	CutF(_g*, bool, bool);
	~CutF();

	void DefVar();
	void DefVarX();	

	void AddObj();

	void AddCons();
	void AddConsOrderTrees();
	void AddConstraintPicknkV(); // pick n-k edges
	void AddConsForEachSubset();
	void AddConsEachToTree(); // assign each edge to a tree
	
	void AddConsTwoEdgesOfVertex();
	void AddConsEdgesOfVertex();
	void AddConsEdgesOfVertex_Improved();
	
	void AddConsAtLeastOneEdge(); // at least one edge
	void AddConsBoundByUB();

	CutF* Run();
	CutF* PrintSol();
	CutF* PrintModel();
	CutF* SetPrintCuts(bool);	
	CutF* SetPrintCycles(bool);
	

	CutF* SetInitSol();
	CutF* ForceSol();
	CutF* InjectSol();
	
	CutF* SetInteger();
	CutF* SetLinear();

	void SaveOpt();

};

