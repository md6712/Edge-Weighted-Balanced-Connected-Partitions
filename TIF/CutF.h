#pragma once
#include "Cplex.h"
#include "_g.h"

class CutF : public Cplex
{

public:
	NumVarMatrix x;	

	CutF(_g*, bool);
	~CutF();


	void DefVar();
	void DefVarX();	

	void AddObj();
	void AddObj2();

	void AddCons();
	void AddConsOrderTrees();
	void AddConstraintPicknkV(); // pick n-k edges
	void AddConsForEachSubset();
	void AddConsEachToTree(); // assign each edge to a tree
	void AddConsTwoEdgesOfVertex();
	void AddConsEdgesOfVertex();
	void AddConsAtLeastOneEdge(); // at least one edge

	CutF* Run();
	CutF* PrintSol();
	CutF* PrintModel();
	CutF* SetPrintCuts(bool);	
	CutF* SetPrintCycles(bool);
	

	CutF* SetInitSol();
	CutF* ForceSol();

	void SaveOpt();

};

