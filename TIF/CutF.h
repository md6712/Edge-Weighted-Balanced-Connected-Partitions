#pragma once
#include "Cplex.h"
#include "_g.h"

class CutF : public Cplex
{

private: 	
	bool second_obj = false;


public:
	NumVarMatrix x;	

	int usercallback_count = 0;
	int lazycallback_count = 0;


	CutF(_g*, bool, bool);
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
	
	CutF* SetInteger();
	CutF* SetLinear();

	CutF* SetSecondObj();

	void SaveOpt();

};

