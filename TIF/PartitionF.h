
#pragma once
#include "Cplex.h"
#include "_g.h"

class PartitionF : public Cplex
{
	public:
	NumVarMatrix x;
	NumVarMatrix y;
	IloNumVar z;
	
	bool printCuts;

	PartitionF(_g*, bool);

	void DefVar();
	void DefVarX();
	void DefVarY();
	void DefVarZ();

	void AddObj();

	void AddCons();
	void AddConsEachVertexInATree();
	void AddConsOrderTrees();
	void AddConsXYRelation();
	void AddConsMinSeperators();

	PartitionF* Run();
	PartitionF* PrintSol();
	PartitionF* PrintModel();
	void SaveOpt();
};



