#pragma once
#include "Cplex.h"
#include "_g.h"


class DiCutF :
    public Cplex
{
private:
	bool second_obj = false;

public:
	NumVarMatrix x;

	DiCutF(_g*, bool);
	~DiCutF();


	void DefVar();
	void DefVarX();

	void AddObj();
	void AddObj2();

	void AddCons();
	void AddConsOrderTrees();
	void AddConsPicknkArcs(); // pick n-k edges	
	void AddConsEachToTree(); // assign each edge to a tree	
	void AddConsArcsOfVertex();	

	DiCutF* Run();
	DiCutF* PrintSol();
	DiCutF* PrintModel();
	DiCutF* SetPrintCuts(bool);
	DiCutF* SetPrintCycles(bool);


	DiCutF* SetInitSol();
	DiCutF* ForceSol();

	DiCutF* SetInteger();
	DiCutF* SetLinear();

	DiCutF* SetSecondObj();

	void SaveOpt();
};

