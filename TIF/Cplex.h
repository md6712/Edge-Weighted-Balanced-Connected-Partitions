#pragma once
#include <stdio.h>
#include <conio.h>
#include <ilcplex/ilocplex.h>
#include "_g.h"

#define TIME_LIMIT 10800 // 3 hours in seconds

using namespace std;

ILOSTLBEGIN

typedef IloArray<IloNumVarArray> NumVarMatrix;
typedef IloArray<NumVarMatrix> NumVarCube;
typedef IloArray<NumVarCube> NumVarCubeArr;
typedef IloArray<IloRangeArray> RangeMatrix;

class Cplex
{
public:

	double gap = 0.01;		// optimality gap
	double opt = 0;			// optimal value
	double lp_bound = 0;	// LP bound
	double root_bound_start = 0;	// root bound start
	double root_bound_end = 0;		//	root bound end

	bool integer = true;		
	bool printCycles;
	bool printCuts;

	bool force_silent = false; // force silent mode

	double time_limit = TIME_LIMIT; // time limit in seconds
	
	char name[50];
	char cut[500];
	char cycle[500];
	_g* instance;
	IloEnv env;
	IloModel model;	
	IloCplex cplex;

	// set Integer vs Linear
	void SetLinear();
	void SetInteger();


	Cplex(_g*);
	~Cplex();
	void CplexSettings();

	Cplex* Run();
	void End();
	virtual void DefVar();
	virtual void AddObj();
	virtual void AddCons();

	void SetPrintCuts(bool);
	void SetPrintCycles(bool);


	// turn off warnings and output
	void* SetQuiet();

	// turn off all cplex cuts
	void* SetNoCuts();

	// turn  cplex cuts to default
	void* SetCutsDefault();

	// set time limit	
	void SetTimeLimit(double time_limit) {
		this->time_limit = time_limit;
		cplex.setParam(IloCplex::TiLim, time_limit);
	}

	// reset time limit
	void ResetTimeLimit() {
		this->time_limit = TIME_LIMIT;
		cplex.setParam(IloCplex::TiLim, TIME_LIMIT);
	}


};
