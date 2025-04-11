#include "Cplex.h"

Cplex::Cplex(_g* instance) {
	this->instance = instance;
	env = IloEnv();
	model = IloModel(env);
	cplex = IloCplex(model);
}

Cplex::~Cplex() {
	cplex.end();
	model.end();
	env.end();
}

void Cplex::CplexSettings() {
	//cplex.setOut(env.getNullStream());

	// turn off cplex warnings
	cplex.setWarning(env.getNullStream());
	
	//cplex.setParam(cplex.NodeFileInd, 3);
	cplex.setParam(cplex.TiLim, 1800);
	//cplex.setParam(cplex.HeurFreq, -1);
	cplex.setParam(cplex.Threads, 1);
	cplex.setParam(IloCplex::EpGap, 0.00001);
	////cplex.setParam(cplex.MIPEmphasis, 4);

	////cplex.setParam(cplex.Cliques, -1);
	////cplex.setParam(cplex.Covers, -1);
	////cplex.setParam(cplex.DisjCuts, -1);
	////cplex.setParam(cplex.FlowCovers, -1);
	////cplex.setParam(cplex.FlowPaths, -1);
	////cplex.setParam(cplex.FracCuts, -1);
	////cplex.setParam(cplex.GUBCovers, -1);
	////cplex.setParam(cplex.ImplBd, -1);
	////cplex.setParam(cplex.MIRCuts, -1);
	////cplex.setParam(cplex.MCFCuts, -1);
	////cplex.setParam(cplex.ZeroHalfCuts, -1);

	////cplex.setParam(cplex.PreInd, 0);
	////cplex.setParam(cplex.AggInd,0);
	////cplex.setParam(cplex.BndStrenInd, 0);
	////cplex.setParam(cplex.CoeRedInd, 0);
	////cplex.setParam(cplex.RelaxPreInd, 0);
	////cplex.setParam(cplex.Reduce, 0);
	////cplex.setParam(cplex.PrePass, 0);
	////cplex.setParam(cplex.RepeatPresolve, 0);
	////cplex.setParam(cplex.HeurFreq, -1);
	////cplex.setParam(cplex.MIPSearch, 1);
}

Cplex* Cplex::Run() {
	CplexSettings();
	if (cplex.solve()) {
		if (cplex.getStatus() == IloAlgorithm::Optimal)
		{
			gap = 0;
			opt = (double)cplex.getObjValue();
		}
		else {
			gap = (double)cplex.getMIPRelativeGap();			
			opt = (double)cplex.getBestObjValue();			
		}
	}
	else {
		if (cplex.getStatus() == IloAlgorithm::Infeasible) {
			gap = -1;
			opt = -1;
			printf("\nInfeasible");
		}
	}
	return this;
}



void Cplex::End() {
	env.end();
}

void Cplex::DefVar() {

}

void Cplex::AddObj() {

}

void Cplex::AddCons() {

}

void Cplex::SetInteger() {
	integer = true;
}	

void Cplex::SetLinear() {
	integer = false;
}

// print cuts
void Cplex::SetPrintCuts(bool printCuts) {
	this->printCuts = printCuts;
}

// print cycles
void Cplex::SetPrintCycles(bool printCycles) {
	this->printCycles = printCycles;
}

void* Cplex::SetQuiet() {
	cplex.setOut(env.getNullStream());
	cplex.setWarning(env.getNullStream());
	return this;
}