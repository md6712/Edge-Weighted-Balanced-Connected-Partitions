#include "Cplex.h"

Cplex::Cplex(_g* instance)
	: instance(instance),  // assign your graph
	env(),               // construct env
	model(env),          // construct model using env
	cplex(model)         // construct cplex using model
{

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
	
	// Use up to ~5 GB for working memory (leaves room for overhead)
	//cplex.setParam(IloCplex::WorkMem, 4 * 1024);

	// Set a large value, but mostly unused with NodeFileInd=1
	//cplex.setParam(IloCplex::NodeFileInd, 1);  // in MB

	// One thread per run to avoid memory spikes from multithreading
	cplex.setParam(IloCplex::Threads, 1);
		
	// Set the time limit to 30 minutes
	cplex.setParam(cplex.TiLim, 1800);
	
	// set number of threads	
	cplex.setParam(cplex.Threads, 1);

	// set the optimality gap
	cplex.setParam(IloCplex::EpGap, 0.00001);

	// set the time limit
	


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

// cplex turn off all cuts
void* Cplex::SetNoCuts() {
	cplex.setParam(cplex.Cliques, -1);
	cplex.setParam(cplex.Covers, -1);
	cplex.setParam(cplex.DisjCuts, -1);
	cplex.setParam(cplex.FlowCovers, -1);
	cplex.setParam(cplex.FlowPaths, -1);
	cplex.setParam(cplex.FracCuts, -1);
	cplex.setParam(cplex.GUBCovers, -1);
	cplex.setParam(cplex.ImplBd, -1);
	cplex.setParam(cplex.MIRCuts, -1);
	cplex.setParam(cplex.MCFCuts, -1);
	cplex.setParam(cplex.ZeroHalfCuts, -1);
	return this;
}

// cplex turn on default cuts
void* Cplex::SetCutsDefault() {
	cplex.setParam(cplex.Cliques, 0);
	cplex.setParam(cplex.Covers, 0);
	cplex.setParam(cplex.DisjCuts, 0);
	cplex.setParam(cplex.FlowCovers, 0);
	cplex.setParam(cplex.FlowPaths, 0);
	cplex.setParam(cplex.FracCuts, 0);
	cplex.setParam(cplex.GUBCovers, 0);
	cplex.setParam(cplex.ImplBd, 0);
	cplex.setParam(cplex.MIRCuts, 0);
	cplex.setParam(cplex.MCFCuts, 0);
	cplex.setParam(cplex.ZeroHalfCuts, 0);
	return this;
}


Cplex* Cplex::Run() {
	CplexSettings();

	if (cplex.solve()) {
		IloAlgorithm::Status status = cplex.getStatus();

		if (status == IloAlgorithm::Optimal) {
			gap = 0.0;
			opt = cplex.getObjValue();
		}
		else if (status == IloAlgorithm::Feasible) {
			gap = cplex.getMIPRelativeGap();
			opt = cplex.getObjValue();
		}
		else {
			std::cerr << "CPLEX returned unexpected status: " << status << std::endl;
			gap = -2;
			opt = -2;
		}
	}
	else {
		IloAlgorithm::Status status = cplex.getStatus();
		if (status == IloAlgorithm::Infeasible) {
			gap = -1;
			opt = -1;
		}
		else {
			std::cerr << "CPLEX failed with status: " << status << std::endl;
			gap = -3;
			opt = -3;
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