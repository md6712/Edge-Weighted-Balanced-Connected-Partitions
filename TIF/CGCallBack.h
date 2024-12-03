// Add column generation callback to CG class
#include <ilcplex/ilocplex.h>
#include <iostream>
#include "CG.h"

ILOSTLBEGIN

class CGColumnGenCallback : public IloCplex::Callback::Function {
public:
    // Constructor to pass necessary data
    CGColumnGenCallback(CG *cg)
        : env(cg->env), obj(cg->obj), Cover(cg->Cover), x(cg->x) {}

    // Override the function operator to define callback logic    

    void invoke(const IloCplex::Callback::Context& context) override {
        // Check if in the appropriate context
        if (!context.inCandidate()) return;

        IloEnv env = context.getEnv();

        IloNumArray reducedCosts(env, x.getSize());        

        //for (int i = 0; i < x.getSize(); ++i) {
        //    if (reducedCosts[i] < -1e-6) { // A small threshold to identify negative reduced costs
        //        // Logic to generate a new column
        //        AddColumn(i, context);
        //    }
        //}

        //reducedCosts.end();
    }

private:
    IloEnv env;
    _g* instance;
    IloObjective obj;
    IloRangeArray Cover;
    IloNumVarArray x;

    //// Function to add a new column
    //void AddColumn(int index, const IloCplex::Callback::Context& context) {
    //    IloEnv env = context.getEnv();
    //    IloNumVar newVar(env, 0, IloInfinity, ILOFLOAT);

    //    // Define the column coefficients (example logic, adjust for your model)
    //    IloNumColumn col = obj(0); // Add to objective
    //    //for (int v = 0; v < instance->num_vertices; v++) {
    //    //    if (checkbin(instance->trees[index]->bin_vertices, v)) {
    //    //        col += CoverUpdate constraint coefficients
    //    //    }
    //    //    else {
    //    //        col += Cover;
    //    //    }
    //    //}

    //    //// Add the new column (variable) to the model
    //    //newVar = col;
    //    //x.add(newVar);
    //    //context.add(newVar);

    //    printf("Added new column for index %d\n", index);
    //}
};