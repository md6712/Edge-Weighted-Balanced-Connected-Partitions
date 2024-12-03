#define _CRT_SECURE_NO_WARNINGS


#include "_g.h"
#include "binary.h"
#include <stdio.h>
#include <conio.h>
#include "PartitionF.h"
#include "timer.h"
#include <iostream>


#include "SetCoverF.h"
#include "CutF.h"
#include "ExF.h"
#include "CG.h"

#include "_tree.h"


FILE* output;

enum RunType {
    t_SETCOVER = 1,
    t_CUT = 2,
    t_EX = 3,
    t_CG = 4
};

void testTree(_g* instance) {
    _tree* tree = new _tree(instance);
    tree->AddEdge(0);
    tree->AddEdge(1);
    tree->AddEdge(2);
    tree->AddEdge(3);
    tree->AddEdge(3);
    tree->AddEdge(4);
	

    tree->PrintTree();
	delete tree;
}

void main()
{    

    RunType runType = RunType::t_CG;

    initBinary();

    int* a = new int[10];
	
    int n[] = { 5, 10, 20, 40, 80, 160, 240, 320, 400 };
    double p[] = { 0.3, 0.4, 0.5 };
    double k[] = { 2 , 3, 4, 5 };
    char nameoutput[80];
    char outputline[1000];

    char str[10000];
    // loop over parameters and print the instance in a string
    for (int i = 1; i < 4; i++) {
        for (int j = 1; j < 3; j++) { 
            for (int l = 1; l < 4; l++) {

                output = fopen("output.txt", "a+");

                //compute number of edges 
                int num_vertices = n[i];  // Number of vertices
                int num_edges = p[j] * (num_vertices) * (num_vertices - 1) / 2;  // Number of edges
                int num_trees = k[l];  // Number of trees

                // initialize the graph
                _g* g = (new _g(num_vertices, num_edges, num_trees));

                // read the graph
                g->readGraph();

                //testTree(g);

                // compute all subsets
                //g->computeSubsets();             

                // sort the edges
                g->SortEdges();

                g->setDiGraph();

                // print the graph            
                g->printGraph();

                // compute the upper bound
                g->computeUB();

                // print min separators
                //g->PrintMinSeperators();

                // generate all trees 

                // run the algorithm
                Timer timer;
                timer.setStartTime();



                if (runType == RunType::t_CUT) {
                    CutF* model =  (new CutF(g, false))
                        //->ForceSol()
                        //->SetPrintCuts(true)
                        //->PrintModel()                    
                    ->Run();
                    delete model;
                }
                
                if (runType == RunType::t_CG) {
                    //g->generateTrees();

                   g->generateSelectTrees();
                    CG* model = (new CG(g, false, false))
                        ->PrintModel()
                        ->Run()
					    ->PrintSol();
                    delete model;
                }

                if (runType == RunType::t_SETCOVER) {
                    g->generateTrees();
					SetCoverF* model = (new SetCoverF(g, false, false))				
                        //->PrintModel()
						->Run();
                    delete model;
				}

                if (runType == RunType::t_EX) {
					ExF* model = (new ExF(g, false, false))
						//->PrintModel()
						->Run();
                    delete model;
				}                

                //PartitionF* model = (new PartitionF(g, false))->Run();
                timer.setEndTime();

                // print the results
                double opt = g->_opt + FLT_EPSILON;
                double gap = g->_gap;
                double eplasedtime = timer.calcElaspedTime_sec();

                sprintf_s(outputline, "\n%40s \t %d \t %d \t %d \t%4.2lf \t %4.2lf \t %4.2lf ", g->getFilename(), num_vertices, num_edges, num_trees, opt, gap, eplasedtime);
                fprintf_s(output, "%s", outputline);
                printf("%s", outputline);   

                fclose(output);

                g->PrintOptEdges();
                //g->outputOPTEdges();

                delete g;
                
            }
        }
    }

   // OutputDebugString(L"-----------_CrtDumpMemoryLeaks ---------");
    _CrtDumpMemoryLeaks();

	return;
}	