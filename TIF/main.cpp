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
#include "FlowF.h"
#include "CG.h"
#include "bb.h"

#include "_tree.h"

#include "_render.h"

FILE* output;

enum RunType {
    t_SETCOVER = 1,
    t_CUT = 2,
    t_DiCUT= 3,
	t_FLOW = 4,
    t_EX = 5,
    t_CG = 6,
    t_BB = 7
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


// create an output file name: 
char* outputFileName(RunType r) {

    char* name = new char[180];
	char* fullname = new char[180];
	
    
        switch (r) {
        case RunType::t_SETCOVER:            
            sprintf_s(name, 180, "SetCover");
            break;
        case RunType::t_CUT:
            sprintf_s(name, 180, "Cut");
            break;
		case RunType::t_DiCUT:
			sprintf_s(name, 180, "DiCut");
			break;
        case RunType::t_FLOW:
            sprintf_s(name, 180, "Flow");
            break;
        case RunType::t_EX:
            sprintf_s(name, 180, "Ex");
            break;
        case RunType::t_CG:
            sprintf_s(name, 180, "CG");
            break;
        case RunType::t_BB:
            sprintf_s(name, 180, "BB");
            break;
        }
    
	// attach a timestamp
	time_t now = time(0);
	tm* ltm = localtime(&now);
	char timestamp[80];
	sprintf_s(timestamp, 180, "%d_%d_%d_%d_%d_%d", 1900 + ltm->tm_year, 1 + ltm->tm_mon, ltm->tm_mday, ltm->tm_hour, ltm->tm_min, ltm->tm_sec);    
	strcat_s(name, 180, timestamp);
	strcat_s(name, 180, ".txt");

	sprintf_s(fullname, 180, "output\\%s", name);
	
	delete name;

	return fullname;
}


void main()
{    

	//_render* render = new _render();    
 //   
 //   // add lines
	//render->addEdge(100, 100, 200, 200);
	//render->addEdge(200, 200, 300, 300);

 //   // add tree nodes
 //   render->addNode(1, 100, 100);
 //   render->addNode(2, 200, 200);
	//render->addNode(3, 300, 300);

	//render->start();

    RunType runType = RunType::t_BB;
    bool linear = false;

    initBinary();

    int* a = new int[10];
	
    int n[] = { 5, 10, 20, 40, 80, 160, 240, 320, 400 };
    double p[] = { 0.3, 0.4, 0.5 };
    double k[] = { 2 , 3, 4, 5 };
    char nameoutput[80];
    char outputline[1000];

    char str[10000];

	const char *outputname = outputFileName(runType);
    output = fopen(outputname, "w");
	fclose(output);

    // loop over parameters and print the instance in a string
    for (int i = 1; i < 2; i++) {
        for (int j = 2; j < 3; j++) { 
            for (int l = 0; l < 4; l++) {

                output = fopen(outputname, "a+");

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

				// compute the directed graph
                g->setDiGraph();

                // print the graph            
                g->printGraph();

                // compute the upper bound
                g->computeUB();

				// compute the lower bound
				g->recomputeLB();

                // compute coords
				g->ForcedDirectedLayout();

				//  print the graph
                g->DrawGraph();

                // print min separators
                //g->PrintMinSeperators();

                // generate all trees 

                // run the algorithm
                Timer timer;
                timer.setStartTime();



                if (runType == RunType::t_CUT) {
                    CutF* model =  (new CutF(g, false, linear))
                        //->ForceSol()
                        //->SetPrintCuts(true)
                        //->PrintModel()                    
                    ->Run();
                    delete model;
                }

				if (runType == RunType::t_DiCUT) {
					CutF* model = (new CutF(g, false, linear))
						->Run();
					delete model;
				}
                
                if (runType == RunType::t_CG) {
                    //g->generateTrees();

                   g->generateSelectTrees();
                    CG* model = (new CG(g, false, linear))
                        ->PrintModel()
                        ->Run()
					    ->PrintSol();
                    delete model;
                }

				if (runType == RunType::t_FLOW) {
                    FlowF* model = (new FlowF(g, false, linear))
                        //->PrintModel()                        
                        ->Run()
                        ->PrintSol();
					delete model;
				}

                if (runType == RunType::t_SETCOVER) {
                    g->generateTrees();
					SetCoverF* model = (new SetCoverF(g, false, linear))
                        //->PrintModel()						
						->Run();
                    delete model;
				}

                if (runType == RunType::t_EX) {
					ExF* model = (new ExF(g, false, linear))
						//->PrintModel()
						->Run();
                    delete model;
				}  

                // enumeration based algorithms 
                if (runType == RunType::t_BB) {
                    bb* model = (new bb(g))
                        ->Run();
                }

                //PartitionF* model = (new PartitionF(g, false))->Run();
                timer.setEndTime();

                // print the results
                double opt = g->_opt + FLT_EPSILON;
                double gap = g->_gap;
                double eplasedtime = timer.calcElaspedTime_sec();

                sprintf_s(outputline, "\n%40s \t %d \t %d \t %d \t%d \t%d \t%4.2lf \t %4.2lf \t %4.2lf  \t %2d", g->getFilename(), num_vertices, num_edges, num_trees, g->LB, g->UB, opt, gap, eplasedtime, g->trees.size());
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