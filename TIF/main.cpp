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

int n[] = { 10, 15, 20, 25, 30, 35, 40, 45, 50 };
double p[] = { 0.1, 0.2, 0.3, 0.4, 0.5 };
int k[] = { 2, 3, 4, 5, 6, 7, 8, 9, 10 };

enum RunType {
    t_SETCOVER = 1,
    t_CUT = 2,
    t_DiCUT= 3,
	t_FLOW = 4,
    t_EX = 5,
    t_CG = 6,
    t_BB = 7
};

struct _settings {
    RunType model;
    int vertex_size;
    int edge_prob;
    int num_trees;
	int duplication_num;
    bool read;
    bool linear;
	bool flow_user_cut;
	bool flow_mutual_exclusion_cycles;
    int min_size;
    int max_size;
    int min_prob;
    int max_prob;
    int min_trees;
    int max_trees;
	int min_duplication;
	int max_duplication;	
    int BB_Heuristic_Time_Limit;
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

char* nameOfModel(RunType r) {
    char* name = new char[180];
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
    return name; 
}

// create an output file name: 
char* outputFileName(_settings s) {

	int len = 280;

	char* name = nameOfModel(s.model);
	char* fullname = new char[len];
    
	// attach a timestamp
	time_t now = time(0);
	tm* ltm = localtime(&now);
	char timestamp[80];
	sprintf_s(timestamp, len, "%d_%d_%d_%d_%d_%d", 1900 + ltm->tm_year, 1 + ltm->tm_mon, ltm->tm_mday, ltm->tm_hour, ltm->tm_min, ltm->tm_sec);
	char setting[100];
	sprintf_s(setting, len, "%d-%d_%2.1lf-%2.1lf_%d-%d_%s", n[s.min_size], n[s.max_size-1], p[s.min_prob], p[s.max_prob-1], k[s.min_trees], k[s.max_trees-1], s.linear ? "linear" : "integer");
    strcat_s(name, len, "_");
    strcat_s(name, len, setting);
	strcat_s(name, len, "_");
	strcat_s(name, len, timestamp);
	strcat_s(name, len, ".txt");

	sprintf_s(fullname, len, "output\\%s", name);
	
	delete name;
	return fullname;
}

// read setting from file
_settings readSettings() {
	// read file called settings.txt
	
    FILE* f = fopen("setting.txt", "r");

	char str[10000];

	if (f == NULL) {
		printf("Error opening file!\n");

        _settings s; 
		s.read = false;
		return s;
	}
	else {        

        // read model
		int model;
		fgets(str, 10000, f);
		fscanf(f, "%d", &model);

		// read vertex size
		int vertex_size;
        fgets(str, 10000, f);
		fgets(str, 10000, f);
        fgets(str, 10000, f);        
		fscanf(f, "%d", &vertex_size);

		// read edge probability
		int edge_prob;
		fgets(str, 10000, f);
		fgets(str, 10000, f);
        fgets(str, 10000, f);
		fscanf(f, "%d", &edge_prob);

		// read num trees
		int num_trees;
		fgets(str, 10000, f);
		fgets(str, 10000, f);
        fgets(str, 10000, f);
		fscanf(f, "%d", &num_trees);

		// read duplication number
		int duplication_num;
		fgets(str, 10000, f);
		fgets(str, 10000, f);
		fgets(str, 10000, f);
		fscanf(f, "%d", &duplication_num);

		RunType type = RunType::t_CUT;
        if (model == 1) {
            type = RunType::t_CUT;
		}
		else if (model == 2) {			
            type = RunType::t_FLOW;
		}
		else if (model == 3) {
            type = RunType::t_SETCOVER;
		}        
        else if (model == 4) {
			type = RunType::t_CG;
        }

        // read linear 
		fgets(str, 10000, f);
		fgets(str, 10000, f);
		fgets(str, 10000, f);
		int linear_int;
		fscanf(f, "%d", &linear_int);


		// read flow user cut
        fgets(str, 10000, f);
        fgets(str, 10000, f);
        fgets(str, 10000, f);
        int flow_user_cut_int;
        fscanf(f, "%d", &flow_user_cut_int);

		// read flow mutual exclusion cycles
		fgets(str, 10000, f);
		fgets(str, 10000, f);
		fgets(str, 10000, f);
		int flow_mutual_exclusion_cycles_int;
		fscanf(f, "%d", &flow_mutual_exclusion_cycles_int);


		// read bb heuristic time limit
		int bb_heuristic_time_limit;
		fgets(str, 10000, f);
		fgets(str, 10000, f);
		fgets(str, 10000, f);
		fscanf(f, "%d", &bb_heuristic_time_limit);
		
		// linear
        bool linear = false;
        if (linear_int > 0) {
            linear = false;
        }
        else {
            linear = true;
        }

		// flow user cut
		bool flow_user_cut = false; 
        if (flow_user_cut_int > 0) {
            flow_user_cut = true;
        }   
		else {
			flow_user_cut = false;
		}

		// flow mutual exclusion cycles
		bool flow_mutual_exclusion_cycles = false;
		if (flow_mutual_exclusion_cycles_int > 0) {
			flow_mutual_exclusion_cycles = true;
		}
		else {
			flow_mutual_exclusion_cycles = false;
		}

        // setting created
        _settings settings;
        settings.model = type;
        settings.vertex_size = vertex_size;
        settings.edge_prob = edge_prob;
        settings.num_trees = num_trees;
        settings.read = true;
		settings.BB_Heuristic_Time_Limit = bb_heuristic_time_limit;


        settings.linear = linear;
		settings.flow_user_cut = flow_user_cut;
		settings.flow_mutual_exclusion_cycles = flow_mutual_exclusion_cycles;

        if (settings.vertex_size == 0) {
            settings.min_size = 0;
            settings.max_size = 9;
        }
        else {
            settings.min_size = settings.vertex_size- 1;
            settings.max_size = settings.vertex_size;
        }

        if (settings.edge_prob == 0) {
            settings.min_prob = 0;
            settings.max_prob = 5;
        }
        else {
            settings.min_prob = settings.edge_prob - 1;
            settings.max_prob = settings.edge_prob;
        }

        if (settings.num_trees == 0) {
            settings.min_trees = 0;
            settings.max_trees = 9;
        }
        else {
            settings.min_trees = settings.num_trees - 1;
            settings.max_trees = settings.num_trees;
        }

        if (duplication_num == 0) {
            settings.min_duplication = 1;
            settings.max_duplication = 5;
        }
		else {
			settings.min_duplication = duplication_num;
			settings.max_duplication = duplication_num + 1;
		}

        // string for vertex range
		char vertex_size_text[100];
        if (settings.min_size + 1 == settings.max_size) {
			sprintf_s(vertex_size_text, "%d", n[settings.min_size]);
        }
        else
        {
			sprintf_s(vertex_size_text, "%d - %d", n[settings.min_size], n[settings.max_size-1]);
        }

		char edge_prob_text[100];
		if (settings.min_prob + 1 == settings.max_prob) {
			sprintf_s(edge_prob_text, "%2.1lf", p[settings.min_prob]);
		}
		else
		{
			sprintf_s(edge_prob_text, "%2.1lf - %2.1lf", p[settings.min_prob], p[settings.max_prob-1]);
		}

		char num_trees_text[100];
		if (settings.min_trees + 1 == settings.max_trees) {
			sprintf_s(num_trees_text, "%d", k[settings.min_trees]);
		}
		else
		{
			sprintf_s(num_trees_text, "%d - %d", k[settings.min_trees], k[settings.max_trees-1]);
		}
		


		// print settings
		char* name = nameOfModel(type);
		printf("Model: %s\n", name);
		delete name;

        printf("Linear: %s\n", linear ? "true" : "false");

		printf("Vertex size: %s\n", vertex_size_text);
		printf("Edge probability: %s\n", edge_prob_text);
		printf("Num trees: %s\n", num_trees_text);
		
		fclose(f);
		return settings;
	}    
}

int main()
{    
	_settings settings = readSettings();
    
    initBinary();
   
    char nameoutput[80];
    char outputline[1000];

    char str[10000];

	const char *outputname = outputFileName(settings);
    output = fopen(outputname, "w");
	fclose(output);

    // loop over parameters and print the instance in a string
    for (int i = settings.min_size; i < settings.max_size; i++) {
        int num_vertices = n[i];  // Number of vertices
        for (int j = settings.min_prob; j < settings.max_prob; j++) {
            double edge_probability = p[j];  // Probability of edge existence
            int num_edges = edge_probability * (num_vertices) * (num_vertices - 1) / 2;  // Number of edges 
            
            if (num_edges < num_vertices) {
                continue;  // Skip if number of edges is less than number of vertices
            }

            for (int l = settings.min_trees; l < settings.max_trees; l++) {
                int num_trees = k[l];  // Number of clusters

                if ((double)num_vertices / 2 < num_trees) {
                    continue;  // Skip if number of trees is greater than number of vertices
                }

                for (int ii = settings.min_duplication; ii < settings.max_duplication; ii++) {
                    // timer
                    Timer timer;
                    timer.setStartTime();

					// open the output file
                    output = fopen(outputname, "a+");
                             
                    // initialize the graph
                    _g* g = (new _g(num_vertices, num_edges, num_trees, ii));

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
                    //g->printGraph();

                    // compute the upper bound
                    //g->computeUB();

                    // compute the lower bound
                    g->recomputeLB();

                    // reset g ub
					g->UB = INT32_MAX;                    

                    // compute ub branching
                    Timer timer_ub;
                    timer_ub.setStartTime();

                    bb* model_1 = (new bb(g))
                        ->Run(settings.BB_Heuristic_Time_Limit);
                    delete model_1;

                    timer_ub.setEndTime();
					double elapsedtime_ub = timer_ub.calcElaspedTime_sec();                    

                    // populate ub trees
                    g->populate_trees_ub_from_select_trees();

                    // compute coords
                    //g->ForcedDirectedLayout();

                    //  print the graph
                    //g->DrawGraph();

                    // print min separators
                    //g->PrintMinSeperators();

                   
                    // shortest path weight
					g->compute_all_shortest_paths_weight();
                                     

					// timer for lower bound
                    Timer timer_lb;                   		
					
					// timer for model
					Timer timer_model;					

					double elapsedtime_lb = 0;  
					double elapsedtime_model = 0;



                    if (settings.model == RunType::t_CUT) {
						timer_lb.setStartTime();
						if (!settings.linear) {							
                            CutF* model = (new CutF(g, false, true))                                
                                ->Run();
                            delete model;

							g->_lp_bound = g->_opt;
						}
						timer_lb.setEndTime();
						elapsedtime_lb = timer_lb.calcElaspedTime_sec();

						timer_model.setStartTime();
                        CutF* model = (new CutF(g, false, settings.linear)); 
                        if (!settings.linear) {
                            model->SetInitSol();
                        }
                        //model->ForceSol();
                        //model->PrintModel();
                        model->Run();
                        delete model;
						timer_model.setEndTime();
						elapsedtime_model = timer_model.calcElaspedTime_sec();
                    }

                    /*if (settings.model == RunType::t_DiCUT) {
                        CutF* model = (new CutF(g, false, settings.linear))
                            ->Run();
                        delete model;
                    }*/

                    if (settings.model == RunType::t_CG) {

                        g->createMWCS();  // create the MWCS instance
                        g->createPCST();  // create the PCST instance                                               

                        // generate more trees
						g->generateTreesForCG();

                        CG* model = ((CG*)(new CG(g, false, settings.linear))
                            ->SetQuiet()
                            )
                            ->Run_BP();
                        delete model;
                    }

                    if (settings.model == RunType::t_FLOW) {
                        timer_lb.setStartTime();
                        if (!settings.linear) {
                            FlowF* model = (new FlowF(g, false, true, settings.flow_mutual_exclusion_cycles));
                            model->SetUserCutsActive(settings.flow_user_cut);
                            model->Run();
                            delete model;

                            g->_lp_bound = g->_opt;
                        }
                        timer_lb.setEndTime();
                        elapsedtime_lb = timer_lb.calcElaspedTime_sec();

                        timer_model.setStartTime();
                        FlowF* model = (new FlowF(g, false, settings.linear, settings.flow_mutual_exclusion_cycles));
                        if (!settings.linear) {
                            model->SetInitSol();                            
                        }
						model->SetUserCutsActive(settings.flow_user_cut);						

                        model->Run();
      
                        delete model;
                        timer_model.setEndTime();
                        elapsedtime_model = timer_model.calcElaspedTime_sec();
                    }

                    if (settings.model == RunType::t_SETCOVER) {                      
                        g->generateTrees();

                        timer_lb.setStartTime();
                        if (!settings.linear) {
                            SetCoverF* model = (new SetCoverF(g, false, true))
                                ->Run();
                            delete model;

                            g->_lp_bound = g->_opt;
                        }
                        timer_lb.setEndTime();
                        elapsedtime_lb = timer_lb.calcElaspedTime_sec();

                        timer_model.setStartTime();
                        SetCoverF* model = (new SetCoverF(g, false, settings.linear));
                        if (!settings.linear) {
							model->SetInitSol();
						}
                        model->Run();
                        delete model;
                        timer_model.setEndTime();
						elapsedtime_model = timer_model.calcElaspedTime_sec();
                    }

                    if (settings.model == RunType::t_EX) {
                        ExF* model = (new ExF(g, false, settings.linear))
                            //->PrintModel()
                            ->Run();
                        delete model;
                    }                   

                    //PartitionF* model = (new PartitionF(g, false))->Run();
                    timer.setEndTime();

                    // print the results
                    double opt = g->_opt + FLT_EPSILON;
                    double gap = g->_gap;
                    double elapsedtime = timer.calcElaspedTime_sec();

					// number of trees depend on the model; for SetCover it is the number of trees in the instance ; for CG it is the number of small trees generated
					int number_of_trees = (int)g->select_trees_for_CG.size();
                    if (settings.model == RunType::t_SETCOVER)
					{
						number_of_trees = g->trees.size();
					}

                    sprintf_s(outputline, 512,
                        "%40s\t%d\t%d\t%d\t%d\t%d\t%d\t%6.2lf\t%6.2lf\t%6.2lf\t%6.2lf\t%6.2lf\t%6.2lf\t%6.2lf\t%6.2lf\t%d\t%d\t%d\t\n",
                        g->getFilename(),       // %40s
                        num_vertices,           // %d
                        num_edges,              // %d
                        num_trees,              // %d                        
                        ii,
                        g->UB_naive,                  // %d
                        g->UB,                  // %d
						g->_lp_bound,           // %6.2lf
                        g->_lb_root,            // %6.2lf
						g->_opt,               // %6.2lf
						g->_gap,               // %6.2lf                                                
						elapsedtime_lb,         // %6.2lf
						elapsedtime_ub,         // %6.2lf
                        elapsedtime_model,      // %6.2lf
                        elapsedtime,            // %6.2lf
                        number_of_trees,    // %d                        
						g->n_user_cuts,         // %d
						g->n_lazy_cuts         // %d
                    );
                    fprintf_s(output, "%s", outputline);
                    printf("%s", outputline);					

                    fclose(output);

                    //g->PrintOptEdges();

                    delete g;
                }
                
            }
        }
    }

    delete outputname;

	return 0;
}	