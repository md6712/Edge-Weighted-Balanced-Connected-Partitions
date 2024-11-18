#define _CRT_SECURE_NO_WARNINGS


#include "_g.h"
#include "binary.h"
#include "CutF.h"
#include <stdio.h>
#include <conio.h>
#include "PartitionF.h"
#include "timer.h"

FILE* output;

void main()
{
    initBinary();
	
    int n[] = { 10, 20, 40, 80, 160, 240, 320, 400 };
    double p[] = { 0.3, 0.4, 0.5 };
    double k[] = { 2 , 3, 4, 5 };
    char nameoutput[80];
    char outputline[1000];

    char str[10000];
    // loop over parameters and print the instance in a string
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 3; j++) { 
            for (int l = 0; l < 4; l++) {

                output = fopen("output.txt", "a+");

                //compute number of edges 
                int num_vertices = n[i];  // Number of vertices
                int num_edges = p[j] * (num_vertices) * (num_vertices - 1) / 2;  // Number of edges
                int num_trees = k[l];  // Number of trees

                // read graph
                _g* g = new _g(num_vertices, num_edges, num_trees);
                g->readGraph();

                // run the algorithm
                Timer timer;
                timer.setStartTime();
                //(new CutF(g, false))->PrintModel()->Run();
                (new PartitionF(g, false))->PrintModel()->Run();
                timer.setEndTime();

                // print the results
                int opt = g->_opt;
                double gap = g->_gap;
                double eplasedtime = timer.calcElaspedTime_sec();

                sprintf_s(outputline, "%40s \t %d \t %d \t %d \t%d \t %4.2lf \t %4.2lf \n", g->getFilename(), num_vertices, num_edges, num_trees, opt, gap, eplasedtime);
                fprintf_s(output, "%s", outputline);
                printf("%s", outputline);
            }
        }
    }

	return;
}	