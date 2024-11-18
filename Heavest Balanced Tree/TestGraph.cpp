#define _CRT_SECURE_NO_WARNINGS

#include <Windows.h>
#include<direct.h>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

typedef boost::adjacency_list<boost::vecS, boost::setS, boost::undirectedS, boost::no_property, boost::property<boost::edge_weight_t, int>> Graph;

int main() {
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;

    // Random number generator setup
    boost::random::mt19937 gen(time(0));
    boost::random::uniform_real_distribution<> dist(0, 1);

    srand(12345);  // Seed for random number generator
    // Parameters

    int n[] = {5 ,10, 20, 40, 80, 160, 240, 320, 400};
    double p[] = {0.3, 0.4, 0.5, 0.7, 0.8};
    double k[] = { 2 , 3, 4, 5};

    char str[10000];
    // loop over parameters and print the instance in a string
    for (int i = 0; i < 9; i++) {
        for (int j = 0; j < 3; j++) {
            for (int l = 0; l < 4; l++) {
                int num_vertices = n[i];  // Number of vertices
                double edge_probability = p[j];  // Probability of edge existence
                int num_edges = edge_probability * (num_vertices) * (num_vertices - 1) / 2;  // Number of edges 
                int num_trees = k[l];  // Number of clusters

                // check if folder exists and if not create it 
                struct stat st = { 0 };
                if (stat("instances", &st) == -1) {
                    _mkdir("instances");
                }

                // generate file name 
                char filename[100];
                sprintf_s(filename, "instances//test%d_%d_%d.txt", num_vertices, num_edges, num_trees);


                // Open a file and write the string to the file
                FILE* f = fopen(filename, "w+");
                if (f == NULL) {
                    printf("Error opening file!\n");
                }
                else {


                    // Create a graph object
                    Graph g;

                    // Generate random graph
                    boost::generate_random_graph(g, num_vertices, num_edges, gen, false, false);

                    // Output the instance in a string
                    fprintf_s(f, "%d\t%d\t%d\n", num_vertices, num_edges, num_trees);

                    // Output edges
                    boost::graph_traits<Graph>::edge_iterator ei, eend;
                    for (boost::tie(ei, eend) = boost::edges(g); ei != eend; ++ei) {
                        int weight = rand() % 100 + 1;
                        fprintf_s(f, "%d\t %d \t %d\n", boost::source(*ei, g), boost::target(*ei, g), weight);
                    }
                }
                fclose(f);
            }
		}
	}
      
    return 0;

}
