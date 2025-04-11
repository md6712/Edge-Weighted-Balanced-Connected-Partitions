#include "_mwcs.h"
#include "_g.h"

// constructor 
_mwcs::_mwcs(void* g) {
	this->instance = g;
	_g* graph = (_g*)g;
	this->num_vertices = graph->num_vertices + graph->num_edges;
	this->num_edges = 2 * graph->num_edges;
	this->vertex_weight = new double[this->num_vertices];
	this->vertex_map = new int[this->num_vertices];
	for (int i = 0; i < this->num_vertices; i++) {
		this->vertex_weight[i] = 0;
		this->vertex_map[i] = -1;
	}
	this->edges = new int[this->num_edges][2];

	for (int e = 0; e < graph->num_edges; e++) {
		this->vertex_map[graph->num_vertices + e] = e;
		this->edges[2 * e][0] = graph->edges[e][0];
		this->edges[2 * e][1] = graph->num_vertices + e;
		this->edges[2 * e + 1][0] = graph->edges[e][1];
		this->edges[2 * e + 1][1] = graph->num_vertices + e;
	}

	// init the string file
	this->stringfile = new char[1000000];


	// create the tree
	this->tree = new _small_tree();

	// set vertices in the tree
	this->vertices_in_tree = new int[this->num_vertices];
}

// destructor 
_mwcs::~_mwcs() {
	delete[] this->vertex_weight;
	delete[] this->vertex_map;
	delete[] this->edges;
	delete[] this->stringfile;
	delete[] this->vertices_in_tree;
	delete this->tree;
}

// set vertex weights
_mwcs* _mwcs::set_vertex_weights(double* vertex_weights, double theta) {
	memcpy(vertex_weight, vertex_weights, sizeof(double) * this->num_vertices);
	this->theta = theta;

	return this;
}

// solve 
_mwcs* _mwcs::solve() {

	// run the solver exe file in R://scipstp.exe	
	system("R://scipstp.exe -q -f R://mwcs.stp -s R://settings/write.set > nul 2>&1");

	return this;
}

// print 
_mwcs* _mwcs::print() {
	std::cout << "---- MWCS graph ----" << std::endl;
	std::cout << "Number of Vertices: " << this->num_vertices << std::endl;
	std::cout << "Number of Edges: " << this->num_edges << std::endl;
	std::cout << "Vertex Weights: ";
	for (int i = 0; i < this->num_vertices; i++) {
		std::cout << this->vertex_weight[i] << std::endl;
	}
	std::cout << std::endl;
	std::cout << "Edges: ";
	for (int e = 0; e < this->num_edges; e++) {
		std::cout << "(" << this->edges[e][0] << "," << this->edges[e][1] << ") " << std::endl;
	}
	std::cout << std::endl;
	return this;
}

// print in file
_mwcs* _mwcs::print_in_file() {

	// create header33
	// D32945 STP File, STP Format Version 1.0
	//SECTION Comment
	//Name	"test"
	//Creator	"CWI_LS_GWK"
	//Problem "Maximum Node Weight Connected Subgraph"
	//END

	// copy header to stringfile
	strcpy(this->stringfile, "33D32945 STP File, STP Format Version 1.0\n\n");
	strcat(this->stringfile, "SECTION Comments\n");
	strcat(this->stringfile, "Name\t\"test\"\n");
	strcat(this->stringfile, "Creator\t\"CWI_LS_GWK\"\n");
	strcat(this->stringfile, "Problem \"Maximum Node Weight Connected Subgraph\"\n");
	strcat(this->stringfile, "END\n\n");

	// copy Section Graph to stringfile
	strcat(this->stringfile, "SECTION Graph\n");
	sprintf(this->stringfile, "%sNodes %d\n", this->stringfile, this->num_vertices);
	sprintf(this->stringfile, "%sEdges %d\n", this->stringfile, this->num_edges);

	// copy edges to stringfile
	for (int e = 0; e < this->num_edges; e++) {
		sprintf(this->stringfile, "%sE %d %d\n", this->stringfile, this->edges[e][0]+1, this->edges[e][1]+1);
	}

	// End the section
	strcat(this->stringfile, "END\n\n");

	// copy Section Terminals to stringfile
	strcat(this->stringfile, "SECTION Terminals\n");
	sprintf(this->stringfile, "%sTerminals %d\n", this->stringfile, this->num_vertices);

	// copy terminals to stringfile
	for (int v = 0; v < this->num_vertices; v++) {		
		sprintf(this->stringfile, "%sT %d %4.6lf\n", this->stringfile, v+1, this->vertex_weight[v]);		
	}

	// End the section
	strcat(this->stringfile, "END\n\n");
		
	// End the file
	strcat(this->stringfile, "EOF\n");

	// write the file

		
	FILE* file;
	fopen_s(&file, "R:\\mwcs.stp", "w");
	fprintf(file, "%s", this->stringfile);
	fclose(file);


	// set start time
	//std::clock_t start = std::clock();

	// test file write speed
	
	/*for (int i = 0; i < 10000; i++) {
		fopen_s(&file, "R:\\mwcs.stp", "w");
		fprintf(file, "%s", this->stringfile);
		fclose(file);
	}*/
	
	// print time 
	//	std::cout << "Time: " << (std::clock() - start) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;

	return this;
}

// read solution from file
_mwcs* _mwcs::read_solution_from_file() {
	// read the file
	FILE* file;
	fopen_s(&file, "R:\\mwcs.stp.sol", "r");
	if (file == NULL) {
		std::cout << "Error reading file" << std::endl;
		return this;
	}
	// read the file
	char line[1000];
	char* token; 
	while (fgets(line, 1000, file) != NULL) {
		// check if the line contains Primal then the number in front is the objective value
		if (strstr(line, "Primal") != NULL) {
			token = strtok(line, " ");
			token = strtok(NULL, " ");	
			objective_value = atof(token) - theta;
			break; 
		}
	}

	int num_vertices_in_tree = 0;

	while (fgets(line, 1000, file) != NULL) {
		// check if the line contains Vertices then the number in front is the number of vertices in the tree
		if (strstr(line, "Vertices") != NULL) {
			token = strtok(line, " ");
			token = strtok(NULL, " ");
			num_vertices_in_tree = atoi(token);
			break;
		}
	}

	// read the vertices in the tree
	int v = 0;
	while (fgets(line, 1000, file) != NULL) {
		if (strstr(line, "V") != NULL) {
			token = strtok(line, " ");
			token = strtok(NULL, " ");
			vertices_in_tree[v++] = atoi(token) - 1;
		}		
	}

	// print objective value
	std::cout << "Objective Value: " << objective_value << "\t";


	// graph
	_g* graph = (_g*)instance;

	// print vertices in the tree
	std::cout << "Vertices in the tree: ";
	for (int i = 0; i < v; i++) {
		if (vertices_in_tree[i] < graph->num_vertices)
			std::cout << vertices_in_tree[i] << " ";
	}

	std::cout << std::endl;

	fclose(file);
	return this;
}


// pipes ready
_mwcs* _mwcs::pipes_ready() {
		
	// convert name
	wchar_t wname[50];
	mbstowcs(wname, pipeName_input, strlen(pipeName_input) + 1);

	// step 1: create named pipe for input
	hNamedPipe_input = CreateNamedPipe(
		wname,						// pipe name 
		PIPE_ACCESS_OUTBOUND,       // write access 
		PIPE_TYPE_BYTE |			// byte type pipe 
		PIPE_READMODE_BYTE |		// byte-read mode 
		PIPE_WAIT,					// blocking mode 
		1,							// number of instances 
		1024,						// output buffer size 
		1024,						// input buffer size 
		0,							// client time-out 
		NULL);						// default security attribute

	// check if handle is okay
	if (hNamedPipe_input == INVALID_HANDLE_VALUE) {
		std::cout << "Error creating named pipe for input" << std::endl;
		return this;
	}

	// convert name output
	mbstowcs(wname, pipeName_output, strlen(pipeName_output) + 1);

	// step 2: create named pipe for output
	hNamedPipe_output = CreateNamedPipe(
		wname,						// pipe name 
		PIPE_ACCESS_INBOUND,        // read access 
		PIPE_TYPE_BYTE |			// byte type pipe 
		PIPE_READMODE_BYTE |		// byte-read mode 
		PIPE_WAIT,					// blocking mode 
		1,							// number of instances 
		1024,						// output buffer size 
		1024,						// input buffer size 
		0,							// client time-out 
		NULL);						// default security attribute

	//check if handle is okay
	if (hNamedPipe_output == INVALID_HANDLE_VALUE) {
		std::cout << "Error creating named pipe for output" << std::endl;
		return this;
	}

	return this;
}




// pipe solve
_mwcs* _mwcs::pipes_solve() {
	
	// security alert setting
	saAttr.nLength = sizeof(SECURITY_ATTRIBUTES);
	saAttr.bInheritHandle = TRUE;
	saAttr.lpSecurityDescriptor = NULL;

	// set command
	std::string command = "R://scipstp.exe -q -f" + 
		std::string(pipeName_input) + " -l " +
		std::string(pipeName_output) + " > nul 2>&1";

	// create the process startup info
	STARTUPINFO si;
	PROCESS_INFORMATION pi;
	ZeroMemory(&si, sizeof(STARTUPINFO));
	si.cb = sizeof(STARTUPINFO);
	si.dwFlags |= STARTF_USESTDHANDLES;

	// convert command
	wchar_t wcommand[1000];
	mbstowcs(wcommand, command.c_str(), strlen(command.c_str()) + 1);

	// create the process
	if (!CreateProcess(NULL, wcommand, NULL, NULL, TRUE, 0, NULL, NULL, &si, &pi)) {
		std::cerr << "Failed to start scipstp.exe." << std::endl;
		CloseHandle(hNamedPipe_input);
		CloseHandle(hNamedPipe_output);
		return this;
	}

	// wait for the solver to open the pipes
	std::cout << "Waiting for solver to open pipes..." << std::endl;
	if (!ConnectNamedPipe(hNamedPipe_input, NULL)) {
		std::cerr << "Solver failed to connect to input pipe." << std::endl;
		CloseHandle(hNamedPipe_input);
		CloseHandle(hNamedPipe_output);
		return this;
	}
	if (!ConnectNamedPipe(hNamedPipe_output, NULL)) {
		std::cerr << "Solver failed to connect to output pipe." << std::endl;
		CloseHandle(hNamedPipe_input);
		CloseHandle(hNamedPipe_output);
		return this;
	}

	return this;

	//LPDWORD bytesWritten;

	//// write input to pipe
	//if (!WriteFile(hNamedPipe_input, this->stringfile, strlen(this->stringfile), &bytesWritten, NULL)) {
	//	std::cerr << "Failed to write STP data to input pipe." << std::endl;
	//}
	//CloseHandle(hNamedPipe_input); // Close pipe to signal EOF
	

	// I gave up this idea.... 
	// remains to be completed in case we need to use PIPE with Daniel's code 


}