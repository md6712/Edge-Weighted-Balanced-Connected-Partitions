#pragma once
#include <SFML/Graphics.h>

// using std vector
#include <vector>

enum _render_object_type
{
	shape,
	text, 
	line
};

class _render_object
{
	

public:
	void* object;
	_render_object_type type;
	sfVector2f position;
};

class _render
{

	sfRenderWindow* window;
	sfVideoMode mode = { 1200, 800, 32 };

public:
	// vector of objects
	std::vector<_render_object *> objects;

	// constructor & destructor
	_render();
	~_render();

	void draw();

	// add a node
	void addNode(int number, int x,  int y);

	
	// add an edge
	void addEdge(int x_from, int y_from, int x_to, int y_to, sfPrimitiveType type);

	// highlight an edge	
	void highlightEdge(int x_from, int y_from, int x_to, int y_to);

	// add think edge
	void addThickEdge(sfVertexArray* vertexArray, sfVector2f start, sfVector2f end, float thickness, sfColor color);
};

