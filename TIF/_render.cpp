#include "_render.h"


// constructor
_render::_render() {
	objects = std::vector<_render_object*>();

	// init the vector of objects
	
}

// destructor 
_render::~_render() {
	
	objects.end();
}

// start the render
void _render::draw() {
	
	window = sfRenderWindow_create(mode, "SFML window", sfResize | sfClose, NULL);
	while (sfRenderWindow_isOpen(window))
	{
		sfEvent event;
		while (sfRenderWindow_pollEvent(window, &event))
		{
			if (event.type == sfEvtClosed)
				sfRenderWindow_close(window);
		}
		sfRenderWindow_clear(window, sfWhite);
		
		// draw all objects
		for (int i = 0; i < objects.size(); i++) {
			if (objects[i]->type == _render_object_type::shape) {
				sfCircleShape* shape = (sfCircleShape*)objects[i]->object;
				sfRenderWindow_drawCircleShape(window, shape, NULL);
			}
			else if (objects[i]->type == _render_object_type::text) {
				sfText* text = (sfText*)objects[i]->object;
				sfRenderWindow_drawText(window, text, NULL);
			}
			else if (objects[i]->type == _render_object_type::line) {
				sfVertexArray* lines = (sfVertexArray*)objects[i]->object;				
				sfRenderWindow_drawVertexArray(window, lines, NULL);
			}
		}

		sfRenderWindow_display(window);
	}
	
	// destroy all objects
	for (int i = 0; i < objects.size(); i++) {
		if (objects[i]->type == _render_object_type::shape) {
			sfCircleShape_destroy((sfCircleShape*)objects[i]->object);
		}
		else if (objects[i]->type == _render_object_type::text) {
			sfText_destroy((sfText*)objects[i]->object);
		}
		else if (objects[i]->type == _render_object_type::line) {
			sfVertexArray_destroy((sfVertexArray*)objects[i]->object);
		}
		delete objects[i];
	}	

	sfRenderWindow_destroy(window);	

	objects.clear();
}


// add a node
void _render::addNode(int number, int x, int y) {
	
	// radius of the node
	int radius = 20;

	// create a label based on number
	char label[20];
	sprintf_s(label, "%d", number);

	// compute approximate width and height of the label with font size 24
	int width = 24 * strlen(label);
	int height = 24;	
	
	// Create a shape
	sfCircleShape* shape = sfCircleShape_create();
	sfCircleShape_setRadius(shape, radius);
	
	// position of node
	sfVector2f vector;
	vector.x = x - radius;
	vector.y = y - radius;

	// position of label
	sfVector2f vector_label;
	vector_label.x = x - width / 2;
	vector_label.y = y - height / 2;
	
	sfCircleShape_setPosition(shape, vector);
	sfCircleShape_setFillColor(shape, sfGreen);
	
	sfFont* font = sfFont_createFromFile("arial.ttf");
	
	sfText* text = sfText_create();
	sfText_setFont(text, font);
	sfText_setCharacterSize(text, 24);  // Set font size
	sfText_setString(text, label);  // Set the text
	sfText_setPosition(text, vector_label);


	// create render_objects
	_render_object* obj = new _render_object();
	obj->object = (void*) shape;
	obj->type = _render_object_type::shape;
	obj->position = vector;

	// create render_objects: text	
	_render_object* obj_text = new _render_object();
	obj_text->object = (void*)text;
	obj_text->type = _render_object_type::text;
	obj_text->position = vector;


	// add the object to the vector
	objects.push_back(obj);
	objects.push_back(obj_text);	

}

// add an edge
void _render::addEdge(int x_from, int y_from, int x_to, int y_to, sfPrimitiveType type) {

	// create a line
	sfVertexArray* lines = sfVertexArray_create();

	// set the primitive type
	sfVertexArray_setPrimitiveType(lines, type);

	// create a line
	sfVertex line[2];
	line[0].position.x = x_from;
	line[0].position.y = y_from;
	line[1].position.x = x_to;
	line[1].position.y = y_to;

	sfVertexArray_append(lines, line[0]);
	sfVertexArray_append(lines, line[1]);
	// create render_objects
	_render_object* obj = new _render_object();
	obj->object = (void*)lines;
	obj->type = _render_object_type::line;
	obj->position = line[0].position;
	// add the object to the vector
	objects.push_back(obj);
}


// highlight an edge, red and thicker
void _render::highlightEdge(int x_from, int y_from, int x_to, int y_to) {
	// create a line
	sfVertexArray* lines = sfVertexArray_create();	
	// set the primitive type
	sfVertexArray_setPrimitiveType(lines, sfQuads);

	sfVector2f v1 = {x_from, y_from};
	sfVector2f v2 = {x_to, y_to};	

	addThickEdge(lines, v1, v2, 3.0f, sfRed);

	_render_object* obj = new _render_object();
	obj->object = (void*)lines;
	obj->type = _render_object_type::line;
	obj->position = { v1.x,v1.y };
	// add the object to the vector
	objects.push_back(obj);
}


void _render::addThickEdge(sfVertexArray* vertexArray, sfVector2f start, sfVector2f end, float thickness, sfColor color) {
	sfVector2f direction = { end.x - start.x, end.y - start.y };
	float length = sqrt(direction.x * direction.x + direction.y * direction.y);

	sfVector2f unit = { direction.x / length, direction.y / length };
	sfVector2f perpendicular = { -unit.y, unit.x };

	sfVector2f offset = { perpendicular.x * thickness / 2, perpendicular.y * thickness / 2 };

	sfVertex v1 = { {start.x - offset.x, start.y - offset.y}, color };
	sfVertex v2 = { {start.x + offset.x, start.y + offset.y}, color };
	sfVertex v3 = { {end.x + offset.x, end.y + offset.y}, color };
	sfVertex v4 = { {end.x - offset.x, end.y - offset.y}, color };

	sfVertexArray_append(vertexArray, v1);
	sfVertexArray_append(vertexArray, v2);
	sfVertexArray_append(vertexArray, v3);
	sfVertexArray_append(vertexArray, v4);
}
	