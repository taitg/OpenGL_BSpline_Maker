/*
B-Spline Maker
Compiled and tested using Visual Studio Community 2015 on Windows 10
*/

#ifdef _WIN32
#include <Windows.h>
#endif

#include <glm/glm.hpp>
#include <GLFW/glfw3.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <iostream>
#include <vector>
#include <math.h>

using namespace std;
using namespace glm;

#define PI		3.14159265359
#define TWOPI	PI*2

// Global vars
GLFWwindow *window;
int w, h;
double mouseX, mouseY;

vector<double> knots;		// knot sequence vector
vector<dvec2> control;		// 2d control points vector
float cRadius = 0.01f;		// control point radius
int selected = -1;			// index of selected control point
bool movePoint = false;		// flag for if a control point is moving
bool showPoints = true;		// flag for showing control points or not
bool showColours = true;	// flag for showing colours or not
bool niceLines = true;		// flag for smoothing/antialiasing lines
bool sRevolution = false;	// flag for surface of revolution generation

int k = 2;				// B-spline order
double uinc = 0.001;	// u parameter increment (for B-spline curve)
double vinc = PI/16;	// v parameter increment (for surface of revolution)

bool rotating = false;		// flag for if view is rotating or not
float yangle = 0.f;			// rotation angle around y axis
float zangle = 0.f;			// rotation angle around z axis
float rotSpeed = 100.f;		// speed factor for rotation
float scale = 1.f;			// scale factor
float zoomSpeed = 0.02f;	// speed factor for scaling

/*
 * Function to generate delta value for a given u
 * Based on algorithm given in lecture
 */
int delta(double u) {
	int m = control.size();
	for (int i = 0; i <= m + k - 1; i++) {
		if (u >= knots[i] && u < knots[i + 1])
			return i;
	}
	return -1;
}

/*
 * Function to generate a point on the B-spline curve from a given u
 * Based on efficient algorithm given in lecture
 */
dvec2 bspline(double u, int d) {

	dvec2 *c = new dvec2[control.size()];
	for (int i = 0; i <= k - 1; i++) {
		c[i] = control[d - i];
	}

	for (int r = k; r >= 2; r--) {
		int i = d;
		for (int s = 0; s <= r - 2; s++) {
			double u_i = knots[i];
			double u_ir1 = knots[i + r - 1];
			double omega = (u - u_i) / (u_ir1 - u_i);
			c[s] = omega * c[s] + (1 - omega) * c[s + 1];
			i--;
		}
	}

	dvec2 result = c[0];
	delete[] c;
	return result;
}

/*
 * Function to generate a standard knot sequence based
 * on the current order and number of control points
 */
void generateKnots() {
	knots.clear();

	for (int i = 0; i < k; i++)
		knots.push_back(0);

	int middle = control.size() - k;
	for (int i = 0; i < middle; i++)
		knots.push_back(double(i+1) / (middle+1));

	for (int i = 0; i < k; i++) 
		knots.push_back(1);
}

/*
 * Rendering function
 */
void render() {
	glEnable(GL_DEPTH_TEST);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// set line quality
	if (niceLines && !sRevolution) {
		glEnable(GL_LINE_SMOOTH);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	}
	else {
		glDisable(GL_LINE_SMOOTH);
		glDisable(GL_BLEND);
	}

	// Functions for changing transformation matrix
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glScalef(scale, scale, scale);			// scale by factor
	glRotatef(yangle, 0.0f, 1.0f, 0.0f);	// rotate by y angle
	glRotatef(zangle, 0.0f, 0.0f, 1.0f);	// rotate by z angle

	// Functions for changing projection matrix
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-1, 1, -1, 1, -10, 10);

	// draw control points if shown
	if (showPoints) {
		glBegin(GL_QUADS);

		for (int i = 0; i < control.size(); i++) {
			float pr = 1.f;
			float pg = 1.f;
			float pb = 1.f;
			if (selected == i) {
				pg = 0.f;
				pb = 0.f;
			}

			glColor3f(pr, pg, pb);
			glVertex3f(control[i].x + cRadius, control[i].y + cRadius, 0);
			glColor3f(pr, pg, pb);
			glVertex3f(control[i].x + cRadius, control[i].y - cRadius, 0);
			glColor3f(pr, pg, pb);
			glVertex3f(control[i].x - cRadius, control[i].y - cRadius, 0);
			glColor3f(pr, pg, pb);
			glVertex3f(control[i].x - cRadius, control[i].y + cRadius, 0);
		}

		glEnd();
	}

	// draw curve
	glLineWidth(3.5);
	glBegin(GL_LINE_STRIP);

	// populate (standard) knot sequence
	generateKnots();

	// generate B-spline
	for (double u = knots[k-1] + uinc; u <= knots[control.size()]; u += uinc) {
		// get delta value for this u
		int d = delta(u);

		// generate vertex/vertices and colour for this u
		if (control.size() >= d) {
			float cr = 1;
			float cg = 1;
			float cb = 1;
			if (showColours) {
				cr = 0.5 * (sin(101 * u) + 1);
				cg = 0.5 + 0.25 * (cos(11 * u) + 1);
				cb = 0.5 * (sin(71 * u) + 1);
			}
			glColor3f(cr, cg, cb);

			dvec2 point = bspline(u, d);
			if (!sRevolution) glVertex3f(point.x, point.y, 0);
			else {
				// generate surface of revolution if shown
				for (float v = 0.f; v < TWOPI; v += vinc) {
					float xr = point.x * cos(v);
					float zr = point.x * sin(v);
					glVertex3f(xr, point.y, zr);
				}
			}
		}
	}

	glEnd();
}

/*
 * Keyboard input callback
 */
void keyboard (GLFWwindow *sender, int key, int scancode, int action, int mods) {
	// close window
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		glfwSetWindowShouldClose(window, GL_TRUE);

	// increase order
	if (key == GLFW_KEY_PAGE_UP && (action == GLFW_PRESS || action == GLFW_REPEAT))
		k++;

	// decrease order
	if (key == GLFW_KEY_PAGE_DOWN && (action == GLFW_PRESS || action == GLFW_REPEAT))
		if (k > 1) k--;

	// increase u parameter resolution
	if (key == GLFW_KEY_UP && (action == GLFW_PRESS || action == GLFW_REPEAT))
		uinc *= 0.9;

	// decrease u parameter resultion
	if (key == GLFW_KEY_DOWN && (action == GLFW_PRESS || action == GLFW_REPEAT))
		uinc *= 1.1;

	// increase v parameter resolution
	if (key == GLFW_KEY_RIGHT && (action == GLFW_PRESS || action == GLFW_REPEAT))
		if (vinc > PI/32) vinc /= 2;

	// decrease v parameter resolution
	if (key == GLFW_KEY_LEFT && (action == GLFW_PRESS || action == GLFW_REPEAT))
		if (vinc < PI) vinc *= 2;

	// toggle surface of revolution
	if (key == GLFW_KEY_SPACE && action == GLFW_PRESS)
		sRevolution = !sRevolution;

	// toggle control points
	if (key == GLFW_KEY_HOME && action == GLFW_PRESS)
		showPoints = !showPoints;

	// toggle nice lines
	if (key == GLFW_KEY_END && action == GLFW_PRESS)
		niceLines = !niceLines;

	// toggle colours
	if (key == GLFW_KEY_DELETE && action == GLFW_PRESS)
		showColours = !showColours;

	// reset view
	if (key == GLFW_KEY_ENTER && action == GLFW_PRESS) {
		scale = 1.f;
		yangle = 0.f;
		zangle = 0.f;
	}

	// reset program
	if (key == GLFW_KEY_BACKSPACE && action == GLFW_PRESS) {
		control.clear();
		k = 2;
		uinc = 0.001;
		showPoints = true;
		sRevolution = false;
		scale = 1.f;
		yangle = 0.f;
		zangle = 0.f;
	}
}

/*
 * Mouse click callback
 */
void mouseClick (GLFWwindow *sender, int button, int action, int mods) {

    if (action == GLFW_PRESS) {
		selected = -1;
		double x = (2 * mouseX / w) - 1;
		double y = (-2 * mouseY / h) + 1;

		// select a control point and begin moving it while mouse button is held
		for (int i = 0; i < control.size(); i++) {
			if (abs(control[i].x - x) <= cRadius && abs(control[i].y - y) <= cRadius) {
				selected = i;
				cout << "Selected point: " << selected << endl;
				movePoint = true;
			}
		}

		// left click action (new)
		if (button == GLFW_MOUSE_BUTTON_LEFT && selected == -1) {
			// make a new control point
			control.push_back(vec2(x, y));
		}

		// right click action (delete)
		else if (button == GLFW_MOUSE_BUTTON_RIGHT && selected >= 0) {
			// delete selected control point
			control.erase(control.begin() + selected);
			cout << "Deleted point: " << selected << endl;
			selected = -1;
		}

		// middle click action (rotate view)
		else if (button == GLFW_MOUSE_BUTTON_MIDDLE) {
			glfwGetCursorPos(window, &mouseX, &mouseY);
			rotating = true;
		}
    }

	// end moving and rotating when mouse button(s) released
	if (action == GLFW_RELEASE) {
		movePoint = false;
		rotating = false;
	}
}

/*
 * Mouse position callback
 */
void mousePos(GLFWwindow *sender, double x, double y) {
	// rotate view
	if (rotating) {
		yangle += rotSpeed * float(x - mouseX) / float(w);
		zangle += rotSpeed * float(y - mouseY) / float(h);
	}

	mouseX = x;
	mouseY = y;

	// move a selected point
	if (movePoint && selected >= 0) {
		double newx = (2 * mouseX / w) - 1;
		double newy = (-2 * mouseY / h) + 1;
		control[selected] = vec2(newx, newy);
	}
}

/*
Mouse scroll event function
*/
void mouseScroll(GLFWwindow *sender, double x, double y) {
	// increase or decrease scale factor
	scale += zoomSpeed * y;
}

/*
 * Main program function
 */
int main () {
	// initialize GLFW
	if (!glfwInit()) return 1;

	// create window
	window = glfwCreateWindow (768, 768, "B-Spline Generator", NULL, NULL);
	if (!window) return 1;
	glfwMakeContextCurrent (window);

	// set input callback functions
	glfwSetKeyCallback (window, keyboard);
	glfwSetMouseButtonCallback (window, mouseClick);
	glfwSetCursorPosCallback (window, mousePos);
	glfwSetScrollCallback(window, mouseScroll);

	// main program loop
	while (!glfwWindowShouldClose (window)) {
		glfwGetFramebufferSize (window, &w, &h);
		glViewport (0, 0, w, h);

		render ();

		glfwSwapBuffers (window);
		glfwPollEvents();
	}

	// clean up and exit
	glfwDestroyWindow (window);
	glfwTerminate();
	return 0;
}


