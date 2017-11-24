
CPSC 589 Assignment 2 (W2017)
Geordie Tait 10013837

This program uses OpenGL and GLFW to draw B-splines, with
an option to generate surfaces of revolution.

Keyboard commands:
-----------------

Up/Down			Increase/decrease resolution of u parameter
				(Decrease/increase u parameter increment)
				Affects B-spline curves

Right/Left		Increase/decrease resolution of v parameter 
				(Decrease/increase v parameter increment)
				Affects surface of revolution curves

Page up/down	Increase/decrease the B-spline order by 1
				(default: 2)

Spacebar		Toggle surface of revolution generation
				(default: off)

Home			Toggle display of control points
				(default: on)

End				Toggle line smoothing/antialiasing
				(disabled during surface of revolution)
				(default: on)

Delete			Toggle colours
				(default: on)

Enter			Reset view (rotation and scaling)

Backspace		Reset program to initial parameters


Mouse commands:
--------------

Left click		Place a new control point or select existing one
				Drag while clicking an existing point to move it

Right click		Delete a control point

Middle click	Drag to rotate in 3D

Mousewheel		Adjust scale/zoom


Compilation:
-----------

The program was compiled on Windows 10 using Visual Studio Community 2015.
To compile on Linux, use the accompanying Makefile.