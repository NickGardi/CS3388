/*            PURPOSE : Building the synthetic camera for 3D viewing

 PREREQUISITES : matrix.h
	author : nicholas gardi 250868721
 */
#include "GL/freeglut.h"
#include <math.h>
#include<iostream>
#include "matrix.h"

#define M_PI        3.141592653589793238462643383279502884

//change these to change the size of the shapes
#define radius 10		//sphere	
#define thickness 5;	//torus
#define width 20;		//torus

#define Ex 60.0			//how far away camera is
#define Ey 60.0
#define Ez 60.0

#define Gx 45.0
#define Gy 0.0
#define Gz 45.0

#define UPx 90.0
#define UPy 90.0
#define UPz 90.0

#define NP 5.0
#define FP 50.0

#define THETA 90.0

#define W  512
#define H  512

#define ASPECT (double)W/(double)H 

const int windowW = 512;
const int windowH = 512;

void OnDisplay();
void OnKeyboard(unsigned char key, int x, int y);
void Draw();

unsigned char frame[windowW * windowH * 3];

dmatrix_t *build_camera_matrix(dmatrix_t *E, dmatrix_t *G) {

	dmatrix_t N; /* Viewing axis */

	N = *dmat_normalize(dmat_sub(E, G));
	N.l = 3;

	dmatrix_t UP;
	dmat_alloc(&UP, 4, 1);
	UP.l = 3;

	UP.m[1][1] = UPx;
	UP.m[2][1] = UPy;
	UP.m[3][1] = UPz;
	UP.m[4][1] = 1.0;

	dmatrix_t U;

	U = *dmat_normalize(dcross_product(&UP, &N));

	dmatrix_t V;
	V = *dcross_product(&N, &U);

	dmatrix_t Mv; /* Build matrix M_v */
	dmat_alloc(&Mv, 4, 4);

	Mv.m[1][1] = U.m[1][1];
	Mv.m[1][2] = U.m[2][1];
	Mv.m[1][3] = U.m[3][1];
	Mv.m[1][4] = -1.0*((*E).m[1][1] * U.m[1][1] + (*E).m[2][1] * U.m[2][1] + (*E).m[3][1] * U.m[3][1]);

	Mv.m[2][1] = V.m[1][1];
	Mv.m[2][2] = V.m[2][1];
	Mv.m[2][3] = V.m[3][1];
	Mv.m[2][4] = -1.0*((*E).m[1][1] * V.m[1][1] + (*E).m[2][1] * V.m[2][1] + (*E).m[3][1] * V.m[3][1]);

	Mv.m[3][1] = N.m[1][1];
	Mv.m[3][2] = N.m[2][1];
	Mv.m[3][3] = N.m[3][1];
	Mv.m[3][4] = -1.0*((*E).m[1][1] * N.m[1][1] + (*E).m[2][1] * N.m[2][1] + (*E).m[3][1] * N.m[3][1]);

	Mv.m[4][1] = 0.0;
	Mv.m[4][2] = 0.0;
	Mv.m[4][3] = 0.0;
	Mv.m[4][4] = 1.0;

	dmatrix_t Mp; /* Build matrix Mp */
	dmat_alloc(&Mp, 4, 4);
	Mp = *dmat_identity(&Mp);

	float a = -1.0*(FP + NP) / (FP - NP);
	float b = -2.0*(FP*NP) / (FP - NP);

	Mp.m[1][1] = NP;
	Mp.m[2][2] = NP;
	Mp.m[3][3] = a;
	Mp.m[3][4] = b;
	Mp.m[4][3] = -1.0;
	Mp.m[4][4] = 0.0;

	/* Build matrices T_1 and S_1 */

	/* Work out coordinates of near plane corners */

	float top = NP * tan(M_PI / 180.0*THETA / 2.0);
	float right = ASPECT * top;
	float bottom = -top;
	float left = -right;

	dmatrix_t T1;
	dmat_alloc(&T1, 4, 4);

	T1 = *dmat_identity(&T1);
	T1.m[1][4] = -(right + left) / 2.0;
	T1.m[2][4] = -(top + bottom) / 2.0;

	dmatrix_t S1;
	dmat_alloc(&S1, 4, 4);

	S1 = *dmat_identity(&S1);
	S1.m[1][1] = 2.0 / (right - left);
	S1.m[2][2] = 2.0 / (top - bottom);

	/* Build matrices T2, S2, and W2 */

	dmatrix_t T2;
	dmatrix_t S2;
	dmatrix_t W2;

	dmat_alloc(&T2, 4, 4);
	dmat_alloc(&S2, 4, 4);
	dmat_alloc(&W2, 4, 4);

	T2 = *dmat_identity(&T2);
	S2 = *dmat_identity(&S2);
	W2 = *dmat_identity(&W2);

	T2.m[1][4] = 1.0;
	T2.m[2][4] = 1.0;

	S2.m[1][1] = W / 2.0;
	S2.m[2][2] = H / 2.0;

	W2.m[2][2] = -1.0;
	W2.m[2][4] = (double)H;

	return dmat_mult(&W2, dmat_mult(&S2, dmat_mult(&T2, dmat_mult(&S1, dmat_mult(&T1, dmat_mult(&Mp, &Mv))))));
}

dmatrix_t *perspective_projection(dmatrix_t *P) {

	(*P).m[1][1] /= (*P).m[4][1];
	(*P).m[2][1] /= (*P).m[4][1];
	(*P).m[3][1] /= (*P).m[4][1];
	(*P).m[4][1] /= (*P).m[4][1];

	return P;
}

int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(windowW, windowH);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("Assignment 1");

	glClearColor(1.0, 1.0, 1.0, 1.0);
	glShadeModel(GL_FLAT);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	glutDisplayFunc(OnDisplay);
	glutKeyboardFunc(OnKeyboard);

	//-- run the program
	glutMainLoop();
	return 0;
}

void OnKeyboard(unsigned char key, int x, int y) {
	switch (key) {
	case 'q':
		exit(0);
		break;
	}
}

void DrawPixel(unsigned int x, unsigned int y, unsigned char r, unsigned char g, unsigned char b) {
	if (x >= windowW || y >= windowH)
		return;

	unsigned int index = 3 * (y * windowW + x);
	frame[index] = r;
	frame[index + 1] = g;
	frame[index + 2] = b;
}

void OnDisplay() {
	memset(frame, 255, windowW * windowH * 3);
	Draw();

	glClear(GL_COLOR_BUFFER_BIT);
	glDrawPixels(windowW, windowH, GL_RGB, GL_UNSIGNED_BYTE, (GLubyte*)frame);
	glutSwapBuffers();
	glFlush();
}

//from assignment 1
void Bresenham(int x1, int y1, int x2, int y2)
{
	int w = x2 - x1;
	int h = y2 - y1;
	int dx1 = 0, dy1 = 0, dx2 = 0, dy2 = 0;
	if (w < 0) dx1 = -1; else if (w > 0) dx1 = 1;
	if (h < 0) dy1 = -1; else if (h > 0) dy1 = 1;
	if (w < 0) dx2 = -1; else if (w > 0) dx2 = 1;
	int longest = abs(w);
	int shortest = abs(h);
	if ((longest <= shortest)) {
		longest = abs(h);
		shortest = abs(w);
		if (h < 0) dy2 = -1; else if (h > 0) dy2 = 1;
		dx2 = 0;
	}
	int numerator = longest >> 1;
	for (int i = 0; i <= longest; i++) {
		DrawPixel(x1, y1, 0, 0, 0);
		numerator += shortest;
		if (!(numerator < longest)) {
			numerator -= longest;
			x1 += dx1;
			y1 += dy1;
		}
		else {
			x1 += dx2;
			y1 += dy2;
		}
	}
}

/*
Sphere Matrix Points Calculator. Creates the matrix points for a sphere with the given parameters
parameters u and v are looped through in sphereWiremeshRender to draw all the lines in the wiremesh
*/
dmatrix_t *sphereMatrix(dmatrix_t *m, double u, double v) {
	//parametric equations for a sphere to create the P matrix
	(*m).m[1][1] = radius * cos(v)*sin(u);
	(*m).m[2][1] = radius * sin(v)*sin(u);
	(*m).m[3][1] = radius * cos(u);
	(*m).m[4][1] = 1;

	return m;
}

/*
Torus Matrix Points Calculator. Creates the matrix points for the torus with the given parameters
parameters a and c are changed to change the size of the shape
parameters u and v are looped through in torusWiremeshRender to draw all the lines in the wiremesh 
*/
dmatrix_t *torusMatrix(dmatrix_t *m, double a, double c, double u, double v) {
	//parametric equations for a sphere to create the P matrix
	(*m).m[1][1] = (c + a * cos(v))*cos(u);
	(*m).m[2][1] = (c + a * cos(v))*sin(u);
	(*m).m[3][1] = a * sin(v);
	(*m).m[4][1] = 1;

	return m;
}

/*
 draws the lines of the sphere wiremesh with bresenham() after calculating and projecting the points
 parameter C is the camera matrix to be used when persepective projecting 
 */
void sphereWiremeshRender(dmatrix_t *C) {

	int x;
	double u, v;
	dmatrix_t P[3];

	//dmat_alloc new matrix P
	dmat_alloc(&P[0], 4, 1);
	dmat_alloc(&P[1], 4, 1);
	dmat_alloc(&P[2], 4, 1);

	//nested for loop to get all the points for a sphere
	for (u = 0; u < 2 * M_PI; u += (2 * M_PI / 36)) {
		for (v = 0; v < 2 * M_PI; v += (M_PI / 18)) {
			//persepctive project the points
			P[0] = *perspective_projection(dmat_mult(C, sphereMatrix(&P[0], u, v)));
			P[1] = *perspective_projection(dmat_mult(C, sphereMatrix(&P[1], u + (2 * M_PI / 36), v)));
			P[2] = *perspective_projection(dmat_mult(C, sphereMatrix(&P[2], u, v + (M_PI / 18))));
			//draw lines with bresenham algorithm
			Bresenham((int)P[0].m[1][1], (int)P[0].m[2][1], (int)P[1].m[1][1], (int)P[1].m[2][1]);
			Bresenham((int)P[0].m[1][1], (int)P[0].m[2][1], (int)P[2].m[1][1], (int)P[2].m[2][1]);
		}
	}
}

/*
 draws the lines of the torus wiremesh with bresenham() after calculating and projecting the points
  parameter C is the camera matrix to be used when persepective projecting
 */
void torusWiremeshRender(dmatrix_t *C) {

	int x;
	double u, v;
	dmatrix_t P[3];
	double a = thickness;
	double c = width;

	//dmat_alloc new matrix P
	dmat_alloc(&P[0], 4, 1);
	dmat_alloc(&P[1], 4, 1);
	dmat_alloc(&P[2], 4, 1);

	//nested for loop to get all the points for a sphere
	for (u = 0; u < 2*M_PI; u += (2 * M_PI / 36)) {
		for (v = 0; v < 2*M_PI; v += (M_PI / 18)) {
			//persepctive project the points
			P[0] = *perspective_projection(dmat_mult(C, torusMatrix(&P[0], a, c, u, v)));
			P[1] = *perspective_projection(dmat_mult(C, torusMatrix(&P[1], a, c, u + (2 * M_PI / 36), v)));
			P[2] = *perspective_projection(dmat_mult(C, torusMatrix(&P[2], a, c, u, v + (M_PI / 18))));
			//draw lines with bresenham algorithm
			Bresenham((int)P[0].m[1][1], (int)P[0].m[2][1], (int)P[1].m[1][1], (int)P[1].m[2][1]);
			Bresenham((int)P[0].m[1][1], (int)P[0].m[2][1], (int)P[2].m[1][1], (int)P[2].m[2][1]);
		}
	}
}

/*
code from main() from camera.c given
*/
void Draw() {

	dmatrix_t E; /* The centre of projection for the camera */

	dmat_alloc(&E, 4, 1);

	E.m[1][1] = Ex;
	E.m[2][1] = Ey;
	E.m[3][1] = Ez;
	E.m[4][1] = 1.0;

	dmatrix_t G; /* Point gazed at by camera */

	dmat_alloc(&G, 4, 1);

	G.m[1][1] = Gx;
	G.m[2][1] = Gy;
	G.m[3][1] = Gz;
	G.m[4][1] = 1.0;

	dmatrix_t C; /* The camera matrix */

	dmat_alloc(&C, 4, 4);
	C = *build_camera_matrix(&E, &G);

	printf("Camera Matrix:\n");
	write_dmatrix(&C);

	dmatrix_t P;

	dmat_alloc(&P, 4, 1);

	P.m[1][1] = 0.0;
	P.m[2][1] = 0.0;
	P.m[3][1] = 0.0;
	P.m[4][1] = 1.0;

	printf("Point (0,0,0) multiplied with camera matrix:\n");
	write_dmatrix(dmat_mult(&C, &P));

	printf("Point (0,0,0) after prespective projection:\n");
	write_dmatrix(perspective_projection(dmat_mult(&C, &P)));

	//calls the render function for the wiremesh of both shapes
	sphereWiremeshRender(&C);
	torusWiremeshRender(&C);
	
}












