/*            PURPOSE : Building the synthetic camera for 3D viewing

nicholas Gardi - 250868721

 PREREQUISITES : matrix.h
 */
#include "GL/freeglut.h"
#include <math.h>
#include<iostream>
#include "matrix.h"

#define M_PI        3.141592653589793238462643383279502884

 //change these to change the size of the shapes
#define sphereRadius 150		//sphere	
#define torusSmallRadius 20		//torus
#define torusLargeRadius 120	//torus
#define coneHeight 200	
#define coneRadius 100


#define Ex 200.0			//how far away camera is
#define Ey 200.0
#define Ez 200.0

#define Gx 0.0
#define Gy 0.0
#define Gz 0.0

#define UPx 0.0
#define UPy 0.0
#define UPz 1.0

#define NP 5.0
#define FP 50.0

#define THETA 90.0

#define W  512
#define H  512


#define ASPECT (double)W/(double)H 

#define X 1
#define Y 2 

const int windowW = 1000;
const int windowH = 800;
void addLighting(dmatrix_t C, dmatrix_t x, dmatrix_t y, dmatrix_t z, dmatrix_t q, dmatrix_t lightMatrix, int red, int green, int blue);
void OnDisplay();
void OnKeyboard(unsigned char key, int x, int y);
void Draw();


void XFillConvexPolygon(dmatrix_t P[], int n, int r, int g, int blue);

unsigned char frame[windowW * windowH * 3];

dmatrix_t *crossProduct(dmatrix_t vect_A, dmatrix_t vect_B) {
	dmatrix_t result;
	dmat_alloc(&result, 4, 1);
	result.m[1][1] = vect_A.m[2][1] * vect_B.m[3][1] - vect_A.m[3][1] * vect_B.m[2][1];
	result.m[2][1] = vect_A.m[1][1] * vect_B.m[3][1] - vect_A.m[3][1] * vect_B.m[1][1];
	result.m[3][1] = vect_A.m[1][1] * vect_B.m[2][1] - vect_A.m[2][1] * vect_B.m[1][1];
	result.m[4][1] = 0;

	return &result;
}

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
	(*m).m[1][1] = sphereRadius * cos(v)*sin(u);
	(*m).m[2][1] = sphereRadius * sin(v)*sin(u);
	(*m).m[3][1] = sphereRadius * cos(u);
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

dmatrix_t *coneMatrix(dmatrix_t *m, double dh, double theta) {

	//parametric equations for a cone to create the P matrix
	(*m).m[1][1] = (coneRadius * ((coneHeight - dh) / coneHeight)*cos(theta));
	(*m).m[2][1] = (coneRadius * ((coneHeight - dh) / coneHeight)*sin(theta));
	(*m).m[3][1] = dh;
	(*m).m[4][1] = 1;

	return m;
}

/*
 draws the colored in cone after calculating and projecting the points and caluclating the light vectors
 parameter C is the camera matrix to be used when persepective projecting
 */
void coneWiremeshRender(dmatrix_t *C) {

	//cone color and values for calcualtions
	int r = 255;
	int g = 215;
	int b = 0;
	double theta = 0;
	double dTheta = 2 * M_PI / 120;
	double h;
	double dh = coneHeight / 100;

	//light
	dmatrix_t lightMatrix;
	dmat_alloc(&lightMatrix, 4, 1);
	lightMatrix.m[1][1] = 1000;
	lightMatrix.m[2][1] = 1000;
	lightMatrix.m[3][1] = 1000;
	lightMatrix.m[4][1] = 1.0;

	//points to draw cone. the buffer here has same idea as in Torus'
	dmatrix_t P[3];
	for (int i = 0; i < 3; i++)
	{
		dmat_alloc(&P[i], 4, 1);
	}

	dmatrix_t x, y, z, q;
	dmat_alloc(&x, 4, 1);
	dmat_alloc(&y, 4, 1);
	dmat_alloc(&z, 4, 1);
	dmat_alloc(&q, 4, 1);

	coneMatrix(&q, dh, theta);

	q = *perspective_projection(dmat_mult(C, &q));


	//nested loop to get all the points for a cone
	for (h = coneHeight; h + dh >= 0; h -= dh)
	{
		for (theta = 0; theta <= 2 * M_PI; theta += dTheta)
		{
			//get cone matrices
			coneMatrix(&x, h, theta);
			coneMatrix(&y, h, theta + dTheta);
			coneMatrix(&z, h - dh, theta);
			coneMatrix(&q, h - dh, theta + dTheta);

			addLighting(*C, x, y, z, q, lightMatrix, r, g, b);

		}
	}
}

/*
 draws the colored in sphere after calculating and projecting the points and caluclating the light vectors
 parameter C is the camera matrix to be used when persepective projecting
 */
void sphereWiremeshRender(dmatrix_t *C) {

	//sphere color and values for calcualtions
	int r = 0;
	int g = 255;
	int b = 255;

	double v = 0;
	double u = 0;
	double dv = 2.0 * M_PI / 72;
	double du = M_PI / 32;

	//light
	dmatrix_t lightMatrix;
	dmat_alloc(&lightMatrix, 4, 1);
	lightMatrix.m[1][1] = 1000;
	lightMatrix.m[2][1] = 1000;
	lightMatrix.m[3][1] = 1000;
	lightMatrix.m[4][1] = 1.0;

	//points to draw sphere. the buffer here has same idea as in Torus'
	dmatrix_t P[3];
	for (int i = 0; i < 3; i++)
	{
		dmat_alloc(&P[i], 4, 1);
	}


	dmatrix_t x, y, z, q;
	dmat_alloc(&x, 4, 1);
	dmat_alloc(&y, 4, 1);
	dmat_alloc(&z, 4, 1);
	dmat_alloc(&q, 4, 1);


	sphereMatrix(&q, u, v);

	q = *perspective_projection(dmat_mult(C, &q));

	//nested loop to get all the points for a sphere
	for (u = 0; u < 2 * M_PI; u += (2 * M_PI / 100))
	{
		for (v = 0; v < 2 * M_PI; v += (M_PI / 40))
		{
			//get sphere matrices
			sphereMatrix(&x, u, v);
			sphereMatrix(&y, u + du, v);
			sphereMatrix(&z, u, v + dv);
			sphereMatrix(&q, u + du, v + dv);

			addLighting(*C, x, y, z, q, lightMatrix, r, g, b);

		}
	}
}

/*
 draws the colored in torus after calculating and projecting the points and caluclating the light vectors
  parameter C is the camera matrix to be used when persepective projecting
 */
void torusWiremeshRender(dmatrix_t *C) {

	//torus color and values for calcualtions
	int r = 138;
	int g = 43;
	int b = 226;
	double v = 0, u = 0;
	double dv = 2.0 * M_PI / 256;
	double du = M_PI / 256;
	double a = torusSmallRadius;
	double c = torusLargeRadius;

	//light
	dmatrix_t lightMatrix;
	dmat_alloc(&lightMatrix, 4, 1);
	lightMatrix.m[1][1] = 0;
	lightMatrix.m[2][1] = 0;
	lightMatrix.m[3][1] = 0;
	lightMatrix.m[4][1] = 1.0;

	dmatrix_t x, y, z, q;
	dmat_alloc(&x, 4, 1);
	dmat_alloc(&y, 4, 1);
	dmat_alloc(&z, 4, 1);
	dmat_alloc(&q, 4, 1);

	torusMatrix(&q, a, c, u, v);

	q = *perspective_projection(dmat_mult(C, &q));

	//nested loop to get all the points for a torus
	for (u = 0; u < 2 * M_PI; u += du)
	{
		for (v = 0; v < M_PI; v += dv)
		{
			//get torus matrices
			torusMatrix(&x, a, c, u, v);
			torusMatrix(&y, a, c, u + du, v);
			torusMatrix(&z, a, c, u, v + dv);
			torusMatrix(&q, a, c, u + du, v + dv);


			addLighting(*C, x, y, z, q, lightMatrix, r, g, b);
		}
	}
}

//this function performs the calculations for adding light to the scene

void addLighting(dmatrix_t C, dmatrix_t x, dmatrix_t y, dmatrix_t z, dmatrix_t q, dmatrix_t lightMatrix, int red, int green, int blue) {
	//Define normal vectors for each vertex of all the objects. These normals determine the
	//orientation of the object relative to the light sources.
	//calculations for light source
	dmatrix_t P[3];
	for (int i = 0; i < 3; i++)
	{
		dmat_alloc(&P[i], 4, 1);
	}

	dmatrix_t XY, YZ, XZ, ZQ, xyzNormal, yzqNormal, xyzLightVector, yzqLightVector, xyzCentrePoint, yzqCentrePoint;
	dmat_alloc(&XY, 4, 1);
	dmat_alloc(&YZ, 4, 1);
	dmat_alloc(&XZ, 4, 1);
	dmat_alloc(&ZQ, 4, 1);
	dmat_alloc(&xyzNormal, 4, 1);
	dmat_alloc(&yzqNormal, 4, 1);
	dmat_alloc(&xyzCentrePoint, 4, 1);
	dmat_alloc(&yzqCentrePoint, 4, 1);

	XY = *dmat_sub(&x, &y);
	YZ = *dmat_sub(&z, &x);
	XZ = *dmat_sub(&z, &y);
	ZQ = *dmat_sub(&q, &z);

	//matrix cross product to get normal
	xyzNormal = *crossProduct(XZ, XY);
	xyzNormal = *dmat_normalize(&xyzNormal);
	yzqNormal = *crossProduct(ZQ, YZ);
	yzqNormal = *dmat_normalize(&yzqNormal);

	xyzCentrePoint.m[1][1] = (x.m[1][1] + y.m[1][1] + z.m[1][1]) / 3; xyzCentrePoint.m[2][1] = (x.m[2][1] + y.m[2][1] + z.m[2][1]) / 3; xyzCentrePoint.m[3][1] = (x.m[3][1] + y.m[3][1] + z.m[3][1]) / 3; xyzCentrePoint.m[4][1] = 1;

	yzqCentrePoint.m[1][1] = (z.m[1][1] + y.m[1][1] + q.m[1][1]) / 3; yzqCentrePoint.m[2][1] = (z.m[2][1] + y.m[2][1] + q.m[2][1]) / 3; yzqCentrePoint.m[3][1] = (z.m[3][1] + y.m[3][1] + q.m[3][1]) / 3; yzqCentrePoint.m[4][1] = 1;

	//prospective project the points
	x = *perspective_projection(dmat_mult(&C, &x));
	y = *perspective_projection(dmat_mult(&C, &y));
	z = *perspective_projection(dmat_mult(&C, &z));
	q = *perspective_projection(dmat_mult(&C, &q));

	dmat_alloc(&xyzLightVector, 4, 1); xyzLightVector = *dmat_sub(&lightMatrix, &xyzCentrePoint);
	dmat_alloc(&yzqLightVector, 4, 1); yzqLightVector = *dmat_sub(&lightMatrix, &yzqCentrePoint);


	//calculate the angle of the light vectors
	double xyz = ddot_product(&xyzLightVector, &xyzNormal) / ((dmat_norm(&xyzLightVector)) * (dmat_norm(&xyzNormal)));
	while (xyz < 0)
		xyz = -xyz;
	double yzq = ddot_product(&yzqLightVector, &yzqNormal) / ((dmat_norm(&yzqLightVector)) * (dmat_norm(&yzqNormal)));
	while (yzq < 0)
		yzq = -yzq;


	//useXFillPolygon to add color
	P[0] = x; P[1] = y; P[2] = z;
	XFillConvexPolygon(P, 3, red * xyz, green * xyz, blue * xyz);
	P[0] = q;
	XFillConvexPolygon(P, 3, red * yzq, green * yzq, blue * yzq);
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
	//torusWiremeshRender(&C);
	//coneWiremeshRender(&C);



}




//code from XFillPolygon.c given and modified slightly
double minimum_coordinate(int coordinate, dmatrix_t P[], int n) {

	int i;
	double min;

	min = P[0].m[coordinate][1];

	for (i = 1; i < n; i++) {
		if (P[i].m[coordinate][1] < min) {
			min = P[i].m[coordinate][1];
		}
	}
	return min;
}

double maximum_coordinate(int coordinate, dmatrix_t P[], int n) {

	int i;
	double max;

	max = P[0].m[coordinate][1];

	for (i = 1; i < n; i++) {
		if (P[i].m[coordinate][1] > max) {
			max = P[i].m[coordinate][1];
		}
	}
	return max;
}

int maximum_intersection(int intersections[], int n) {

	int i, max;

	max = intersections[0];
	for (i = 1; i < n; i++) {
		if (intersections[i] > max) {
			max = intersections[i];
		}
	}
	return max;
}

int minimum_intersection(int intersections[], int n) {

	int i, min;

	min = intersections[0];
	for (i = 1; i < n; i++) {
		if (intersections[i] < min) {
			min = intersections[i];
		}
	}
	return min;
}

void XFillConvexPolygon(dmatrix_t P[], int n, int r, int g, int blue) {

	int i, j;
	int y, y_min, y_max, min_int, max_int;
	double m, b;
	int *active, *horizontal, *intersections;

	horizontal = (int *)malloc(n * sizeof(int)); /* Allocate horizontal segment table */
	active = (int *)malloc(n * sizeof(int)); /* Allocate active segment table */
	intersections = (int *)malloc(n * sizeof(int)); /* Allocate intersection table */

	y_min = (int)minimum_coordinate(Y, P, n); /* Determine number of scan lines */
	y_max = (int)maximum_coordinate(Y, P, n);

	for (i = 0; i < n; i++) {
		horizontal[i] = (int)P[i].m[Y][1] == (int)P[(i + 1) % n].m[Y][1]; /* Find horizontal segments */
	}

	for (y = y_min; y <= y_max; y++) { /* For each scan line y */
		for (i = 0; i < n; i++) {  /* Update segment table */
			if (!horizontal[i]) {
				active[i] = (y >= (int)P[i].m[Y][1] && y <= (int)P[(i + 1) % n].m[Y][1]) || (y <= (int)P[i].m[Y][1] && y >= (int)P[(i + 1) % n].m[Y][1]);
			}
		}
		j = 0;
		for (i = 0; i < n; i++) { /* find intersection x-value. The y-value is given by the scan line */
			if (active[i] && !horizontal[i]) {
				if ((int)P[i].m[X][1] == (int)P[(i + 1) % n].m[X][1]) { /* Vertical segment */
					intersections[j++] = (int)P[i].m[X][1];
				}
				else {
					m = (double)((int)P[(i + 1) % n].m[Y][1] - (int)P[i].m[Y][1]) / (double)((int)P[(i + 1) % n].m[X][1] - (int)P[i].m[X][1]); /* Compute slope and intercept */
					b = (double)((int)P[i].m[Y][1]) - m * (double)((int)P[i].m[X][1]);
					intersections[j++] = (int)(((double)y - b) / m); /* Compute intersection */
				}
			}
		}
		min_int = minimum_intersection(intersections, j);
		max_int = maximum_intersection(intersections, j) + 1;
		for (i = min_int; i < max_int; i++) { /* Tracing from minimum to maximum intersection */
			DrawPixel(i, y, r, g, blue);
		}
	}
	free(horizontal);
	free(active);
	free(intersections);
}
