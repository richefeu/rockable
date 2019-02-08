#ifndef GLTOOLS_HPP_6D6BEBAF
#define GLTOOLS_HPP_6D6BEBAF

#include <cstdarg> // for va_start and va_end
#include <cstring> // for strcpy
#include <cmath>

#include "OBB.hpp"

class facetSphere
{
protected:
	
	static void normalize(GLfloat *a) {
		GLfloat d = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
		a[0] /= d;
		a[1] /= d;
		a[2] /= d;
	}

	static void drawtri(GLfloat *a, GLfloat *b, GLfloat *c, int div, float r) {
		if (div <= 0) {
			glNormal3fv(a);
			glVertex3f(-a[0]*r, -a[1]*r, -a[2]*r);
			glNormal3fv(b);
			glVertex3f(-b[0]*r, -b[1]*r, -b[2]*r);
			glNormal3fv(c);
			glVertex3f(-c[0]*r, -c[1]*r, -c[2]*r);
		}
		else {
			GLfloat ab[3], ac[3], bc[3];
			for (int i = 0; i < 3; i++) {
				ab[i] = (a[i] + b[i]) / 2;
				ac[i] = (a[i] + c[i]) / 2;
				bc[i] = (b[i] + c[i]) / 2;
			}
			normalize(ab);
			normalize(ac);
			normalize(bc);
			drawtri(a, ab, ac, div - 1, r);
			drawtri(b, bc, ab, div - 1, r);
			drawtri(c, ac, bc, div - 1, r);
			drawtri(ab, bc, ac, div - 1, r);
		}
	}
	
public:
	
	static GLfloat vdata[12][3];
	static GLuint tindices[20][3];
	
	static void draw(int ndiv, float radius)
	{
		glBegin (GL_TRIANGLES);
		for (int i = 0 ; i < 20 ; ++i)
			drawtri(vdata[tindices[i][0]], vdata[tindices[i][1]], vdata[tindices[i][2]], ndiv, radius);
		glEnd();
	}
};

#define X .525731112119133606
#define Z .850650808352039932
GLfloat facetSphere::vdata[12][3] = {
	{ -X, 0.0, Z}, {X, 0.0, Z}, { -X, 0.0, -Z}, {X, 0.0, -Z},
	{0.0, Z, X}, {0.0, Z, -X}, {0.0, -Z, X}, {0.0, -Z, -X},
	{Z, X, 0.0}, { -Z, X, 0.0}, {Z, -X, 0.0}, { -Z, -X, 0.0}
};
#undef X
#undef Z

GLuint facetSphere::tindices[20][3] = {
	{0, 4, 1}, {0, 9, 4}, {9, 5, 4}, {4, 5, 8}, {4, 8, 1},
	{8, 10, 1}, {8, 3, 10}, {5, 3, 8}, {5, 2, 3}, {2, 7, 3},
	{7, 10, 3}, {7, 6, 10}, {7, 11, 6}, {11, 0, 6}, {0, 1, 6},
	{6, 1, 10}, {9, 0, 11}, {9, 11, 2}, {9, 2, 5}, {7, 2, 11}
};

class glutShape {
public:
	static void sphere(float radius, int ndiv) {
		facetSphere::draw(ndiv, radius);
	}
	
	static void drawArrow(const vec3r & orig, const vec3r & arrow, double arrowSize = -1.0, double arrowAngle = 0.7) {
		vec3r dest = orig + arrow;

		glLineWidth(2.0f);
		glBegin (GL_LINES);
		glVertex3f (orig.x, orig.y, orig.z);
		glVertex3f (dest.x, dest.y, dest.z);
		glEnd ();

		vec3r v = arrow;
		double len = v.normalize();
		if (arrowSize <= 0.0) arrowSize = 0.04 * len;
		vec3r vmz(v.x, v.y, v.z - 1.0); // v - z

		vec3r a;
		if (norm2(vmz) > 0.1) a.set(v.y, -v.x, 0.0);
		else a.set(-v.z, 0.0, v.x);
		a.normalize();
		vec3r b = cross(v, a);

		vec3r head = dest - arrowSize * v;
		vec3r c;
		double r = arrowSize * tan(0.5 * arrowAngle);
		glBegin (GL_TRIANGLE_FAN);
		glVertex3f (dest.x, dest.y, dest.z);
		for (double angle = 0.0 ; angle <= 2.0 * M_PI ; angle += 0.2 * M_PI) {
			c = cos(angle) * a + sin(angle) * b;
			glNormal3f (c.x, c.y, c.z); // Pas tout à fait juste (!)
			c = head + r * c;
			glVertex3f (c.x, c.y, c.z);
		}
		glEnd();
	}
	
	static void drawTube(vec3r & orig, vec3r & arrow, double diam) {
		vec3r dest = orig + arrow;
		vec3r v = arrow;
		v.normalize();
		vec3r vmz(v.x, v.y, v.z - 1.0); // v - z

		vec3r a;
		if (norm2(vmz) > 0.1) a.set(v.y, -v.x, 0.0);
		else a.set(-v.z, 0.0, v.x);
		if (a.isnull()) a.set(1.0, 1.0, 0.0); // for the rare cases v = (0, 0, -1)
		a.normalize();
		vec3r b = cross(v, a);

		vec3r c1, c2, n;
		double r = 0.5 * diam;
		glBegin (GL_TRIANGLE_STRIP);
		for (double angle = 0.0 ; angle <= 2.0 * M_PI ; angle += 0.05 * M_PI) {
			n =  cos(angle) * a + sin(angle) * b;
			glNormal3f(n.x, n.y, n.z);
			n *= r;
			c1 = orig + n;
			c2 = dest + n;
			glVertex3f (c1.x, c1.y, c1.z);
			glVertex3f (c2.x, c2.y, c2.z);
		}
		glEnd();
	}
	
	static void drawObb(OBB & obb)
	{
		glDisable(GL_LIGHTING);

		glLineWidth(1.0f);

		vec3r orig = obb.center;
		vec3r e0 = obb.extent[0] * obb.e[0];
		vec3r e1 = obb.extent[1] * obb.e[1];
		vec3r e2 = obb.extent[2] * obb.e[2];
	
		vec3r p0 = orig - e0 + e1 + e2;
		vec3r p1 = orig - e0 - e1 + e2;
		vec3r p2 = orig - e0 - e1 - e2;
		vec3r p3 = orig - e0 + e1 - e2;
	
		vec3r p4 = orig + e0 + e1 + e2;
		vec3r p5 = orig + e0 - e1 + e2;
		vec3r p6 = orig + e0 - e1 - e2;
		vec3r p7 = orig + e0 + e1 - e2;

		glBegin (GL_LINE_LOOP);
		glVertex3f (p0.x, p0.y, p0.z);
		glVertex3f (p1.x, p1.y, p1.z);
		glVertex3f (p2.x, p2.y, p2.z);
		glVertex3f (p3.x, p3.y, p3.z);
		glEnd();
	
		glBegin(GL_LINE_LOOP);
		glVertex3f (p4.x, p4.y, p4.z);
		glVertex3f (p5.x, p5.y, p5.z);
		glVertex3f (p6.x, p6.y, p6.z);
		glVertex3f (p7.x, p7.y, p7.z);
		glEnd();
	
		glBegin (GL_LINES);
		glVertex3f (p0.x, p0.y, p0.z);
		glVertex3f (p4.x, p4.y, p4.z);

		glVertex3f (p1.x, p1.y, p1.z);
		glVertex3f (p5.x, p5.y, p5.z);

		glVertex3f (p2.x, p2.y, p2.z);
		glVertex3f (p6.x, p6.y, p6.z);

		glVertex3f (p3.x, p3.y, p3.z);
		glVertex3f (p7.x, p7.y, p7.z);
		glEnd ();
	}
	
};

// Usage:
// switch2D::go(w, h);
// Draw things with the coordinate (0, 0) sets at the bottom left of the window
// switch2D::back(); // All attributes abd matrix are set back
class switch2D
{
public:
	static void go(int w, int h) {
		// (x,y) is from the bottom left of the window
	    glMatrixMode(GL_PROJECTION);
	    glPushMatrix();
	    glLoadIdentity();
	    glOrtho(0, w, 0, h, -1.0f, 1.0f);
	    glMatrixMode(GL_MODELVIEW);
	    glPushMatrix();
	    glLoadIdentity();
	    glPushAttrib(GL_DEPTH_TEST);
	    glDisable(GL_DEPTH_TEST);
		glPushAttrib(GL_LIGHTING);
		glDisable(GL_LIGHTING);
	}
	
	static void back() {
	    glPopAttrib();
	    glMatrixMode(GL_PROJECTION);
	    glPopMatrix();
	    glMatrixMode(GL_MODELVIEW);
	    glPopMatrix();
	}
};


class glText
{
public:
	
	// usage example: glText::print(GLUT_BITMAP_8_BY_13, 15, 100, "the value is %d", 12);
	static void print(void * font, int x, int y, const char * fmt, ...) {
		va_list args;
		va_start(args, fmt);
		char buffer[128];
		vsnprintf(buffer, 127, fmt, args);
		va_end(args);
		
	    glRasterPos2i(x, y);
		for (size_t i = 0; buffer[i]; ++i) glutBitmapCharacter(font, buffer[i]);
	}
	
	static void print3D(void * font, float x, float y, float z, const char * fmt, ...) {
		va_list args;
		va_start(args, fmt);
		char buffer[128];
		vsnprintf(buffer, 127, fmt, args);
		va_end(args);
		
	    glRasterPos3f(x, y, z);
		for (size_t i = 0; buffer[i]; ++i) glutBitmapCharacter(font, buffer[i]);
	}
};


class glTextZone
{
#define NB_LINE_MAX 40

protected:
		
	char textzone[NB_LINE_MAX][128];
	int nbLine;
	int * width;
	int * height;
	
	void drawBackground() {
		switch2D::go(*width, *height);
		
		glColor4f(1.0f, 1.0f, 1.0f, 0.6f);
		int h = nbLine * 16;
		glBegin( GL_QUADS );
		glVertex2i(0, 0);
		glVertex2i(*width, 0);
		glVertex2i(*width, h);
		glVertex2i(0, h);
		glEnd();
		
		switch2D::back();
	}
	
	void printText(int x, int y, const char *str) {
		switch2D::go(*width, *height);
		
	    glRasterPos2i(x, y);
		for (size_t i = 0; str[i]; ++i) glutBitmapCharacter(GLUT_BITMAP_9_BY_15, str[i]);
	    
		switch2D::back();
	}
	
public:
	
	// Ctors
	glTextZone(int * W, int * H): nbLine(NB_LINE_MAX), width(W), height(H) { }
	glTextZone(int n, int * W, int * H): nbLine(n), width(W), height(H) { }
	
	void increase_nbLine() {
		nbLine++;
		if (nbLine >= NB_LINE_MAX) nbLine = NB_LINE_MAX - 1;
	}
	
	void decrease_nbLine() {
		nbLine--;
		if (nbLine < 1) nbLine = 1;
	}

	void reset() {
		for (size_t i = 0 ; i < NB_LINE_MAX ; i++) textzone[i][0] = '\0';
	}

	void draw() {
		drawBackground();
		
		glColor3i(0, 0, 0);
		for (int i = 0 ; i < nbLine ; ++i) {
			printText(4, 4 + i * 16, textzone[i]);
		}
	}

	void addLine(const char * fmt, ...) {
		va_list args;
		va_start(args, fmt);
		char buffer[128];
		vsnprintf(buffer, 127, fmt, args);
		for (int i = NB_LINE_MAX - 1 ; i > 0 ; i--) {
			strcpy(textzone[i], textzone[i-1]);
		}
		sprintf ((char *) textzone[0], "%s", buffer);
		va_end(args);
	}
	
#undef NB_LINE_MAX
};


class glutTools {
public:
	
	// This funtion can be called before drawing anything.
	// It clears the screen with eventually a color gradient from bottom to top.
	// Default is light blue to white.
	static void clearBackground (bool grad, 
	int bottomRed = 135, int bottomGreen = 206, int bottomBlue = 250, 
	int topRed = 255, int topGreen = 255, int topBlue = 255)
	{
		glClearColor (1.0f, 1.0f, 1.0f, 1.0f);
		glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		if (!grad) return;

		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

		glDisable(GL_LIGHTING);
		glDisable(GL_DEPTH_TEST);

		glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
		glBegin(GL_QUADS);
		glColor3ub(bottomRed, bottomGreen, bottomBlue); // Bottom color
		glVertex2f(-1.0f, -1.0f);
		glVertex2f(1.0f, -1.0f);
		glColor3ub(topRed, topGreen, topBlue); // Top color
		glVertex2f(1.0f, 1.0f);
		glVertex2f(-1.0f, 1.0f);
		glEnd();

		glEnable(GL_LIGHTING);
		glEnable(GL_DEPTH_TEST);

		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
	}
};

#endif /* end of include guard: GLTOOLS_HPP_6D6BEBAF */
