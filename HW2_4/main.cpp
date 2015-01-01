#include <gl/freeglut.h>
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>

class A {
public:
	void DrawFractalStar()
	{
		for (int i = 0; i < 5; ++i)
		{
			DrawStarRecursively();
			glRotatef(72.f, 0, 0, 1);
		}
	}
private:
	void DrawStarRecursively(int i = 0)
	{
		if (i > 10)
		{
			return;
		}
		DrawStar();
		glPushMatrix();
		glTranslatef(0, 1, 0);
		float pentagramRatio = 2 * COS(72);
		glScalef(pentagramRatio, pentagramRatio, pentagramRatio);
		DrawStarRecursively(++i);
		glPopMatrix();
	}
	void DrawStar()
	{
		glBegin(GL_LINE_LOOP);
		float x = 54, y = 54;
		glVertex2f(COS(x), SIN(y));
		glVertex2f(COS(x += 72), SIN(y += 72));
		glVertex2f(COS(x += 72), SIN(y += 72));
		glVertex2f(COS(x += 72), SIN(y += 72));
		glVertex2f(COS(x += 72), SIN(y += 72));
		glVertex2f(COS(x += 72), SIN(y += 72));
		glVertex2f(COS(x += 144), SIN(y += 144));
		glVertex2f(COS(x += 144), SIN(y += 144));
		glVertex2f(COS(x += 144), SIN(y += 144));
		glVertex2f(COS(x += 144), SIN(y += 144));
		glEnd();
	}
	float ANGLE(float angle) {
		return static_cast<float>(angle * M_PI / 180.f);
	}
	float COS(float angle) {
		return cos(ANGLE(angle));
	}
	float SIN(float angle) {
		return sin(ANGLE(angle));
	}
};

void display() {
	glLoadIdentity();
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(0.0f, 0.0f, 0.0f);
	A a;
	glScalef(0.35f, 0.35f, 0.35f);
	a.DrawFractalStar();
	glFinish();
}

void reshape(int width, int height)
{
	if (width > height)
		glViewport((width - height) / 2, 0, height, height);
	else
		glViewport(0, (height - width) / 2, width, width);
}

int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutCreateWindow("OpenGL");
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutMainLoop();
	return 0;
}