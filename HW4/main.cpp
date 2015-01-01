#include <GL/glut.h>
#include "curve.h"
#include <vector>
#include <list>
#include <algorithm>
#include <iostream>
#include <limits>
#include <array>
#define _USE_MATH_DEFINES
#include <math.h>
#include "Matrices.h"
#include "texture.h"

const static int ROUND_DIV_COUNT = 32;
const static int RES = 1024;
const static int RADIUS = 50;
const static float MIN_RADIUS = 10;
const static float MAX_RADIUS = 210;
CubicBezierCurve curve;
GLsizei width = 640, height = 800;
GLsizei subWidth = width / 2, subHeight = height / 2;
int edit_ctrlpts_idx = -1;
bool isDrawControlMesh = true;
typedef std::vector<std::array<float, 4>> ConstantsDictType;
ConstantsDictType BeizerConstants;
typedef std::vector<std::array<Vector3, ROUND_DIV_COUNT>> CirclePointsType;
typedef std::vector<Vector3> CenterPointsType;
CirclePointsType Points;
CenterPointsType CenterPoints;
typedef std::vector<float> RadiusType;
RadiusType Radius;
GLint window, xyzwin, xywin, zxwin, yzwin;
void redisplay_all();
GLuint texture;
int texture_width, texture_height;
std::string texture_filename = "grace_probe.png";
int mouse_x, mouse_y;
BYTE mouse_plane = 0;
const static Vector3 x_axis(1, 0, 0), y_axis(0, 1, 0), z_axis(0, 0, 1);
bool isDampled = false;
int damplingPointIndex = -1;
Vector3 damplingInitPoint;
Vector3 damplingVector;
float damplingSecond = 0;

int hit_index(CubicBezierCurve* curve, int x, int y, int type)
{
	if (type == 0)
	{
		float new_x = 4 * x - 2 * static_cast<float>(subWidth);
		float new_y = 4 * y - 2 * static_cast<float>(subHeight);
		for (int i = 0; i < 4; i++)
		{
			REAL tx = curve->control_pts[i][0] - new_x;
			REAL ty = curve->control_pts[i][1] - new_y;
			if ((tx * tx + ty * ty) < 300) return i;
		}
		return -1;
	}
	else if (type == 1)
	{
		float new_y = 4 * x - 2 * static_cast<float>(subWidth);
		float new_z = 4 * y - 2 * static_cast<float>(subHeight);
		for (int i = 0; i < 4; i++)
		{
			REAL tx = curve->control_pts[i][1] - new_y;
			REAL ty = curve->control_pts[i][2] - new_z;
			if ((tx * tx + ty * ty) < 300) return i;
		}
		return -1;
	}
	else if (type == 2)
	{
		float new_x = 4 * x - 2 * static_cast<float>(subWidth);
		float new_z = 2 * static_cast<float>(subHeight)-4 * y;
		for (int i = 0; i < 4; i++)
		{
			REAL tx = curve->control_pts[i][0] - new_x;
			REAL ty = curve->control_pts[i][2] - new_z;
			if ((tx * tx + ty * ty) < 300) return i;
		}
		return -1;
	}
	else
	{
		return -1;
	}
}

/**
* get a point-normal
*/
Vector3 GetNormalVector(const Vector3& vector)
{
	if (vector.y != 0 || vector.z != 0)
	{
		Vector3 temp(1, 0, 0);
		return vector.cross(temp);
	}
	else
	{
		Vector3 temp(0, 1, 0);
		return vector.cross(temp);
	}
}

void AssignCirclePoints(size_t index, const Vector3& vi, const Vector3& ri, const Vector3& ui)
{
	auto& circle = Points[index];
	for (int j = 0; j < ROUND_DIV_COUNT; ++j)
	{
		Vector3 diffR = ri * Radius[index] * cosf(2 * static_cast<float>(M_PI)* j / ROUND_DIV_COUNT);
		Vector3 diffU = ui * Radius[index] * sinf(2 * static_cast<float>(M_PI)* j / ROUND_DIV_COUNT);
		circle[j].x += diffR.x + diffU.x;
		circle[j].y += diffR.y + diffU.y;
		circle[j].z += diffR.z + diffU.z;
	}
}

void AssignPoints()
{
	for (int i = 0; i <= RES; ++i)
	{
		Vector3& center = CenterPoints[i];
		center.set(0, 0, 0);
		for (int j = 0; j < 4; ++j)
		{
			center[0] += BeizerConstants[i][j] * curve.control_pts[j][0];
			center[1] += BeizerConstants[i][j] * curve.control_pts[j][1];
			center[2] += BeizerConstants[i][j] * curve.control_pts[j][2];
		}
		auto& circle = Points[i];
		circle.fill(center);
	}
	for (int i = 0; i <= RES; ++i)
	{
		Radius[i] = 0;
		for (int j = 0; j < 4; ++j)
		{
			Radius[i] += BeizerConstants[i][j] * curve.Radius[j];
		}
	}
	for (int i = 1; i < RES; ++i)
	{
		Vector3 tangent = CenterPoints[i + 1] - CenterPoints[i];
		Vector3 vi = tangent + z_axis;
		Matrix3 householder = householder.identity() - (2 / vi.dot(vi)) * Matrix3(vi.x * vi.x, vi.x * vi.y, vi.x * vi.z, vi.y * vi.x, vi.y * vi.y, vi.y * vi.z, vi.z * vi.x, vi.z * vi.y, vi.z * vi.z);
		Vector3 ui = householder * x_axis;
		Vector3 ri = householder * y_axis;
		if (i == 1)
		{
			AssignCirclePoints(0, vi, ri, ui);
		}
		AssignCirclePoints(i, vi, ri, ui);
		if (i == RES - 1)
		{
			AssignCirclePoints(RES, vi, ri, ui);
		}
	}
}

void InitControlPoints()
{
	SET_PT3(curve.control_pts[0], 50, 100, 50);
	SET_PT3(curve.control_pts[1], 200, 300, 100);
	SET_PT3(curve.control_pts[2], 400, 300, 150);
	SET_PT3(curve.control_pts[3], 550, 100, 200);
	for (size_t i = 0; i < curve.Radius.size(); ++i)
	{
		curve.Radius[i] = RADIUS;
	}
}

void init()
{
	InitControlPoints();
	BeizerConstants = ConstantsDictType(RES + 1);
	Points = CirclePointsType(RES + 1);
	CenterPoints = CenterPointsType(RES + 1);
	Radius = RadiusType(RES + 1);
	for (int i = 0; i <= RES; ++i)
	{
		float t = static_cast<float>(i) / RES;
		float one_minus_t = 1 - t;
		BeizerConstants[i][0] = one_minus_t * one_minus_t * one_minus_t;
		BeizerConstants[i][1] = 3 * one_minus_t * one_minus_t * t;
		BeizerConstants[i][2] = 3 * one_minus_t * t * t;
		BeizerConstants[i][3] = t * t * t;
	}
	AssignPoints();
}

void SetDisplayXYZ(GLint, GLint)
{
	glViewport(0, 0, subWidth, subHeight);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-width, width, -height, height, -1000, 1000);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	Vector3 base(0, 0, 0);
	Vector3 camera(1, 1, 1);
	Vector3 up = GetNormalVector(camera - base);
	gluLookAt(camera.x, camera.y, camera.z, base.x, base.y, base.z, up.x, up.y, up.z);
}

void SetDisplayXY(GLint, GLint)
{
	glViewport(0, 0, subWidth, subHeight);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-width, width, -height, height, -1000, 1000);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	Vector3 base(0, 0, 0);
	Vector3 camera(0, 0, 1);
	Vector3 up = GetNormalVector(camera - base);
	gluLookAt(camera.x, camera.y, camera.z, base.x, base.y, base.z, up.x, up.y, up.z);
}

void SetDisplayYZ(GLint, GLint)
{
	glViewport(0, 0, subWidth, subHeight);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-width, width, -height, height, -1000, 1000);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	Vector3 base(0, 0, 0);
	Vector3 camera(1, 0, 0);
	Vector3 up = GetNormalVector(camera - base);
	gluLookAt(camera.x, camera.y, camera.z, base.x, base.y, base.z, up.x, up.y, up.z);
}

void SetDisplayZX(GLint, GLint)
{
	glViewport(0, 0, subWidth, subHeight);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-width, width, -height, height, -1000, 1000);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	Vector3 base(0, 0, 0);
	Vector3 camera(0, 1, 0);
	Vector3 up = GetNormalVector(camera - base);
	gluLookAt(camera.x, camera.y, camera.z, base.x, base.y, base.z, up.x, up.y, up.z);
}

void DisplayXYZ()
{
	glClear(GL_COLOR_BUFFER_BIT);
	/* curve */
	Vector3 base(0, 0, 0);
	Vector3 camera(1, 1, 1);
	Vector3 v_prime = (base - camera).normalize();
	for (int i = 0; i < RES; i++)
	{
		int next_i = i + 1;
		glColor3f(1.0f, 1.0f, 1.0f);
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, texture);
		glBegin(GL_QUADS);
		for (int j = 0; j < ROUND_DIV_COUNT; ++j)
		{
			int next_j = (j + 1) % ROUND_DIV_COUNT;
			{
				Vector3 n_prime = (CenterPoints[i] - Points[i][j]).normalize();
				Vector3 r = v_prime - 2 * (n_prime.dot(v_prime)) * n_prime;
				float m = std::sqrt(r.x * r.x + r.y * r.y + (r.z + 1) * (r.z + 1));
				float u = r.x / (2 * m) + 0.5f;
				float v = r.y / (2 * m) + 0.5f;
				glTexCoord2f(u, v);
				glVertex3f(Points[i][j][0], Points[i][j][1], Points[i][j][2]);
			}
			{
				Vector3 n_prime = (CenterPoints[i] - Points[i][next_j]).normalize();
				Vector3 r = v_prime - 2 * (n_prime.dot(v_prime)) * n_prime;
				float m = std::sqrt(r.x * r.x + r.y * r.y + (r.z + 1) * (r.z + 1));
				float u = r.x / (2 * m) + 0.5f;
				float v = r.y / (2 * m) + 0.5f;
				glTexCoord2f(u, v);
				glVertex3f(Points[i][next_j][0], Points[i][next_j][1], Points[i][next_j][2]);
			}
			{
				Vector3 n_prime = (CenterPoints[next_i] - Points[next_i][next_j]).normalize();
				Vector3 r = v_prime - 2 * (n_prime.dot(v_prime)) * n_prime;
				float m = std::sqrt(r.x * r.x + r.y * r.y + (r.z + 1) * (r.z + 1));
				float u = r.x / (2 * m) + 0.5f;
				float v = r.y / (2 * m) + 0.5f;
				glTexCoord2f(u, v);
				glVertex3f(Points[next_i][next_j][0], Points[next_i][next_j][1], Points[next_i][next_j][2]);
			}
			{
				Vector3 n_prime = (CenterPoints[next_i] - Points[next_i][j]).normalize();
				Vector3 r = v_prime - 2 * (n_prime.dot(v_prime)) * n_prime;
				float m = std::sqrt(r.x * r.x + r.y * r.y + (r.z + 1) * (r.z + 1));
				float u = r.x / (2 * m) + 0.5f;
				float v = r.y / (2 * m) + 0.5f;
				glTexCoord2f(u, v);
				glVertex3f(Points[next_i][j][0], Points[next_i][j][1], Points[next_i][j][2]);
			}
		}
		glEnd();
	};

	/* control mesh */
	if (isDrawControlMesh)
	{
		glColor3ub(255, 0, 0);
		glBegin(GL_LINE_STRIP);
		for (int i = 0; i < 4; i++)
		{
			REAL *pt = curve.control_pts[i];
			glVertex3f(pt[0], pt[1], pt[2]);
		}
		glEnd();
	}

	/* control pts */
	glColor3ub(0, 0, 255);
	glPointSize(10.0);
	glBegin(GL_POINTS);
	for (int i = 0; i < 4; i++)
	{
		REAL *pt = curve.control_pts[i];
		glVertex3f(pt[0], pt[1], pt[2]);
	}
	glEnd();
	glutPostRedisplay();
	glutSwapBuffers();
}

void DisplayXY()
{
	glClear(GL_COLOR_BUFFER_BIT);
	/* curve */
	glColor3ub(0, 0, 0);
	glBegin(GL_TRIANGLE_STRIP);
	for (int i = 0; i < RES; i++)
	{
		for (int j = 0; j < ROUND_DIV_COUNT; ++j)
		{
			glVertex3f(Points[i][j][0], Points[i][j][1], Points[i][j][2]);
			glVertex3f(Points[i + 1][j][0], Points[i + 1][j][1], Points[i + 1][j][2]);
		}
		int j = 0;
		glVertex3f(Points[i][j][0], Points[i][j][1], Points[i][j][2]);
		glVertex3f(Points[i + 1][j][0], Points[i + 1][j][1], Points[i + 1][j][2]);
	};
	glEnd();

	/* control mesh */
	if (isDrawControlMesh)
	{
		glColor3ub(255, 0, 0);
		glBegin(GL_LINE_STRIP);
		for (int i = 0; i < 4; i++)
		{
			REAL *pt = curve.control_pts[i];
			glVertex3f(pt[0], pt[1], pt[2]);
		}
		glEnd();
	}

	/* control pts */
	glColor3ub(0, 0, 255);
	glPointSize(10.0);
	glBegin(GL_POINTS);
	for (int i = 0; i < 4; i++)
	{
		REAL *pt = curve.control_pts[i];
		glVertex3f(pt[0], pt[1], pt[2]);
	}
	glEnd();

	if (isDampled && damplingPointIndex != -1 && edit_ctrlpts_idx == -1)
	{
		if (damplingInitPoint.distance(Vector3(curve.control_pts[damplingPointIndex][0], curve.control_pts[damplingPointIndex][1], curve.control_pts[damplingPointIndex][2])) < 10)
		{
			damplingPointIndex = -1;
		}
		else
		{
			Vector3 damplingNow = damplingInitPoint + damplingVector * exp(-damplingSecond) * cosf(2 * static_cast<float>(M_PI)* damplingSecond);
			curve.control_pts[damplingPointIndex][0] = damplingNow.x;
			curve.control_pts[damplingPointIndex][1] = damplingNow.y;
			curve.control_pts[damplingPointIndex][2] = damplingNow.z;
			damplingSecond += 0.06f;
			AssignPoints();
		}
	}

	glutSwapBuffers();
	redisplay_all();
}

void DisplayYZ()
{
	glClear(GL_COLOR_BUFFER_BIT);
	/* curve */
	glColor3ub(0, 0, 0);
	for (int i = 0; i < RES; i++)
	{
		glBegin(GL_QUAD_STRIP);
		for (int j = 0; j < ROUND_DIV_COUNT; ++j)
		{
			glVertex3f(Points[i][j][0], Points[i][j][1], Points[i][j][2]);
			glVertex3f(Points[i + 1][j][0], Points[i + 1][j][1], Points[i + 1][j][2]);
		}
		int j = 0;
		glVertex3f(Points[i][j][0], Points[i][j][1], Points[i][j][2]);
		glVertex3f(Points[i + 1][j][0], Points[i + 1][j][1], Points[i + 1][j][2]);
		glEnd();
	};

	/* control mesh */
	if (isDrawControlMesh)
	{
		glColor3ub(255, 0, 0);
		glBegin(GL_LINE_STRIP);
		for (int i = 0; i < 4; i++)
		{
			REAL *pt = curve.control_pts[i];
			glVertex3f(pt[0], pt[1], pt[2]);
		}
		glEnd();
	}

	/* control pts */
	glColor3ub(0, 0, 255);
	glPointSize(10.0);
	glBegin(GL_POINTS);
	for (int i = 0; i < 4; i++)
	{
		REAL *pt = curve.control_pts[i];
		glVertex3f(pt[0], pt[1], pt[2]);
	}
	glEnd();

	if (isDampled && damplingPointIndex != -1 && edit_ctrlpts_idx == -1)
	{
		if (damplingInitPoint.distance(Vector3(curve.control_pts[damplingPointIndex][0], curve.control_pts[damplingPointIndex][1], curve.control_pts[damplingPointIndex][2])) < 10)
		{
			damplingPointIndex = -1;
		}
		else
		{
			Vector3 damplingNow = damplingInitPoint + damplingVector * exp(-damplingSecond) * cosf(2 * static_cast<float>(M_PI)* damplingSecond);
			curve.control_pts[damplingPointIndex][0] = damplingNow.x;
			curve.control_pts[damplingPointIndex][1] = damplingNow.y;
			curve.control_pts[damplingPointIndex][2] = damplingNow.z;
			damplingSecond += 0.06f;
			AssignPoints();
		}
	}

	glutSwapBuffers();
	redisplay_all();
}

void DisplayZX()
{
	glClear(GL_COLOR_BUFFER_BIT);
	/* curve */
	glColor3ub(0, 0, 0);
	for (int i = 0; i < RES; i++)
	{
		glBegin(GL_QUAD_STRIP);
		for (int j = 0; j < ROUND_DIV_COUNT; ++j)
		{
			glVertex3f(Points[i][j][0], Points[i][j][1], Points[i][j][2]);
			glVertex3f(Points[i + 1][j][0], Points[i + 1][j][1], Points[i + 1][j][2]);
		}
		int j = 0;
		glVertex3f(Points[i][j][0], Points[i][j][1], Points[i][j][2]);
		glVertex3f(Points[i + 1][j][0], Points[i + 1][j][1], Points[i + 1][j][2]);
		glEnd();
	};

	/* control mesh */
	if (isDrawControlMesh)
	{
		glColor3ub(255, 0, 0);
		glBegin(GL_LINE_STRIP);
		for (int i = 0; i < 4; i++)
		{
			REAL *pt = curve.control_pts[i];
			glVertex3f(pt[0], pt[1], pt[2]);
		}
		glEnd();
	}

	/* control pts */
	glColor3ub(0, 0, 255);
	glPointSize(10.0);
	glBegin(GL_POINTS);
	for (int i = 0; i < 4; i++)
	{
		REAL *pt = curve.control_pts[i];
		glVertex3f(pt[0], pt[1], pt[2]);
	}
	glEnd();

	if (isDampled && damplingPointIndex != -1 && edit_ctrlpts_idx == -1)
	{
		if (damplingInitPoint.distance(Vector3(curve.control_pts[damplingPointIndex][0], curve.control_pts[damplingPointIndex][1], curve.control_pts[damplingPointIndex][2])) < 10)
		{
			damplingPointIndex = -1;
		}
		else
		{
			Vector3 damplingNow = damplingInitPoint + damplingVector * exp(-damplingSecond) * cosf(2 * static_cast<float>(M_PI)* damplingSecond);
			curve.control_pts[damplingPointIndex][0] = damplingNow.x;
			curve.control_pts[damplingPointIndex][1] = damplingNow.y;
			curve.control_pts[damplingPointIndex][2] = damplingNow.z;
			damplingSecond += 0.06f;
			AssignPoints();
		}
	}

	glutSwapBuffers();
	redisplay_all();
}

// void glutKeyboardFunc(void (*func)(unsigned char key, int x, int y));
void keyboard_callback(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'i': case 'I':
		InitControlPoints();
		AssignPoints();
		break;
	case 'c': case 'C':
		isDrawControlMesh ^= true;
		break;
	case 'b': case 'B':
	{
				  int edit_ctrlpts_idx = hit_index(&curve, mouse_x, subHeight - mouse_y, mouse_plane);
				  if (edit_ctrlpts_idx != -1)
				  {
					  curve.Radius[edit_ctrlpts_idx] = std::min(curve.Radius[edit_ctrlpts_idx] + 10, MAX_RADIUS);
					  AssignPoints();
				  }
	}
		break;
	case 's': case 'S':
	{
				  int edit_ctrlpts_idx = hit_index(&curve, mouse_x, subHeight - mouse_y, mouse_plane);
				  if (edit_ctrlpts_idx != -1)
				  {
					  curve.Radius[edit_ctrlpts_idx] = std::max(curve.Radius[edit_ctrlpts_idx] - 10, MIN_RADIUS);
					  AssignPoints();
				  }
	}
		break;
	case 'd': case'D':
		isDampled ^= true;
		break;
	case (27) : exit(0); break;
	default: break;
	}
	glutPostRedisplay();
}

void xy_mouse(int button, int action, int x, int y)
{
	if (GLUT_LEFT_BUTTON == button)
	{
		switch (action)
		{
		case GLUT_DOWN:
			edit_ctrlpts_idx = hit_index(&curve, x, subHeight - y, 0);
			if (isDampled && edit_ctrlpts_idx != -1)
			{
				damplingPointIndex = edit_ctrlpts_idx;
				damplingInitPoint = Vector3(curve.control_pts[edit_ctrlpts_idx][0], curve.control_pts[edit_ctrlpts_idx][1], curve.control_pts[edit_ctrlpts_idx][2]);
			}
			break;
		case GLUT_UP:
			if (isDampled && edit_ctrlpts_idx != -1)
			{
				Vector3 damplingEndPoint = Vector3(curve.control_pts[edit_ctrlpts_idx][0], curve.control_pts[edit_ctrlpts_idx][1], curve.control_pts[edit_ctrlpts_idx][2]);
				damplingVector = damplingEndPoint - damplingInitPoint;
				damplingSecond = 0;
			}
			edit_ctrlpts_idx = -1;
			break;
		default: break;
		}
	}
	else if (button == 3) // wheel up
	{
		edit_ctrlpts_idx = hit_index(&curve, x, subHeight - y, 0);
		if (edit_ctrlpts_idx != -1)
		{
			curve.Radius[edit_ctrlpts_idx] = std::min(curve.Radius[edit_ctrlpts_idx] + 10, MAX_RADIUS);
			AssignPoints();
		}
	}
	else if (button == 4) // wheel down
	{
		edit_ctrlpts_idx = hit_index(&curve, x, subHeight - y, 0);
		if (edit_ctrlpts_idx != -1)
		{
			curve.Radius[edit_ctrlpts_idx] = std::max(curve.Radius[edit_ctrlpts_idx] - 10, MIN_RADIUS);
			AssignPoints();
		}
	}
	redisplay_all();
}

void xy_motion(int x, int y)
{
	mouse_plane = 0;
	mouse_x = x;
	mouse_y = y;
	if (edit_ctrlpts_idx != -1)
	{
		int new_x = 4 * x - 2 * subWidth;
		int new_y = 2 * subHeight - 4 * y;
		curve.control_pts[edit_ctrlpts_idx][0] = static_cast<REAL>(new_x);
		curve.control_pts[edit_ctrlpts_idx][1] = static_cast<REAL>(new_y);
	}
	if (isDampled == false)
	{
		AssignPoints();
	}
	redisplay_all();
}

void yz_mouse(int button, int action, int x, int y)
{
	if (GLUT_LEFT_BUTTON == button)
	{
		switch (action)
		{
		case GLUT_DOWN:
			edit_ctrlpts_idx = hit_index(&curve, x, subHeight - y, 1);
			if (isDampled && edit_ctrlpts_idx != -1)
			{
				damplingPointIndex = edit_ctrlpts_idx;
				damplingInitPoint = Vector3(curve.control_pts[edit_ctrlpts_idx][0], curve.control_pts[edit_ctrlpts_idx][1], curve.control_pts[edit_ctrlpts_idx][2]);
			}
			break;
		case GLUT_UP:
			if (isDampled && edit_ctrlpts_idx != -1)
			{
				Vector3 damplingEndPoint = Vector3(curve.control_pts[edit_ctrlpts_idx][0], curve.control_pts[edit_ctrlpts_idx][1], curve.control_pts[edit_ctrlpts_idx][2]);
				damplingVector = damplingEndPoint - damplingInitPoint;
				damplingSecond = 0;
			}
			edit_ctrlpts_idx = -1;
			break;
		default: break;
		}
	}
	else if (button == 3) // wheel up
	{
		edit_ctrlpts_idx = hit_index(&curve, x, subHeight - y, 1);
		if (edit_ctrlpts_idx != -1)
		{
			curve.Radius[edit_ctrlpts_idx] = std::min(curve.Radius[edit_ctrlpts_idx] + 10, MAX_RADIUS);
			AssignPoints();
		}
	}
	else if (button == 4) // wheel down
	{
		edit_ctrlpts_idx = hit_index(&curve, x, subHeight - y, 1);
		if (edit_ctrlpts_idx != -1)
		{
			curve.Radius[edit_ctrlpts_idx] = std::max(curve.Radius[edit_ctrlpts_idx] - 10, MIN_RADIUS);
			AssignPoints();
		}
	}
	redisplay_all();
}

void yz_motion(int x, int y)
{
	mouse_plane = 1;
	mouse_x = x;
	mouse_y = y;
	if (edit_ctrlpts_idx != -1)
	{
		int new_y = 4 * x - 2 * subWidth;
		int new_z = 2 * subHeight - 4 * y;
		curve.control_pts[edit_ctrlpts_idx][1] = static_cast<REAL>(new_y);
		curve.control_pts[edit_ctrlpts_idx][2] = static_cast<REAL>(new_z);
	}
	if (isDampled == false)
	{
		AssignPoints();
	}
	redisplay_all();
}

void zx_mouse(int button, int action, int x, int z)
{
	if (GLUT_LEFT_BUTTON == button)
	{
		switch (action)
		{
		case GLUT_DOWN:
			edit_ctrlpts_idx = hit_index(&curve, x, subHeight - z, 2);
			if (isDampled && edit_ctrlpts_idx != -1)
			{
				damplingPointIndex = edit_ctrlpts_idx;
				damplingInitPoint = Vector3(curve.control_pts[edit_ctrlpts_idx][0], curve.control_pts[edit_ctrlpts_idx][1], curve.control_pts[edit_ctrlpts_idx][2]);
			}
			break;
		case GLUT_UP:
			if (isDampled && edit_ctrlpts_idx != -1)
			{
				Vector3 damplingEndPoint = Vector3(curve.control_pts[edit_ctrlpts_idx][0], curve.control_pts[edit_ctrlpts_idx][1], curve.control_pts[edit_ctrlpts_idx][2]);
				damplingVector = damplingEndPoint - damplingInitPoint;
				damplingSecond = 0;
			}
			edit_ctrlpts_idx = -1;
			break;
		default: break;
		}
	}
	else if (button == 3) // wheel up
	{
		edit_ctrlpts_idx = hit_index(&curve, x, subHeight - z, 2);
		if (edit_ctrlpts_idx != -1)
		{
			curve.Radius[edit_ctrlpts_idx] = std::min(curve.Radius[edit_ctrlpts_idx] + 10, MAX_RADIUS);
			AssignPoints();
		}
	}
	else if (button == 4) // wheel down
	{
		edit_ctrlpts_idx = hit_index(&curve, x, subHeight - z, 2);
		if (edit_ctrlpts_idx != -1)
		{
			curve.Radius[edit_ctrlpts_idx] = std::max(curve.Radius[edit_ctrlpts_idx] - 10, MIN_RADIUS);
			AssignPoints();
		}
	}
	redisplay_all();
}

void zx_motion(int x, int y)
{
	mouse_plane = 2;
	mouse_x = x;
	mouse_y = y;
	if (edit_ctrlpts_idx != -1)
	{
		int new_x = 4 * x - 2 * subWidth;
		int new_y = 4 * y - 2 * subHeight;
		curve.control_pts[edit_ctrlpts_idx][0] = static_cast<REAL>(new_x);
		curve.control_pts[edit_ctrlpts_idx][2] = static_cast<REAL>(new_y);
	}
	if (isDampled == false)
	{
		AssignPoints();
	}
	redisplay_all();
}

void main_reshape(int newWidth, int newHeight)
{
	width = newWidth;
	height = newHeight;
	subWidth = width / 2;
	subHeight = height / 2;

	glutSetWindow(xyzwin);
	glutPositionWindow(subWidth, 0);
	glutReshapeWindow(subWidth, subHeight);

	glutSetWindow(xywin);
	glutPositionWindow(0, 0);
	glutReshapeWindow(subWidth, subHeight);

	glutSetWindow(yzwin);
	glutPositionWindow(0, subHeight);
	glutReshapeWindow(subWidth, subHeight);

	glutSetWindow(zxwin);
	glutPositionWindow(subWidth, subHeight);
	glutReshapeWindow(subWidth, subHeight);
}

void main_display()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glutSwapBuffers();
}

void redisplay_all()
{
	glutSetWindow(xywin);
	SetDisplayXY(subWidth, subHeight);
	glutPostRedisplay();
	glutSetWindow(yzwin);
	SetDisplayYZ(subWidth, subHeight);
	glutPostRedisplay();
	glutSetWindow(zxwin);
	SetDisplayZX(subWidth, subHeight);
	glutPostRedisplay();
}

int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowSize(width, height);
	glutInitWindowPosition(50, 50);

	init();

	window = glutCreateWindow("Bezier Curve in 3D");
	glClearColor(1, 1, 1, 1);
	glutReshapeFunc(main_reshape);
	glutDisplayFunc(main_display);

	xyzwin = glutCreateSubWindow(window, subWidth, 0, subWidth, subHeight);
	initPNG(&texture, texture_filename.c_str(), texture_width, texture_height); // should be called after glutCreateWindow
	glClearColor(1, 1, 1, 1);
	glutReshapeFunc(SetDisplayXYZ);
	glutDisplayFunc(DisplayXYZ);
	glutKeyboardFunc(keyboard_callback);

	xywin = glutCreateSubWindow(window, 0, 0, subWidth, subHeight);
	glClearColor(1, 1, 1, 1);
	glutReshapeFunc(SetDisplayXY);
	glutDisplayFunc(DisplayXY);
	glutKeyboardFunc(keyboard_callback);
	glutMouseFunc(xy_mouse);
	glutMotionFunc(xy_motion);
	glutPassiveMotionFunc(xy_motion);

	yzwin = glutCreateSubWindow(window, 0, subHeight, subWidth, subHeight);
	glClearColor(1, 1, 1, 1);
	glutReshapeFunc(SetDisplayYZ);
	glutDisplayFunc(DisplayYZ);
	glutKeyboardFunc(keyboard_callback);
	glutMouseFunc(yz_mouse);
	glutMotionFunc(yz_motion);
	glutPassiveMotionFunc(yz_motion);

	zxwin = glutCreateSubWindow(window, subWidth, subHeight, subWidth, subHeight);
	glClearColor(1, 1, 1, 1);
	glutReshapeFunc(SetDisplayZX);
	glutDisplayFunc(DisplayZX);
	glutKeyboardFunc(keyboard_callback);
	glutMouseFunc(zx_mouse);
	glutMotionFunc(zx_motion);
	glutPassiveMotionFunc(zx_motion);

	glutMainLoop();
	return 0;
}