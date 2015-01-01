#include <GL/glut.h>
#include "curve.h"
#include <vector>
#include <list>
#include <algorithm>
#include <iostream>
#include <limits>

CubicBezierCurve curve[2];
GLsizei width = 640, height = 800;
int edit_ctrlpts_idx = -1;
bool isDrawControlMesh = true;
bool isDottedLine = false;
bool isAABB = false;
bool isOBB = false;
bool isAABBIntersection = false;
bool isOBBIntersection = false;
bool isAABBSelfIntersection = false;
bool isOBBSelfIntersection = false;
typedef std::vector<std::vector<double>> ConstantsDictType;
ConstantsDictType BeizerConstants;
typedef std::pair<double, double> PointType;
typedef std::vector<PointType> PointsType;
PointsType Points[2];
#define RES 1024

int hit_index(CubicBezierCurve* curve, int x, int y)
{
	for (int i = 0; i < 4; i++)
	{
		REAL tx = curve[0].control_pts[i][0] - x;
		REAL ty = curve[0].control_pts[i][1] - y;
		if ((tx * tx + ty * ty) < 30) return i;
	}
	for (int i = 0; i < 4; i++)
	{
		REAL tx = curve[1].control_pts[i][0] - x;
		REAL ty = curve[1].control_pts[i][1] - y;
		if ((tx * tx + ty * ty) < 30) return i + 4;
	}
	return -1;
}

void init()
{
	SET_PT2(curve[0].control_pts[0], 50, 100);
	SET_PT2(curve[0].control_pts[1], 200, 300);
	SET_PT2(curve[0].control_pts[2], 400, 300);
	SET_PT2(curve[0].control_pts[3], 550, 100);
	SET_PT2(curve[1].control_pts[0], 50, 400);
	SET_PT2(curve[1].control_pts[1], 200, 600);
	SET_PT2(curve[1].control_pts[2], 400, 600);
	SET_PT2(curve[1].control_pts[3], 550, 400);
	BeizerConstants = std::vector<std::vector<double>>(1025);
	Points[0] = std::vector<PointType>(1025);
	Points[1] = std::vector<PointType>(1025);
	for (int i = 0; i <= RES; ++i)
	{
		double t = static_cast<double>(i) / RES;
		double one_minus_t = 1 - t;
		BeizerConstants[i] = std::vector<double>(4);
		BeizerConstants[i][0] = one_minus_t * one_minus_t * one_minus_t;
		BeizerConstants[i][1] = 3 * one_minus_t * one_minus_t * t;
		BeizerConstants[i][2] = 3 * one_minus_t * t * t;
		BeizerConstants[i][3] = t * t * t;
	}
	for (int i = 0; i <= RES; ++i)
	{
		Points[0][i] = PointType();
		for (int j = 0; j < 4; ++j)
		{
			Points[0][i].first += BeizerConstants[i][j] * curve[0].control_pts[j][0];
			Points[0][i].second += BeizerConstants[i][j] * curve[0].control_pts[j][1];
		}
	}
	for (int i = 0; i <= RES; ++i)
	{
		Points[1][i] = PointType();
		for (int j = 0; j < 4; ++j)
		{
			Points[1][i].first += BeizerConstants[i][j] * curve[1].control_pts[j][0];
			Points[1][i].second += BeizerConstants[i][j] * curve[1].control_pts[j][1];
		}
	}
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glMatrixMode(GL_PROJECTION);
	gluOrtho2D(0, width, 0, height);
}

void reshape_callback(GLint nw, GLint nh)
{
	width = nw;
	height = nh;
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0, width, 0, height);
}

void draw_AABB(const PointsType& points, const unsigned int low, const unsigned int high)
{
	if (high <= low)
	{
		return;
	}

	double left = points[low].first;
	double right = points[low].first;
	double top = points[low].second;
	double bottom = points[low].second;
	for (unsigned int i = low + 1; i <= high; ++i)
	{
		if (points[i].first < left)
		{
			left = points[i].first;
		}
		if (right < points[i].first)
		{
			right = points[i].first;
		}
		if (top < points[i].second)
		{
			top = points[i].second;
		}
		if (points[i].second < bottom)
		{
			bottom = points[i].second;
		}
	}
	glColor3ub(0, 255, 0);
	glBegin(GL_LINE_LOOP);
	glVertex2d(left, top);
	glVertex2d(left, bottom);
	glVertex2d(right, bottom);
	glVertex2d(right, top);
	glEnd();

	if (high - low == 1)
	{
		return;
	}

	unsigned int average = (low + high) / 2;
	draw_AABB(points, low, average);
	draw_AABB(points, average, high);
}

inline double InnerProduct(const PointType& lhs, const PointType& rhs)
{
	return lhs.first * rhs.first + lhs.second * rhs.second;
}

PointsType GetOBBPoints(const PointsType& points, const unsigned int low, const unsigned int high)
{
	if (low == high)
	{
		PointsType OBBPoints;
		OBBPoints.push_back(points[low]);
		OBBPoints.push_back(points[low]);
		OBBPoints.push_back(points[low]);
		OBBPoints.push_back(points[low]);
		return OBBPoints;
	}
	PointType vect = std::make_pair(points[high].first - points[low].first, points[high].second - points[low].second);
	PointType orthoVect = std::make_pair(vect.second, -vect.first);
	unsigned int leftIndex = low;
	unsigned int rightIndex = low;
	unsigned int topIndex = low;
	unsigned int bottomIndex = low;
	double left = InnerProduct(points[low], vect);
	double right = InnerProduct(points[low], vect);
	double top = InnerProduct(points[low], orthoVect);
	double bottom = InnerProduct(points[low], orthoVect);
	for (unsigned int i = low + 1; i <= high; ++i)
	{
		double vectValue = InnerProduct(points[i], vect);
		double orthoVectValue = InnerProduct(points[i], orthoVect);
		if (vectValue < left)
		{
			left = vectValue;
			leftIndex = i;
		}
		if (right < vectValue)
		{
			right = vectValue;
			rightIndex = i;
		}
		if (top < orthoVectValue)
		{
			top = orthoVectValue;
			topIndex = i;
		}
		if (orthoVectValue < bottom)
		{
			bottom = orthoVectValue;
			bottomIndex = i;
		}
	}
	double x1 = (orthoVect.second * vect.first * points[leftIndex].first + orthoVect.second * vect.second * points[leftIndex].second - vect.second * orthoVect.first * points[topIndex].first - vect.second * orthoVect.second * points[topIndex].second) / (orthoVect.second * vect.first - vect.second * orthoVect.first);
	double y1 = (orthoVect.first * vect.first * points[leftIndex].first + orthoVect.first * vect.second * points[leftIndex].second - vect.first * orthoVect.first * points[topIndex].first - vect.first * orthoVect.second * points[topIndex].second) / (orthoVect.first * vect.second - vect.first * orthoVect.second);
	double x2 = (orthoVect.second * vect.first * points[leftIndex].first + orthoVect.second * vect.second * points[leftIndex].second - vect.second * orthoVect.first * points[bottomIndex].first - vect.second * orthoVect.second * points[bottomIndex].second) / (orthoVect.second * vect.first - vect.second * orthoVect.first);
	double y2 = (orthoVect.first * vect.first * points[leftIndex].first + orthoVect.first * vect.second * points[leftIndex].second - vect.first * orthoVect.first * points[bottomIndex].first - vect.first * orthoVect.second * points[bottomIndex].second) / (orthoVect.first * vect.second - vect.first * orthoVect.second);
	double x3 = (orthoVect.second * vect.first * points[rightIndex].first + orthoVect.second * vect.second * points[rightIndex].second - vect.second * orthoVect.first * points[bottomIndex].first - vect.second * orthoVect.second * points[bottomIndex].second) / (orthoVect.second * vect.first - vect.second * orthoVect.first);
	double y3 = (orthoVect.first * vect.first * points[rightIndex].first + orthoVect.first * vect.second * points[rightIndex].second - vect.first * orthoVect.first * points[bottomIndex].first - vect.first * orthoVect.second * points[bottomIndex].second) / (orthoVect.first * vect.second - vect.first * orthoVect.second);
	double x4 = (orthoVect.second * vect.first * points[rightIndex].first + orthoVect.second * vect.second * points[rightIndex].second - vect.second * orthoVect.first * points[topIndex].first - vect.second * orthoVect.second * points[topIndex].second) / (orthoVect.second * vect.first - vect.second * orthoVect.first);
	double y4 = (orthoVect.first * vect.first * points[rightIndex].first + orthoVect.first * vect.second * points[rightIndex].second - vect.first * orthoVect.first * points[topIndex].first - vect.first * orthoVect.second * points[topIndex].second) / (orthoVect.first * vect.second - vect.first * orthoVect.second);
	PointsType OBBPoints;
	OBBPoints.push_back(PointType(x1, y1));
	OBBPoints.push_back(PointType(x2, y2));
	OBBPoints.push_back(PointType(x3, y3));
	OBBPoints.push_back(PointType(x4, y4));
	return OBBPoints;
}

void draw_OBB(const PointsType& points, const unsigned int low, const unsigned int high)
{
	if (high <= low)
	{
		return;
	}
	PointsType OBBPoints = GetOBBPoints(points, low, high);
	glColor3ub(0, 255, 0);
	glBegin(GL_LINE_LOOP);
	for (int i = 0; i < 4; ++i)
	{
		glVertex2d(OBBPoints[i].first, OBBPoints[i].second);
	}
	glEnd();

	if (high - low == 1)
	{
		return;
	}

	unsigned int average = (low + high) / 2;
	draw_OBB(points, low, average);
	draw_OBB(points, average, high);
}

bool GetIntersectPoint(const PointType& point1, const PointType& point2, const PointType& point3, const PointType& point4, PointType& intersectPoint)
{
	double t;
	double s;
	double under = (point4.second - point3.second) * (point2.first - point1.first) - (point4.first - point3.first) * (point2.second - point1.second);
	if (under == 0)
	{
		return false;
	}
	double _t = (point4.first - point3.first) * (point1.second - point3.second) - (point4.second - point3.second) * (point1.first - point3.first);
	double _s = (point2.first - point1.first) * (point1.second - point3.second) - (point2.second - point1.second) * (point1.first - point3.first);
	t = _t / under;
	s = _s / under;
	if (t < 0.0 || t > 1.0 || s < 0.0 || s > 1.0)
	{
		return false;
	}
	else if (_t == 0 && _s == 0)
	{
		return false;
	}
	intersectPoint.first = point1.first + t * (point2.first - point1.first);
	intersectPoint.second = point1.second + t * (point2.second - point1.second);
	return true;
}

inline double GetRange(const PointType& lhs, const PointType& rhs)
{
	return std::sqrt(std::pow(rhs.first - lhs.first, 2) + std::pow(rhs.second - lhs.second, 2));
}

void draw_AABB_intersection(const PointsType& lhs, const PointsType& rhs, const unsigned int lhs_low, const unsigned int lhs_high, const unsigned int rhs_low, const unsigned int rhs_high)
{
	if (lhs_high <= lhs_low)
	{
		return;
	}
	else if (rhs_high <= rhs_low)
	{
		return;
	}

	double lhs_left = lhs[lhs_low].first;
	double lhs_right = lhs[lhs_low].first;
	double lhs_top = lhs[lhs_low].second;
	double lhs_bottom = lhs[lhs_low].second;
	for (unsigned int i = lhs_low + 1; i <= lhs_high; ++i)
	{
		if (lhs[i].first < lhs_left)
		{
			lhs_left = lhs[i].first;
		}
		if (lhs_right < lhs[i].first)
		{
			lhs_right = lhs[i].first;
		}
		if (lhs_top < lhs[i].second)
		{
			lhs_top = lhs[i].second;
		}
		if (lhs[i].second < lhs_bottom)
		{
			lhs_bottom = lhs[i].second;
		}
	}

	glColor3ub(0, 255, 0);
	glBegin(GL_LINE_LOOP);
	glVertex2d(lhs_left, lhs_top);
	glVertex2d(lhs_left, lhs_bottom);
	glVertex2d(lhs_right, lhs_bottom);
	glVertex2d(lhs_right, lhs_top);
	glEnd();

	double rhs_left = rhs[rhs_low].first;
	double rhs_right = rhs[rhs_low].first;
	double rhs_top = rhs[rhs_low].second;
	double rhs_bottom = rhs[rhs_low].second;
	for (unsigned int i = rhs_low + 1; i <= rhs_high; ++i)
	{
		if (rhs[i].first < rhs_left)
		{
			rhs_left = rhs[i].first;
		}
		if (rhs_right < rhs[i].first)
		{
			rhs_right = rhs[i].first;
		}
		if (rhs_top < rhs[i].second)
		{
			rhs_top = rhs[i].second;
		}
		if (rhs[i].second < rhs_bottom)
		{
			rhs_bottom = rhs[i].second;
		}
	}

	glBegin(GL_LINE_LOOP);
	glVertex2d(rhs_left, rhs_top);
	glVertex2d(rhs_left, rhs_bottom);
	glVertex2d(rhs_right, rhs_bottom);
	glVertex2d(rhs_right, rhs_top);
	glEnd();

	if (lhs_right < rhs_left)
	{
		return;
	}
	else if (rhs_right < lhs_left)
	{
		return;
	}
	else if (lhs_top < rhs_bottom)
	{
		return;
	}
	else if (rhs_top < lhs_bottom)
	{
		return;
	}

	// Bounding Volume °ãÄ§
	if ((lhs_high - lhs_low) == 1 && (rhs_high - rhs_low) == 1)
	{
		const PointType& x1_y1 = lhs[lhs_low];
		const PointType& x2_y2 = lhs[lhs_high];
		const PointType& x3_y3 = rhs[rhs_low];
		const PointType& x4_y4 = rhs[rhs_high];
		double x1 = x1_y1.first, y1 = x1_y1.second;
		double x2 = x2_y2.first, y2 = x2_y2.second;
		double x3 = x3_y3.first, y3 = x3_y3.second;
		double x4 = x4_y4.first, y4 = x4_y4.second;
		double denominator = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
		if (denominator == 0)
		{
			// ÆòÇà
			glColor3ub(255, 0, 0);
			glPointSize(10.0);
			glBegin(GL_POINTS);
			glVertex2d((x1 + x2 + x3 + x4) / 4, (y1 + y2 + y3 + y4) / 4);
			glEnd();
		}
		else
		{
			PointType intersectPoint;
			if (GetIntersectPoint(x1_y1, x2_y2, x3_y3, x4_y4, intersectPoint))
			{
				glColor3ub(255, 0, 0);
				glPointSize(10.0);
				glBegin(GL_POINTS);
				glVertex2d(intersectPoint.first, intersectPoint.second);
				glEnd();
			}
		}
		return;
	}

	if ((lhs_high - lhs_low) > (rhs_high - rhs_low))
	{
		unsigned int average = (lhs_low + lhs_high) / 2;
		draw_AABB_intersection(lhs, rhs, lhs_low, average, rhs_low, rhs_high);
		draw_AABB_intersection(lhs, rhs, average, lhs_high, rhs_low, rhs_high);
	}
	else
	{
		unsigned int average = (rhs_low + rhs_high) / 2;
		draw_AABB_intersection(lhs, rhs, lhs_low, lhs_high, rhs_low, average);
		draw_AABB_intersection(lhs, rhs, lhs_low, lhs_high, average, rhs_high);
	}
}

bool IsIntersected(const PointType& basePoint1, const PointType& basePoint2, const PointType& basePoint3, const PointsType& checkPoints)
{
	PointType vect = PointType(basePoint3.first - basePoint2.first, basePoint3.second - basePoint2.second);
	double vectLength = std::sqrt(InnerProduct(vect, vect));
	PointType point = PointType((basePoint1.first + basePoint2.first) / 2, (basePoint1.second + basePoint2.second) / 2);
	double r = GetRange(basePoint1, point);
	double prev = 0;
	for (int i = 0; i < 4; ++i)
	{
		double curr = vect.second * checkPoints[i].first - vect.first * checkPoints[i].second + vect.first * point.second - vect.second * point.first;
		if (std::abs(curr) / vectLength <= r)
		{
			return true;
		}
		else if (prev * curr < 0)
		{
			return true;
		}
		prev = curr;
	}
	return false;
}

bool IsIntersected(const PointsType& lhsOBBPoints, const PointsType& rhsOBBPoints)
{
	double range1 = GetRange(lhsOBBPoints[0], lhsOBBPoints[1]);
	double range2 = GetRange(lhsOBBPoints[1], lhsOBBPoints[2]);
	double range3 = GetRange(rhsOBBPoints[0], rhsOBBPoints[1]);
	double range4 = GetRange(rhsOBBPoints[1], rhsOBBPoints[2]);
	std::vector<double> ranges;
	ranges.push_back(range1);
	ranges.push_back(range2);
	ranges.push_back(range3);
	ranges.push_back(range4);
	std::sort(ranges.begin(), ranges.end());
	for (int i = 0; i < 4; ++i)
	{
		double range = ranges[i];
		if (range == range1)
		{
			if (IsIntersected(lhsOBBPoints[0], lhsOBBPoints[1], lhsOBBPoints[2], rhsOBBPoints) == false)
			{
				return false;
			}
		}
		else if (range == range2)
		{
			if (IsIntersected(lhsOBBPoints[1], lhsOBBPoints[2], lhsOBBPoints[0], rhsOBBPoints) == false)
			{
				return false;
			}
		}
		else if (range == range3)
		{
			if (IsIntersected(rhsOBBPoints[0], rhsOBBPoints[1], rhsOBBPoints[2], lhsOBBPoints) == false)
			{
				return false;
			}
		}
		else if (range == range4)
		{
			if (IsIntersected(rhsOBBPoints[1], rhsOBBPoints[2], rhsOBBPoints[0], lhsOBBPoints) == false)
			{
				return false;
			}
		}
	}
	return true;
}

void DrawOBBIntersection(const PointsType& lhs, const PointsType& rhs, const unsigned int lhsLow, const unsigned int lhsHigh, const unsigned int rhsLow, const unsigned int rhsHigh)
{
	if (lhsHigh <= lhsLow)
	{
		return;
	}
	else if (rhsHigh <= rhsLow)
	{
		return;
	}

	PointsType lhsOBBPoints = GetOBBPoints(lhs, lhsLow, lhsHigh);
	glColor3ub(0, 255, 0);
	glBegin(GL_LINE_LOOP);
	for (int i = 0; i < 4; ++i)
	{
		glVertex2d(lhsOBBPoints[i].first, lhsOBBPoints[i].second);
	}
	glEnd();

	PointsType rhsOBBPoints = GetOBBPoints(rhs, rhsLow, rhsHigh);
	glBegin(GL_LINE_LOOP);
	for (int i = 0; i < 4; ++i)
	{
		glVertex2d(rhsOBBPoints[i].first, rhsOBBPoints[i].second);
	}
	glEnd();

	// Bounding Volumne °ãÄ¡´ÂÁö È®ÀÎ
	if (IsIntersected(lhsOBBPoints, rhsOBBPoints) == false)
	{
		return;
	}

	// Bounding Volume °ãÄ§
	if ((lhsHigh - lhsLow) == 1 && (rhsHigh - rhsLow) == 1)
	{
		const PointType& x1_y1 = lhs[lhsLow];
		const PointType& x2_y2 = lhs[lhsHigh];
		const PointType& x3_y3 = rhs[rhsLow];
		const PointType& x4_y4 = rhs[rhsHigh];
		double x1 = x1_y1.first, y1 = x1_y1.second;
		double x2 = x2_y2.first, y2 = x2_y2.second;
		double x3 = x3_y3.first, y3 = x3_y3.second;
		double x4 = x4_y4.first, y4 = x4_y4.second;
		double denominator = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
		if (denominator == 0)
		{
			// ÆòÇà
			glColor3ub(255, 0, 0);
			glPointSize(10.0);
			glBegin(GL_POINTS);
			glVertex2d((x1 + x2 + x3 + x4) / 4, (y1 + y2 + y3 + y4) / 4);
			glEnd();
		}
		else
		{
			PointType intersectPoint;
			if (GetIntersectPoint(x1_y1, x2_y2, x3_y3, x4_y4, intersectPoint))
			{
				glColor3ub(255, 0, 0);
				glPointSize(10.0);
				glBegin(GL_POINTS);
				glVertex2d(intersectPoint.first, intersectPoint.second);
				glEnd();
			}
		}
		return;
	}

	if ((lhsHigh - lhsLow) > (rhsHigh - rhsLow))
	{
		unsigned int average = (lhsLow + lhsHigh) / 2;
		DrawOBBIntersection(lhs, rhs, lhsLow, average, rhsLow, rhsHigh);
		DrawOBBIntersection(lhs, rhs, average, lhsHigh, rhsLow, rhsHigh);
	}
	else
	{
		unsigned int average = (rhsLow + rhsHigh) / 2;
		DrawOBBIntersection(lhs, rhs, lhsLow, lhsHigh, rhsLow, average);
		DrawOBBIntersection(lhs, rhs, lhsLow, lhsHigh, average, rhsHigh);
	}
}

void DrawAABBSelfIntersection(const PointsType& points, const unsigned int low, const unsigned int high)
{
	if (high - low <= 1)
	{
		return;
	}
	unsigned int average = (low + high) / 2;
	DrawAABBSelfIntersection(points, low, average);
	DrawAABBSelfIntersection(points, average, high);
	draw_AABB_intersection(points, points, low, average - 1, average + 1, high);
}

void DrawAABBSelfIntersection(const PointsType& lhs, const PointsType& rhs)
{
	draw_AABB_intersection(lhs, rhs, 0, RES, 0, RES);
	DrawAABBSelfIntersection(lhs, 0, RES);
	DrawAABBSelfIntersection(rhs, 0, RES);
}

void DrawOBBSelfIntersection(const PointsType& points, const unsigned int low, const unsigned int high)
{
	if (high - low <= 1)
	{
		return;
	}
	unsigned int average = (low + high) / 2;
	DrawOBBSelfIntersection(points, low, average);
	DrawOBBSelfIntersection(points, average, high);
	DrawOBBIntersection(points, points, low, average - 1, average + 1, high);
}

void DrawOBBSelfIntersection(const PointsType& lhs, const PointsType& rhs)
{
	DrawOBBIntersection(lhs, rhs, 0, RES, 0, RES);
	DrawOBBSelfIntersection(lhs, 0, RES);
	DrawOBBSelfIntersection(rhs, 0, RES);
}

void display_callback()
{
	/* curve */
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3ub(0, 0, 0);
	if (isDottedLine)
		glBegin(GL_LINES);
	else
		glBegin(GL_LINE_STRIP);
	for (int i = 0; i <= RES; i++)
	{
		glVertex2d(Points[0][i].first, Points[0][i].second);
	}
	glEnd();
	if (isDottedLine)
		glBegin(GL_LINES);
	else
		glBegin(GL_LINE_STRIP);
	for (int i = 0; i <= RES; i++)
	{
		glVertex2d(Points[1][i].first, Points[1][i].second);
	}
	glEnd();

	// draw AABB
	if (isAABB)
	{
		draw_AABB(Points[0], 0, RES);
		draw_AABB(Points[1], 0, RES);
	}

	// draw OBB
	if (isOBB)
	{
		draw_OBB(Points[0], 0, RES);
		draw_OBB(Points[1], 0, RES);
	}

	if (isAABBIntersection)
	{
		draw_AABB_intersection(Points[0], Points[1], 0, RES, 0, RES);
	}

	if (isOBBIntersection)
	{
		DrawOBBIntersection(Points[0], Points[1], 0, RES, 0, RES);
	}

	if (isAABBSelfIntersection)
	{
		DrawAABBSelfIntersection(Points[0], Points[1]);
	}

	if (isOBBSelfIntersection)
	{
		DrawOBBSelfIntersection(Points[0], Points[1]);
	}

	/* control mesh */
	if (isDrawControlMesh)
	{
		glColor3ub(255, 0, 0);
		glBegin(GL_LINE_STRIP);
		for (int i = 0; i < 4; i++)
		{
			REAL *pt = curve[0].control_pts[i];
			glVertex2d(pt[0], pt[1]);
		}
		glEnd();
		glBegin(GL_LINE_STRIP);
		for (int i = 0; i < 4; i++)
		{
			REAL *pt = curve[1].control_pts[i];
			glVertex2d(pt[0], pt[1]);
		}
		glEnd();
	}

	/* control pts */
	glColor3ub(0, 0, 255);
	glPointSize(10.0);
	glBegin(GL_POINTS);
	for (int i = 0; i < 4; i++)
	{
		REAL *pt = curve[0].control_pts[i];
		glVertex2d(pt[0], pt[1]);
	}
	glEnd();
	glBegin(GL_POINTS);
	for (int i = 0; i < 4; i++)
	{
		REAL *pt = curve[1].control_pts[i];
		glVertex2d(pt[0], pt[1]);
	}
	glEnd();
	glutSwapBuffers();
}

// void glutMouseFunc(void (*func)(int button, int state, int x, int y));
void mouse_callback(GLint button, GLint action, GLint x, GLint y)
{
	if (GLUT_LEFT_BUTTON == button)
	{
		switch (action)
		{
		case GLUT_DOWN:
			edit_ctrlpts_idx = hit_index(curve, x, height - y);
			break;
		case GLUT_UP:
			edit_ctrlpts_idx = -1;
			break;
		default: break;
		}
	}
	glutPostRedisplay();
}

// void glutMotionFunc(void (*func)(int x, int y));
void mouse_move_callback(GLint x, GLint y)
{
	if (edit_ctrlpts_idx != -1)
	{
		int curve_index = edit_ctrlpts_idx < 4 ? 0 : 1;
		int control_pts_index = edit_ctrlpts_idx < 4 ? edit_ctrlpts_idx : edit_ctrlpts_idx - 4;
		curve[curve_index].control_pts[control_pts_index][0] = static_cast<double>(x);
		curve[curve_index].control_pts[control_pts_index][1] = static_cast<double>(height - y);
		for (int i = 0; i <= RES; ++i)
		{
			Points[curve_index][i] = PointType();
			for (int j = 0; j < 4; ++j)
			{
				Points[curve_index][i].first += BeizerConstants[i][j] * curve[curve_index].control_pts[j][0];
				Points[curve_index][i].second += BeizerConstants[i][j] * curve[curve_index].control_pts[j][1];
			}
		}
	}
	glutPostRedisplay();
}

// void glutKeyboardFunc(void (*func)(unsigned char key, int x, int y));
void keyboard_callback(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'i': case 'I':
		SET_PT2(curve[0].control_pts[0], 50, 100);
		SET_PT2(curve[0].control_pts[1], 200, 300);
		SET_PT2(curve[0].control_pts[2], 400, 300);
		SET_PT2(curve[0].control_pts[3], 550, 100);
		SET_PT2(curve[1].control_pts[0], 50, 400);
		SET_PT2(curve[1].control_pts[1], 200, 600);
		SET_PT2(curve[1].control_pts[2], 400, 600);
		SET_PT2(curve[1].control_pts[3], 550, 400);
		for (int i = 0; i <= RES; ++i)
		{
			Points[0][i] = PointType();
			for (int j = 0; j < 4; ++j)
			{
				Points[0][i].first += BeizerConstants[i][j] * curve[0].control_pts[j][0];
				Points[0][i].second += BeizerConstants[i][j] * curve[0].control_pts[j][1];
			}
		}
		for (int i = 0; i <= RES; ++i)
		{
			Points[1][i] = PointType();
			for (int j = 0; j < 4; ++j)
			{
				Points[1][i].first += BeizerConstants[i][j] * curve[1].control_pts[j][0];
				Points[1][i].second += BeizerConstants[i][j] * curve[1].control_pts[j][1];
			}
		}
		break;
	case 'l': case 'L':
		isDottedLine ^= true;
		break;
	case 'c': case 'C':
		isDrawControlMesh ^= true;
		break;
	case 'a':
	case 'A':
		isAABB ^= true;
		break;
	case 'o':
	case 'O':
		isOBB ^= true;
		break;
	case 's':
	case 'S':
		isAABBIntersection ^= true;
		break;
	case 'b':
	case 'B':
		isOBBIntersection ^= true;
		break;
	case 'p':
	case 'P':
		isAABBSelfIntersection ^= true;
		break;
	case 'q':
	case 'Q':
		isOBBSelfIntersection ^= true;
		break;
	case (27) : exit(0); break;
	default: break;
	}
	glutPostRedisplay();
}

int main(int argc, char *argv[])
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowSize(width, height);
	glutCreateWindow("Beizer Editor");

	init();
	glutReshapeFunc(reshape_callback);
	glutMouseFunc(mouse_callback);
	glutMotionFunc(mouse_move_callback);
	glutDisplayFunc(display_callback);
	glutKeyboardFunc(keyboard_callback);
	glutMainLoop();

	return 0;
}
