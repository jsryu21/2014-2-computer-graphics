#include <GL/freeglut.h>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <assert.h>


// lineintersection
void IntersectLineSegments(const float A1x, const float A1y, const float A1z,
                           const float A2x, const float A2y, const float A2z, 
                           const float B1x, const float B1y, const float B1z,
                           const float B2x, const float B2y, const float B2z,
                           bool infinite_lines, float epsilon, float &PointOnSegAx,
                           float &PointOnSegAy, float &PointOnSegAz, float &PointOnSegBx,
                           float &PointOnSegBy, float &PointOnSegBz, float &NearestPointX,
                           float &NearestPointY, float &NearestPointZ, float &NearestVectorX,
                           float &NearestVectorY, float &NearestVectorZ, bool &true_intersection);

void FindNearestPointOnLineSegment(const float A1x, const float A1y, const float A1z,
                                   const float Lx, const float Ly, const float Lz,
                                   const float Bx, const float By, const float Bz,
                                   bool infinite_line, float epsilon_squared, float &NearestPointX,
                                   float &NearestPointY, float &NearestPointZ,
                                   float &parameter);

void FindNearestPointOfParallelLineSegments(float A1x, float A1y, float A1z,
                                            float A2x, float A2y, float A2z,
                                            float Lax, float Lay, float Laz,
                                            float B1x, float B1y, float B1z,
                                            float B2x, float B2y, float B2z,
                                            float Lbx, float Lby, float Lbz,
                                            bool infinite_lines, float epsilon_squared,
                                            float &PointOnSegAx, float &PointOnSegAy, float &PointOnSegAz,
                                            float &PointOnSegBx, float &PointOnSegBy, float &PointOnSegBz);

void AdjustNearestPoints(float A1x, float A1y, float A1z,
                         float Lax, float Lay, float Laz,
                         float B1x, float B1y, float B1z,
                         float Lbx, float Lby, float Lbz,
                         float epsilon_squared, float s, float t,
                         float &PointOnSegAx, float &PointOnSegAy, float &PointOnSegAz,
                         float &PointOnSegBx, float &PointOnSegBy, float &PointOnSegBz);



// uncomment the following line to have the code check intermediate results
//#define CHECK_ANSWERS

// uncomment the following line to use Cramer's rule instead of Gaussian elimination
//#define USE_CRAMERS_RULE

#define FMAX(a,b) ((a) > (b) ? (a) : (b))
#define FMIN(a,b) ((a) > (b) ? (b) : (a))
#define FABS(a) ((a) < 0.0f ? -(a) : (a))
#define OUT_OF_RANGE(a) ((a) < 0.0f || (a) > 1.f)

// pragma to get rid of math.h inline function removal warnings.
#pragma warning(disable:4514)

/**************************************************************************
|
|     Method: IntersectLineSegments
|
|    Purpose: Find the nearest point between two finite length line segments
|             or two infinite lines in 3-dimensional space. The function calculates
|             the point on each line/line segment that is closest to the other
|             line/line segment, the midpoint between the nearest points, and
|             the vector between these two points. If the two nearest points
|             are close within a tolerance, a flag is set indicating the lines
|             have a "true" intersection.
|
| Parameters: Input:
|             ------
|             A1x, A1y, A1z   - Coordinates of first defining point of line/segment A
|             A2x, A2y, A2z   - Coordinates of second defining point of line/segment A
|             B1x, B1y, B1z   - Coordinates of first defining point of line/segment B
|             B2x, B2y, B2z   - Coordinates of second defining point of line/segment B
|             infinite_lines  - set to true if lines are to be treated as infinite
|             epsilon         - tolerance value to be used to check for degenerate
|                               and parallel lines, and to check for true intersection.
|
|             Output:
|             -------
|             PointOnSegAx,   - Coordinates of the point on segment A that are nearest
|             PointOnSegAy,     to segment B. This corresponds to point C in the text.
|             PointOnSegAz
|             PointOnSegBx,   - Coordinates of the point on segment B that are nearest
|             PointOnSegBy,     to segment A. This corresponds to point D in the text.
|             PointOnSegBz
|             NearestPointX,  - Midpoint between the two nearest points. This can be
|             NearestPointY,    treated as *the* intersection point if nearest points
|             NearestPointZ     are sufficiently close. This corresponds to point P
|                               in the text.
|             NearestVectorX, - Vector between the nearest point on A to the nearest
|                               point on segment B. This vector is normal to both
|                               lines if the lines are infinite, but is not guaranteed
|                               to be normal to both lines if both lines are finite
|                               length.
|           true_intersection - true if the nearest points are close within a small
|                               tolerance.
**************************************************************************/
void IntersectLineSegments(const float A1x, const float A1y, const float A1z,
                           const float A2x, const float A2y, const float A2z, 
                           const float B1x, const float B1y, const float B1z,
                           const float B2x, const float B2y, const float B2z,
                           bool infinite_lines, float epsilon, float &PointOnSegAx,
                           float &PointOnSegAy, float &PointOnSegAz, float &PointOnSegBx,
                           float &PointOnSegBy, float &PointOnSegBz, float &NearestPointX,
                           float &NearestPointY, float &NearestPointZ, float &NearestVectorX,
                           float &NearestVectorY, float &NearestVectorZ, bool &true_intersection)
{
  float temp = 0.f;
  float epsilon_squared = epsilon * epsilon;

// Compute parameters from Equations (1) and (2) in the text
  float Lax = A2x - A1x;
  float Lay = A2y - A1y;
  float Laz = A2z - A1z;
  float Lbx = B2x - B1x;
  float Lby = B2y - B1y;
  float Lbz = B2z - B1z;
// From Equation (15)
  float L11 =  (Lax * Lax) + (Lay * Lay) + (Laz * Laz);
  float L22 =  (Lbx * Lbx) + (Lby * Lby) + (Lbz * Lbz);

// Line/Segment A is degenerate ---- Special Case #1
  if (L11 < epsilon_squared)
  {
    PointOnSegAx = A1x;
    PointOnSegAy = A1y;
    PointOnSegAz = A1z;
    FindNearestPointOnLineSegment(B1x, B1y, B1z, Lbx, Lby, Lbz, A1x, A1y, A1z,
                                  infinite_lines, epsilon, PointOnSegBx, PointOnSegBy,
                                  PointOnSegBz, temp);
  }
// Line/Segment B is degenerate ---- Special Case #1
  else if (L22 < epsilon_squared)
  {
    PointOnSegBx = B1x;
    PointOnSegBy = B1y;
    PointOnSegBz = B1z;
    FindNearestPointOnLineSegment(A1x, A1y, A1z, Lax, Lay, Laz, B1x, B1y, B1z,
                                  infinite_lines, epsilon, PointOnSegAx, PointOnSegAy,
                                  PointOnSegAz, temp);
  }
// Neither line/segment is degenerate
  else
  {
// Compute more parameters from Equation (3) in the text.
    float ABx = B1x - A1x;
    float ABy = B1y - A1y;
    float ABz = B1z - A1z;

// and from Equation (15).
    float L12 = -(Lax * Lbx) - (Lay * Lby) - (Laz * Lbz);

    float DetL = L11 * L22 - L12 * L12;
// Lines/Segments A and B are parallel ---- special case #2.
    if (FABS(DetL) < epsilon)
    {
      FindNearestPointOfParallelLineSegments(A1x, A1y, A1z, A2x, A2y, A2z,
                                             Lax, Lay, Laz,
                                             B1x, B1y, B1z, B2x, B2y, B2z,
                                             Lbx, Lby, Lbz,
                                             infinite_lines, epsilon,
                                             PointOnSegAx, PointOnSegAy, PointOnSegAz,
                                             PointOnSegBx, PointOnSegBy, PointOnSegBz);
    }
// The general case
    else
    {
// from Equation (15)
      float ra = Lax * ABx + Lay * ABy + Laz * ABz;
      float rb = -Lbx * ABx - Lby * ABy - Lbz * ABz;

      float t = (L11 * rb - ra * L12)/DetL; // Equation (12)

#ifdef USE_CRAMERS_RULE
      float s = (L22 * ra - rb * L12)/DetL;
#else
      float s = (ra-L12*t)/L11;             // Equation (13)
#endif // USE_CRAMERS_RULE

#ifdef CHECK_ANSWERS
      float check_ra = s*L11 + t*L12;
      float check_rb = s*L12 + t*L22;
      assert(FABS(check_ra-ra) < epsilon);
      assert(FABS(check_rb-rb) < epsilon);
#endif // CHECK_ANSWERS

// if we are dealing with infinite lines or if parameters s and t both
// lie in the range [0,1] then just compute the points using Equations
// (1) and (2) from the text.
      PointOnSegAx = (A1x + s * Lax);
      PointOnSegAy = (A1y + s * Lay);
      PointOnSegAz = (A1z + s * Laz);
      PointOnSegBx = (B1x + t * Lbx);
      PointOnSegBy = (B1y + t * Lby);
      PointOnSegBz = (B1z + t * Lbz);
// otherwise, at least one of s and t is outside of [0,1] and we have to
// handle this case.
      if (false == infinite_lines && (OUT_OF_RANGE(s) || OUT_OF_RANGE(t)))
      {
        AdjustNearestPoints(A1x, A1y, A1z, Lax, Lay, Laz,
                            B1x, B1y, B1z, Lbx, Lby, Lbz,
                            epsilon, s, t,
                            PointOnSegAx, PointOnSegAy, PointOnSegAz,
                            PointOnSegBx, PointOnSegBy, PointOnSegBz);
      }
    }
  }

  NearestPointX = 0.5f * (PointOnSegAx + PointOnSegBx);
  NearestPointY = 0.5f * (PointOnSegAy + PointOnSegBy);
  NearestPointZ = 0.5f * (PointOnSegAz + PointOnSegBz);

  NearestVectorX = PointOnSegBx - PointOnSegAx;
  NearestVectorY = PointOnSegBy - PointOnSegAy;
  NearestVectorZ = PointOnSegBz - PointOnSegAz;

// optional check to indicate if the lines truly intersect
  true_intersection = (FABS(NearestVectorX) +
                       FABS(NearestVectorY) +
                       FABS(NearestVectorZ)) < epsilon ? true : false;
}

/**************************************************************************
|
|     Method: FindNearestPointOnLineSegment
|
|    Purpose: Given a line (segment) and a point in 3-dimensional space,
|             find the point on the line (segment) that is closest to the
|             point.
|
| Parameters: Input:
|             ------
|             A1x, A1y, A1z   - Coordinates of first defining point of the line/segment
|             Lx, Ly, Lz      - Vector from (A1x, A1y, A1z) to the second defining point
|                               of the line/segment.
|             Bx, By, Bz      - Coordinates of the point
|             infinite_lines  - set to true if lines are to be treated as infinite
|             epsilon_squared - tolerance value to be used to check for degenerate
|                               and parallel lines, and to check for true intersection.
|
|             Output:
|             -------
|             NearestPointX,  - Point on line/segment that is closest to (Bx, By, Bz)
|             NearestPointY,
|             NearestPointZ
|             parameter       - Parametric coordinate of the nearest point along the
|                               line/segment. parameter = 0 at (A1x, A1y, A1z) and
|                               parameter = 1 at the second defining point of the line/
|                               segmetn
**************************************************************************/
void FindNearestPointOnLineSegment(const float A1x, const float A1y, const float A1z,
                                   const float Lx, const float Ly, const float Lz,
                                   const float Bx, const float By, const float Bz,
                                   bool infinite_line, float epsilon_squared, float &NearestPointX,
                                   float &NearestPointY, float &NearestPointZ,
                                   float &parameter)
{
// Line/Segment is degenerate --- special case #1
  float D = Lx * Lx + Ly * Ly + Lz * Lz;
  if (D < epsilon_squared)
  {
    NearestPointX = A1x;
    NearestPointY = A1y;
    NearestPointZ = A1z;
    return;
  }

  float ABx = Bx - A1x;
  float ABy = By - A1y;
  float ABz = Bz - A1z;

// parameter is computed from Equation (20).
  parameter = (Lx * ABx + Ly * ABy + Lz * ABz) / D;

  if (false == infinite_line) parameter = FMAX(0.0f, FMIN(1.0f, parameter));

  NearestPointX = A1x + parameter * Lx;
  NearestPointY = A1y + parameter * Ly;
  NearestPointZ = A1z + parameter * Lz;
  return;
}

/**************************************************************************
|
|     Method: FindNearestPointOfParallelLineSegments
|
|    Purpose: Given two lines (segments) that are known to be parallel, find
|             a representative point on each that is nearest to the other. If
|             the lines are considered to be finite then it is possible that there
|             is one true point on each line that is nearest to the other. This
|             code properly handles this case.
|
|             This is the most difficult line intersection case to handle, since
|             there is potentially a family, or locus of points on each line/segment
|             that are nearest to the other.
| Parameters: Input:
|             ------
|             A1x, A1y, A1z   - Coordinates of first defining point of line/segment A
|             A2x, A2y, A2z   - Coordinates of second defining point of line/segment A
|             Lax, Lay, Laz   - Vector from (A1x, A1y, A1z) to the (A2x, A2y, A2z).
|             B1x, B1y, B1z   - Coordinates of first defining point of line/segment B
|             B2x, B2y, B2z   - Coordinates of second defining point of line/segment B
|             Lbx, Lby, Lbz   - Vector from (B1x, B1y, B1z) to the (B2x, B2y, B2z).
|             infinite_lines  - set to true if lines are to be treated as infinite
|             epsilon_squared - tolerance value to be used to check for degenerate
|                               and parallel lines, and to check for true intersection.
|
|             Output:
|             -------
|             PointOnSegAx,   - Coordinates of the point on segment A that are nearest
|             PointOnSegAy,     to segment B. This corresponds to point C in the text.
|             PointOnSegAz
|             PointOnSegBx,   - Coordinates of the point on segment B that are nearest
|             PointOnSegBy,     to segment A. This corresponds to point D in the text.
|             PointOnSegBz

**************************************************************************/
void FindNearestPointOfParallelLineSegments(float A1x, float A1y, float A1z,
                                            float A2x, float A2y, float A2z,
                                            float Lax, float Lay, float Laz,
                                            float B1x, float B1y, float B1z,
                                            float B2x, float B2y, float B2z,
                                            float Lbx, float Lby, float Lbz,
                                            bool infinite_lines, float epsilon_squared,
                                            float &PointOnSegAx, float &PointOnSegAy, float &PointOnSegAz,
                                            float &PointOnSegBx, float &PointOnSegBy, float &PointOnSegBz)
{
  float s[2], temp;
  FindNearestPointOnLineSegment(A1x, A1y, A1z, Lax, Lay, Laz, B1x, B1y, B1z,
                                true, epsilon_squared, PointOnSegAx, PointOnSegAy, PointOnSegAz, s[0]);
  if (true == infinite_lines)
  {
    PointOnSegBx = B1x;
    PointOnSegBy = B1y;
    PointOnSegBz = B1z;
  }
  else
  {
    float tp[3];
    FindNearestPointOnLineSegment(A1x, A1y, A1z, Lax, Lay, Laz, B2x, B2y, B2z,
                                  true, epsilon_squared, tp[0], tp[1], tp[2], s[1]);
    if (s[0] < 0.f && s[1] < 0.f)
    {
      PointOnSegAx = A1x;
      PointOnSegAy = A1y;
      PointOnSegAz = A1z;
      if (s[0] < s[1])
      {
        PointOnSegBx = B2x;
        PointOnSegBy = B2y;
        PointOnSegBz = B2z;
      }
      else
      {
        PointOnSegBx = B1x;
        PointOnSegBy = B1y;
        PointOnSegBz = B1z;
      }
    }
    else if (s[0] > 1.f && s[1] > 1.f)
    {
      PointOnSegAx = A2x;
      PointOnSegAy = A2y;
      PointOnSegAz = A2z;
      if (s[0] < s[1])
      {
        PointOnSegBx = B1x;
        PointOnSegBy = B1y;
        PointOnSegBz = B1z;
      }
      else
      {
        PointOnSegBx = B2x;
        PointOnSegBy = B2y;
        PointOnSegBz = B2z;
      }
    }
    else
    {
      temp = 0.5f*(FMAX(0.0f, FMIN(1.0f, s[0])) + FMAX(0.0f, FMIN(1.0f, s[1])));
      PointOnSegAx = (A1x + temp * Lax);
      PointOnSegAy = (A1y + temp * Lay);
      PointOnSegAz = (A1z + temp * Laz);
      FindNearestPointOnLineSegment(B1x, B1y, B1z, Lbx, Lby, Lbz,
                                    PointOnSegAx, PointOnSegAy, PointOnSegAz, true,
                                    epsilon_squared, PointOnSegBx, PointOnSegBy, PointOnSegBz, temp);
    }
  }
}

/**************************************************************************
|
|     Method: AdjustNearestPoints
|
|    Purpose: Given nearest point information for two infinite lines, adjust
|             to model finite line segments.
|
| Parameters: Input:
|             ------
|             A1x, A1y, A1z   - Coordinates of first defining point of line/segment A
|             Lax, Lay, Laz   - Vector from (A1x, A1y, A1z) to the (A2x, A2y, A2z).
|             B1x, B1y, B1z   - Coordinates of first defining point of line/segment B
|             Lbx, Lby, Lbz   - Vector from (B1x, B1y, B1z) to the (B2x, B2y, B2z).
|             epsilon_squared - tolerance value to be used to check for degenerate
|                               and parallel lines, and to check for true intersection.
|             s               - parameter representing nearest point on infinite line A
|             t               - parameter representing nearest point on infinite line B
|
|             Output:
|             -------
|             PointOnSegAx,   - Coordinates of the point on segment A that are nearest
|             PointOnSegAy,     to segment B. This corresponds to point C in the text.
|             PointOnSegAz
|             PointOnSegBx,   - Coordinates of the point on segment B that are nearest
|             PointOnSegBy,     to segment A. This corresponds to point D in the text.
|             PointOnSegBz
**************************************************************************/
void AdjustNearestPoints(float A1x, float A1y, float A1z,
                         float Lax, float Lay, float Laz,
                         float B1x, float B1y, float B1z,
                         float Lbx, float Lby, float Lbz,
                         float epsilon_squared, float s, float t,
                         float &PointOnSegAx, float &PointOnSegAy, float &PointOnSegAz,
                         float &PointOnSegBx, float &PointOnSegBy, float &PointOnSegBz)
{
// handle the case where both parameter s and t are out of range
  if (OUT_OF_RANGE(s) && OUT_OF_RANGE(t))
  {
    s = FMAX(0.0f, FMIN(1.0f, s));
    PointOnSegAx = (A1x + s * Lax);
    PointOnSegAy = (A1y + s * Lay);
    PointOnSegAz = (A1z + s * Laz);
    FindNearestPointOnLineSegment(B1x, B1y, B1z, Lbx, Lby, Lbz, PointOnSegAx,
                                  PointOnSegAy, PointOnSegAz, true, epsilon_squared,
                                  PointOnSegBx, PointOnSegBy, PointOnSegBz, t);
    if (OUT_OF_RANGE(t))
    {
      t = FMAX(0.0f, FMIN(1.0f, t));
      PointOnSegBx = (B1x + t * Lbx);
      PointOnSegBy = (B1y + t * Lby);
      PointOnSegBz = (B1z + t * Lbz);
      FindNearestPointOnLineSegment(A1x, A1y, A1z, Lax, Lay, Laz, PointOnSegBx,
                                    PointOnSegBy, PointOnSegBz, false, epsilon_squared,
                                    PointOnSegAx, PointOnSegAy, PointOnSegAz, s);
      FindNearestPointOnLineSegment(B1x, B1y, B1z, Lbx, Lby, Lbz, PointOnSegAx,
                                    PointOnSegAy, PointOnSegAz, false, epsilon_squared,
                                    PointOnSegBx, PointOnSegBy, PointOnSegBz, t);
    }
  }
// otherwise, handle the case where the parameter for only one segment is
// out of range
  else if (OUT_OF_RANGE(s))
  {
    s = FMAX(0.0f, FMIN(1.0f, s));
    PointOnSegAx = (A1x + s * Lax);
    PointOnSegAy = (A1y + s * Lay);
    PointOnSegAz = (A1z + s * Laz);
    FindNearestPointOnLineSegment(B1x, B1y, B1z, Lbx, Lby, Lbz, PointOnSegAx,
                                  PointOnSegAy, PointOnSegAz, false, epsilon_squared,
                                  PointOnSegBx, PointOnSegBy, PointOnSegBz, t);
  }
  else if (OUT_OF_RANGE(t))
  {
    t = FMAX(0.0f, FMIN(1.0f, t));
    PointOnSegBx = (B1x + t * Lbx);
    PointOnSegBy = (B1y + t * Lby);
    PointOnSegBz = (B1z + t * Lbz);
    FindNearestPointOnLineSegment(A1x, A1y, A1z, Lax, Lay, Laz, PointOnSegBx,
                                  PointOnSegBy, PointOnSegBz, false, epsilon_squared,
                                  PointOnSegAx, PointOnSegAy, PointOnSegAz, s);
  }
  else
  {
    assert(0);
  }
}

//-------------------------------------------------------------------------------------------


//  For priority Queue -----------------------------------------------------------------------
#pragma once
#include <exception>
namespace handmade
{
 class Empty : public std::exception
 {
 public:
      const char* what() const
      {
           return "empty Queue";
      }
 };
 static Empty ExceptionEmpty;
 template<typename T> class priorityQueue
 {
 private:
      class Node
      {
      public:
           Node() : next(0)
           {
           }
           ~Node()
           {
           }
           T     data;
           Node *next;
  };
 private:
          Node * topNode;
 public:
          priorityQueue() : topNode(0)
          {
          }
          ~priorityQueue()
          {
          }
 public:
      void push( T data )
      {
           Node *newNode = new Node();
           newNode->data = data;
           Node **cursorNode = &topNode;  
           while( *cursorNode != 0 )
           {   
                if( (*cursorNode)->data < data )
                {
                     newNode->next = *cursorNode;
                     break;
                }
                else
                {
                     cursorNode = &(*cursorNode)->next;
                }
           }
           *cursorNode = newNode;
      }
      void pop()
      {
           if( topNode != 0 )
           {
                Node *removeNode = topNode;
                topNode                 = topNode->next;
                delete removeNode;    
           }
      }
      T top()
      {
           if( topNode == 0 )
           {
                throw ExceptionEmpty;
           }
           return topNode->data;
      }
      bool empty()
      {
           return topNode == 0;   
      }
  };
}

//list for samplings
struct Sample {
	float d_lower;
	int t00;
	int t01;
	
	int t10;
	int t11;
	
	Sample() {}
	Sample(float d, int tt00, int tt01, int tt10, int tt11): d_lower(d), t00(tt00), t01(tt01), t10(tt10), t11(tt11) {}

	friend bool operator < (const Sample& s1, const Sample& s2);
};

bool operator < (const Sample& s1, const Sample& s2)
{
	return s1.d_lower > s2.d_lower;
}

handmade::priorityQueue<Sample> volume_queue;
//---------------------------------------------------------------

//4 subwindows
GLint window, xyzwin, xywin, zxwin, yzwin;
void redisplay_all();



//for distance
float Intersection[] = {0.0, 0.0, 0.0};
float Vector[] = {0.0, 0.0, 0.0};
bool true_intersection = false;;
float ftemp[6];
bool infinite_lines = false;

//Camera rotation
double theta =0;
bool rotate_hit = true;
float speed = 0.02f;

//d^bar
float d_upper;

//points for shortest line
float short_x0, short_y0, short_x1, short_y1; 

//variables for pi
const float pi = 3.141592f;

//matrix for curves
float CurvePts0[1025][3];
float CurvePts1[1025][3];

//matrix for epsilon
float Epsil[2][12];

//initialize window size
GLsizei winWidth = 800, winHeight = 800;
GLsizei sub_width = 400, sub_height = 400;

//control points for a bezier curve
GLfloat ctrlPts0[4][3] = {{-30, 20, 20.}, {-10, 40, 20.}, {10, 0, 20.}, {30, 20, 20.}};  
GLfloat ctrlPts1[4][3] = {{-30, -20, -20.0}, {-10, 0, -20.0}, {10, -40, -20.0}, {30, -20, -20.0}};   

//cheking mouse click
int mouse_hit;
int hit_point(float x,float y, int mode)
{
	//mode 0: xy, mode1: yz, mode2: zx
	switch(mode) {
	case 0:
		for (int i=0; i<4; i++) {
			GLfloat tx = ctrlPts0[i][0] - x;
			GLfloat ty = ctrlPts0[i][1] - y;
			if ((tx*tx + ty*ty) < 3.3) return i;
		}

		for (int i=0; i<4; i++) {
			GLfloat tx = ctrlPts1[i][0] - x;
			GLfloat ty = ctrlPts1[i][1] - y;
			if ((tx*tx + ty*ty) < 3.3) return i+4;
		}
		break;
	case 1:
		for (int i=0; i<4; i++) {
			GLfloat tx = ctrlPts0[i][1] - x;
			GLfloat ty = ctrlPts0[i][2] - y;
			if ((tx*tx + ty*ty) < 3.3) return i;
		}

		for (int i=0; i<4; i++) {
			GLfloat tx = ctrlPts1[i][1] - x;
			GLfloat ty = ctrlPts1[i][2] - y;
			if ((tx*tx + ty*ty) < 3.3) return i+4;
		}
		break;
	case 2:
		for (int i=0; i<4; i++) {
			GLfloat tx = ctrlPts0[i][2] - x;
			GLfloat ty = ctrlPts0[i][0] - y;
			if ((tx*tx + ty*ty) < 3.3) return i;
		}

		for (int i=0; i<4; i++) {
			GLfloat tx = ctrlPts1[i][2] - x;
			GLfloat ty = ctrlPts1[i][0] - y;
			if ((tx*tx + ty*ty) < 3.3) return i+4;
		}
		break;
	}

	return -1;
}


float epsilon0(int h)
{
	return .75f * (1/pow(4.f,h)) * std::max(sqrt(pow((ctrlPts0[0][0] - 2*ctrlPts0[1][0] + ctrlPts0[2][0]),2) + pow((ctrlPts0[0][1] - 2*ctrlPts0[1][1] + ctrlPts0[2][1]),2) + pow((ctrlPts0[0][2] - 2*ctrlPts0[1][2] + ctrlPts0[2][2]),2)) , sqrt(pow((ctrlPts0[1][0] - 2*ctrlPts0[2][0] + ctrlPts0[3][0]),2) + pow((ctrlPts0[1][1] - 2*ctrlPts0[2][1] + ctrlPts0[3][1]),2) + pow((ctrlPts0[1][2] - 2*ctrlPts0[2][2] + ctrlPts0[3][2]),2)));
}

float epsilon1(int h)
{
	return .75f * (1/pow(4.f,h)) * std::max(sqrt(pow((ctrlPts1[0][0] - 2*ctrlPts1[1][0] + ctrlPts1[2][0]),2) + pow((ctrlPts1[0][1] - 2*ctrlPts1[1][1] + ctrlPts1[2][1]),2) + pow((ctrlPts1[0][2] - 2*ctrlPts1[1][2] + ctrlPts1[2][2]),2)) , sqrt(pow((ctrlPts1[1][0] - 2*ctrlPts1[2][0] + ctrlPts1[3][0]),2) + pow((ctrlPts1[1][1] - 2*ctrlPts1[2][1] + ctrlPts1[3][1]),2) + pow((ctrlPts1[1][2] - 2*ctrlPts1[2][2] + ctrlPts1[3][2]),2)));
}


//get C(t)
float cx0(GLfloat t)
{
	return pow((1-t),3) * ctrlPts0[0][0] + 3*pow((1-t),2)*t*ctrlPts0[1][0] + 3*(1-t)*pow(t,2)*ctrlPts0[2][0] + pow(t,3)*ctrlPts0[3][0];
}

float cy0(GLfloat t)
{
	return pow((1-t),3) * ctrlPts0[0][1] + 3*pow((1-t),2)*t*ctrlPts0[1][1] + 3*(1-t)*pow(t,2)*ctrlPts0[2][1] + pow(t,3)*ctrlPts0[3][1];
}

float cz0(GLfloat t)
{
	return pow((1-t),3) * ctrlPts0[0][2] + 3*pow((1-t),2)*t*ctrlPts0[1][2] + 3*(1-t)*pow(t,2)*ctrlPts0[2][2] + pow(t,3)*ctrlPts0[3][2];
}

float cx1(GLfloat t)
{
	return pow((1-t),3) * ctrlPts1[0][0] + 3*pow((1-t),2)*t*ctrlPts1[1][0] + 3*(1-t)*pow(t,2)*ctrlPts1[2][0] + pow(t,3)*ctrlPts1[3][0];
}

float cy1(GLfloat t)
{
	return pow((1-t),3) * ctrlPts1[0][1] + 3*pow((1-t),2)*t*ctrlPts1[1][1] + 3*(1-t)*pow(t,2)*ctrlPts1[2][1] + pow(t,3)*ctrlPts1[3][1];
}

float cz1(GLfloat t)
{
	return pow((1-t),3) * ctrlPts1[0][2] + 3*pow((1-t),2)*t*ctrlPts1[1][2] + 3*(1-t)*pow(t,2)*ctrlPts1[2][2] + pow(t,3)*ctrlPts1[3][2];
}





//----- Part 3 ----------------------------------------------------------------------------------------------------

int mylog(int k)
{
	switch(k) {
	case 1024:
		return 0;
	case 512:
		return 1;
	case 256:
		return 2;
	case 128:
		return 3;
	case 64:
		return 4;
	case 32:
		return 5;
	case 16:
		return 6;
	case 8:
		return 7;
	case 4:
		return 8;
	case 2:
		return 9;
	case 1:
		return 10;
	default:
		return -1;
	}
}


//---- Distance between 2 volumes -------------------
float dist_line(float x00, float y00, float z00, float x01, float y01, float z01, float x10, float y10, float z10 , float x11, float y11, float z11)
{
	IntersectLineSegments(x00, y00, z00, x01, y01, z01,
					 x10, y10, z10, x11, y11, z11,
                    infinite_lines, 1.e-6, ftemp[0], ftemp[1], ftemp[2],
                    ftemp[3], ftemp[4], ftemp[5],
                    Intersection[0], Intersection[1], Intersection[2],
                    Vector[0], Vector[1], Vector[2], true_intersection);
	if (true_intersection)
		return 0.f;
	else {
		float d = sqrt( (ftemp[0] - ftemp[3])*(ftemp[0] - ftemp[3]) + (ftemp[1] - ftemp[4])*(ftemp[1] - ftemp[4]) + (ftemp[2] - ftemp[5])*(ftemp[2] - ftemp[5]));
		return d;
	}
}


float dist_lss(GLfloat x00, GLfloat y00, GLfloat z00, GLfloat x01, GLfloat y01, GLfloat z01, GLfloat r0, GLfloat x10, GLfloat y10, GLfloat z10, GLfloat x11, GLfloat y11, GLfloat z11, GLfloat r1)
{
	float d = dist_line(x00,y00,z00,x01,y01,z01,x10,y10,z10,x11,y11,z11);
	return std::max(0.f, d-r0-r1);
}

float dist_point(float x0, float y0, float z0, float x1, float y1, float z1)
{
	return sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1) + (z0-z1)*(z0-z1));
}


// main part
void shortest_line()
{
	d_upper = std::min(dist_point(CurvePts0[0][0],CurvePts0[0][1],CurvePts0[0][2], CurvePts1[0][0],CurvePts1[0][1],CurvePts1[0][2]), dist_point(CurvePts0[0][0],CurvePts0[0][1],CurvePts0[0][2],CurvePts1[1024][0],CurvePts1[1024][1],CurvePts1[1024][2]));
	d_upper = std::min(d_upper, dist_point(CurvePts0[1024][0],CurvePts0[1024][1],CurvePts0[1024][2], CurvePts1[0][0], CurvePts1[0][1],CurvePts1[0][2]));
	d_upper = std::min(d_upper, dist_point(CurvePts0[1024][0],CurvePts0[1024][1],CurvePts0[1024][2], CurvePts1[1024][0],CurvePts1[1024][1],CurvePts1[1024][2]));

	int t;
	int t00, t01, t10, t11;
	int dt0, dt1;

	while (true) {
		t00 = volume_queue.top().t00;
		t01 = volume_queue.top().t01;
		t10 = volume_queue.top().t10;
		t11 = volume_queue.top().t11;
		dt0 = t01 - t00;
		dt1 = t11 - t10;

		if(dt0 < 2 && dt1 < 2)
			break;
		else {
			if (dt0 >= dt1) {
				t = (t01 + t00) / 2;

				int depth0 = mylog(dt0) + 1;
				int depth1 = mylog(dt1);
	
				float d0_lower = dist_lss(CurvePts0[t00][0], CurvePts0[t00][1], CurvePts0[t00][2], CurvePts0[t][0], CurvePts0[t][1], CurvePts0[t][2], Epsil[0][depth0], CurvePts1[t10][0], CurvePts1[t10][1], CurvePts1[t10][2], CurvePts1[t11][0], CurvePts1[t11][1], CurvePts1[t11][2], Epsil[1][depth1]);
				float d1_lower = dist_lss(CurvePts0[t][0], CurvePts0[t][1], CurvePts0[t][2], CurvePts0[t01][0], CurvePts0[t01][1], CurvePts0[t01][2], Epsil[0][depth0], CurvePts1[t10][0], CurvePts1[t10][1], CurvePts1[t10][2], CurvePts1[t11][0], CurvePts1[t11][1], CurvePts1[t11][2], Epsil[1][depth1]);
			
				d_upper = std::min(d_upper, dist_point(CurvePts0[t][0],CurvePts0[t][1], CurvePts0[t][2], CurvePts1[t10][0],CurvePts1[t10][1], CurvePts1[t10][2]));
				d_upper = std::min(d_upper, dist_point(CurvePts0[t][0],CurvePts0[t][1], CurvePts0[t][2], CurvePts1[t11][0],CurvePts1[t11][1], CurvePts1[t11][2]));

				volume_queue.pop();
				if(d0_lower < d_upper)
					volume_queue.push(Sample(d0_lower, t00, t, t10, t11));
				if(d1_lower < d_upper)
					volume_queue.push(Sample(d1_lower, t, t01, t10, t11));
			} else {
				t = (t10 + t11) / 2;

				int depth0 = mylog(dt0);
				int depth1 = mylog(dt1) + 1;
				
				float d0_lower = dist_lss(CurvePts0[t00][0], CurvePts0[t00][1], CurvePts0[t00][2], CurvePts0[t01][0], CurvePts0[t01][1], CurvePts0[t01][2], Epsil[0][depth0], CurvePts1[t10][0], CurvePts1[t10][1], CurvePts1[t10][2], CurvePts1[t][0], CurvePts1[t][1], CurvePts1[t][2], Epsil[1][depth1]);
				float d1_lower = dist_lss(CurvePts0[t00][0], CurvePts0[t00][1], CurvePts0[t00][2], CurvePts0[t01][0], CurvePts0[t01][1], CurvePts0[t01][2], Epsil[0][depth0], CurvePts1[t][0], CurvePts1[t][1], CurvePts1[t][2], CurvePts1[t11][0], CurvePts1[t11][1], CurvePts1[t11][2], Epsil[1][depth1]);
			
				d_upper = std::min(d_upper, dist_point(CurvePts0[t00][0],CurvePts0[t00][1], CurvePts0[t00][2], CurvePts1[t][0],CurvePts1[t][1], CurvePts1[t][2]));
				d_upper = std::min(d_upper, dist_point(CurvePts0[t01][0],CurvePts0[t01][1], CurvePts0[t01][2], CurvePts1[t][0],CurvePts1[t][1], CurvePts1[t][2]));

				volume_queue.pop();
				if(d0_lower < d_upper)
					volume_queue.push(Sample(d0_lower,t00,t01,t10,t));
				if(d1_lower < d_upper)
					volume_queue.push(Sample(d1_lower,t00,t01,t,t11));
			}
		}
	}

	t00 = volume_queue.top().t00;
	t01 = volume_queue.top().t01;
	
	t10 = volume_queue.top().t10;
	t11 = volume_queue.top().t11;

	IntersectLineSegments(CurvePts0[t00][0], CurvePts0[t00][1], CurvePts0[t00][2], CurvePts0[t01][0], CurvePts0[t01][1], CurvePts0[t01][2],
					 CurvePts1[t10][0], CurvePts1[t10][1], CurvePts1[t10][2], CurvePts1[t11][0], CurvePts1[t11][1], CurvePts1[t11][2],
                    infinite_lines, 1.e-6, ftemp[0], ftemp[1], ftemp[2],
                    ftemp[3], ftemp[4], ftemp[5],
                    Intersection[0], Intersection[1], Intersection[2],
                    Vector[0], Vector[1], Vector[2], true_intersection);

	while(!volume_queue.empty()) 
		volume_queue.pop();
}

//------------------------------------------------------------------------------
void init()
{
	//For part 3
	for (int i = 0; i<=1024; i++) {
		CurvePts0[i][0] = cx0(float(i) /1024.f);
		CurvePts0[i][1] = cy0(float(i) /1024.f);
		CurvePts0[i][2] = cz0(float(i) /1024.f);


		CurvePts1[i][0] = cx1(float(i) /1024.f);
		CurvePts1[i][1] = cy1(float(i) /1024.f);
		CurvePts1[i][2] = cz1(float(i) /1024.f);
	}

	for (int j = 0; j<12; j++) {
		Epsil[0][j] = epsilon0(j);
		Epsil[1][j] = epsilon1(j);
	}

	volume_queue.push(Sample(0.f,0,1024,0,1024));

	shortest_line();
}

void xyz_display() {	
	//Curve0
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMap1f(GL_MAP1_VERTEX_3, 0.0, 1.0, 3, 4, *ctrlPts0);
	glEnable(GL_MAP1_VERTEX_3);

	GLint k;

	glColor3f(0.0, 0.0, 0.0);
	glBegin (GL_LINE_STRIP);
		for (k = 0; k<=50; k++)
			glEvalCoord1f(GLfloat(k)/50.f);
	glEnd();


	// Curve1
	glMap1f(GL_MAP1_VERTEX_3, 0.0, 1.0, 3, 4, *ctrlPts1);
	glEnable(GL_MAP1_VERTEX_3);

	glColor3f(0.0, 0.0, 0.0);
	glBegin (GL_LINE_STRIP);
		for (k = 0; k<=50; k++)
			glEvalCoord1f(GLfloat(k)/50.f);
	glEnd();

	//coordinates
	glBegin(GL_LINES);
		glColor3f (0.0, 1.0, 0.0);
		glVertex3f(0.f,0.f, 0.f);
		glVertex3f(0.f,33.f, 0.f);

		glColor3f (1.0, 0.0, 0.0);
		glVertex3f(0.f,0.f, 0.f);
		glVertex3f(33.f,0.f, 0.f);
		
		glColor3f (0.0, 0.0, 1.0);
		glVertex3f(0.f,0.f, 0.f);
		glVertex3f(0.f,0.f, 33.f);
	glEnd();


	// rotate a camera
	if (theta >= 360)
		theta = 0;
	if (rotate_hit)
		theta += speed;

	double x = 100*cos(theta);
	double z = 100*sin(theta);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(x, 100, z, 0, 0., 0., 0., 1., 0.);

	// shortest line
	glColor3f(0.0, 1.0, 0.0);
	glBegin (GL_LINES);
		glVertex3f(ftemp[0],ftemp[1],ftemp[2]);
		glVertex3f(ftemp[3],ftemp[4],ftemp[5]);
	glEnd();


	glutPostRedisplay();

	glutSwapBuffers();
	glFlush();
}

void xyz_reshape(GLint subWidth, GLint subHeight)
{
	//set the color of display window
	glClearColor (1.f,1.f,1.f, 0.f);

	glViewport(0, 0, subWidth, subHeight);

	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();

	glFrustum(-subWidth * 0.06 ,subWidth * 0.06, -subHeight * 0.06 ,subHeight * 0.06, 50, 500.);

}

void zx_display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// coordinate
	glBegin(GL_LINES);
		glColor3f (0.0, 1.0, 0.0);
		glVertex3f(0.f,0.f, 0.f);
		glVertex3f(0.f,33.f, 0.f);

		glColor3f (1.0, 0.0, 0.0);
		glVertex3f(0.f,0.f, 0.f);
		glVertex3f(33.f,0.f, 0.f);
		
		glColor3f (0.0, 0.0, 1.0);
		glVertex3f(0.f,0.f, 0.f);
		glVertex3f(0.f,0.f, 33.f);
	glEnd();

		//Curve0
	glMap1f(GL_MAP1_VERTEX_3, 0.0, 1.0, 3, 4, *ctrlPts0);
	glEnable(GL_MAP1_VERTEX_3);

	GLint k;

	glColor3f(0.0, 0.0, 0.0);
	glBegin (GL_LINE_STRIP);
		for (k = 0; k<=50; k++)
			glEvalCoord1f(GLfloat(k)/50.f);
	glEnd();

	glColor3f (1.0, 0.0, 0.0);
	glPointSize(7.0);
	glBegin(GL_POINTS);
		for(k=0; k<4; k++) 
			glVertex3fv(&ctrlPts0[k][0]);
	glEnd();

	// Curve1
	glMap1f(GL_MAP1_VERTEX_3, 0.0, 1.0, 3, 4, *ctrlPts1);
	glEnable(GL_MAP1_VERTEX_3);

	glColor3f(0.0, 0.0, 0.0);
	glBegin (GL_LINE_STRIP);
		for (k = 0; k<=50; k++)
			glEvalCoord1f(GLfloat(k)/50.f);
	glEnd();

	glColor3f (1.0, 0.5, 0.0);
	glPointSize(7.0);
	glBegin(GL_POINTS);
		for(k=0; k<4; k++) 
			glVertex3fv(&ctrlPts1[k][0]);
	glEnd();

	glColor3f(0.0, 1.0, 0.0);
	glBegin (GL_LINES);
		glVertex3f(ftemp[0],ftemp[1],ftemp[2]);
		glVertex3f(ftemp[3],ftemp[4],ftemp[5]);
	glEnd();

	glutSwapBuffers();
}

void zx_reshape(int subWidth, int subHeight)
{
	//set the color of display window
	glClearColor (224.f/255.f, 1, 224.f/255.f, 0.f);

	glViewport(0, 0, subWidth, subHeight);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

	glOrtho(-subWidth * 0.15 ,subWidth * 0.15, -subHeight * 0.15 ,subHeight * 0.15, -1000, 1000);
    glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0, 0, 0, 0, -1, 0, 0, 0, -1);
}

void zx_mouse(int button, int action, int x, int z)
{
	z = sub_height - z;

	GLfloat new_x = (x - (GLfloat)sub_width/2.f) * 0.3f;
	GLfloat new_z = ((GLfloat)sub_height/2.f - z) * 0.3f;
	
	if (button == GLUT_LEFT_BUTTON && action == GLUT_DOWN)		//마우스 왼쪽 키를 누르는 이벤트가 발생할 경우
		mouse_hit = hit_point(new_z,new_x,2);
	else if (button == GLUT_RIGHT_BUTTON) {
		mouse_hit = -1;

		ctrlPts0[0][0] = -30;
		ctrlPts0[0][1] = 20;
		ctrlPts0[0][2] = 20;
		ctrlPts0[1][0] = -10;
		ctrlPts0[1][1] = 40;
		ctrlPts0[1][2] = 20;
		ctrlPts0[2][0] = 10;
		ctrlPts0[2][1] = 0;
		ctrlPts0[2][2] = 20;
		ctrlPts0[3][0] = 30;
		ctrlPts0[3][1] = 20.0;
		ctrlPts0[3][2] = 20.0;
		
		ctrlPts1[0][0] = -30.0;
		ctrlPts1[0][1] = -20.0;
		ctrlPts1[0][2] = -20.0;
		ctrlPts1[1][0] = -10.0;
		ctrlPts1[1][1] = 0.0;
		ctrlPts1[1][2] = -20.0;
		ctrlPts1[2][0] = 10.0;
		ctrlPts1[2][1] = -40.0;
		ctrlPts1[2][2] = -20.0;
		ctrlPts1[3][0] = 30.0;
		ctrlPts1[3][1] = -20.0;
		ctrlPts1[3][2] = -20.0;

		//For part 3
		for (int i = 0; i<=1024; i++) {
			CurvePts0[i][0] = cx0(float(i) /1024.f);
			CurvePts0[i][1] = cy0(float(i) /1024.f);
			CurvePts0[i][2] = cz0(float(i) /1024.f);


			CurvePts1[i][0] = cx1(float(i) /1024.f);
			CurvePts1[i][1] = cy1(float(i) /1024.f);
			CurvePts1[i][2] = cz1(float(i) /1024.f);
		}

		for (int j = 0; j<12; j++) {
			Epsil[0][j] = epsilon0(j);
			Epsil[1][j] = epsilon1(j);
		}

		volume_queue.push(Sample(0.f,0,1024,0,1024));

		shortest_line();
	}
	else mouse_hit = -1;
	redisplay_all();
}

void zx_motion(int x, int z) {
	z = sub_height - z;

	GLfloat new_x = (x - (GLfloat)sub_width/2.f) * 0.3f;
	GLfloat new_z = ((GLfloat)sub_height/2.f - z) * 0.3f;

	if (mouse_hit != -1) {
		if(mouse_hit < 4) {
			ctrlPts0[mouse_hit][2] = new_z;
			ctrlPts0[mouse_hit][0] = new_x;
		} else {
			ctrlPts1[mouse_hit-4][2] = new_z;
			ctrlPts1[mouse_hit-4][0] = new_x;
		}
	}

	//For part 3
	for (int i = 0; i<=1024; i++) {
		CurvePts0[i][0] = cx0(float(i) /1024.f);
		CurvePts0[i][1] = cy0(float(i) /1024.f);
		CurvePts0[i][2] = cz0(float(i) /1024.f);


		CurvePts1[i][0] = cx1(float(i) /1024.f);
		CurvePts1[i][1] = cy1(float(i) /1024.f);
		CurvePts1[i][2] = cz1(float(i) /1024.f);
	}

	for (int j = 0; j<12; j++) {
		Epsil[0][j] = epsilon0(j);
		Epsil[1][j] = epsilon1(j);
	}

	volume_queue.push(Sample(0.f,0,1024,0,1024));

	shortest_line();

	redisplay_all();
}



void xy_display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	//coordinate
	glBegin(GL_LINES);
		glColor3f (0.0, 1.0, 0.0);
		glVertex3f(0.f,0.f, 0.f);
		glVertex3f(0.f,33.f, 0.f);

		glColor3f (1.0, 0.0, 0.0);
		glVertex3f(0.f,0.f, 0.f);
		glVertex3f(33.f,0.f, 0.f);
		
		glColor3f (0.0, 0.0, 1.0);
		glVertex3f(0.f,0.f, 0.f);
		glVertex3f(0.f,0.f, 33.f);
	glEnd();

	
	//Curve0
	glMap1f(GL_MAP1_VERTEX_3, 0.0, 1.0, 3, 4, *ctrlPts0);
	glEnable(GL_MAP1_VERTEX_3);

	GLint k;

	glColor3f(0.0, 0.0, 0.0);
	glBegin (GL_LINE_STRIP);
		for (k = 0; k<=50; k++)
			glEvalCoord1f(GLfloat(k)/50.f);
	glEnd();

	glColor3f (1.0, 0.0, 0.0);
	glPointSize(7.0);
	glBegin(GL_POINTS);
		for(k=0; k<4; k++) 
			glVertex3fv(&ctrlPts0[k][0]);
	glEnd();

	// Curve1
	glMap1f(GL_MAP1_VERTEX_3, 0.0, 1.0, 3, 4, *ctrlPts1);
	glEnable(GL_MAP1_VERTEX_3);

	glColor3f(0.0, 0.0, 0.0);
	glBegin (GL_LINE_STRIP);
		for (k = 0; k<=50; k++)
			glEvalCoord1f(GLfloat(k)/50.f);
	glEnd();

	glColor3f (1.0, 0.5, 0.0);
	glPointSize(7.0);
	glBegin(GL_POINTS);
		for(k=0; k<4; k++) 
			glVertex3fv(&ctrlPts1[k][0]);
	glEnd();

	glColor3f(0.0, 1.0, 0.0);
	glBegin (GL_LINES);
		glVertex3f(ftemp[0],ftemp[1],ftemp[2]);
		glVertex3f(ftemp[3],ftemp[4],ftemp[5]);
	glEnd();

	 glutSwapBuffers();
}

void xy_reshape(int subWidth, int subHeight)
{
	//set the color of display window
	glClearColor (224.f/255.f, 1, 1, 0.f);

	glViewport(0, 0, subWidth, subHeight);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    glOrtho(-subWidth * 0.15 ,subWidth * 0.15, -subHeight * 0.15 ,subHeight * 0.15, -1000, 1000);
    glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void xy_mouse(int button, int action, int x, int y)
{
	y = sub_height - y;
	
	GLfloat new_x = (x - (GLfloat)sub_width/2.f) * 0.3f;
	GLfloat new_y = (y - (GLfloat)sub_height/2.f) * 0.3f;

	if (button == GLUT_LEFT_BUTTON && action == GLUT_DOWN)		//마우스 왼쪽 키를 누르는 이벤트가 발생할 경우
		mouse_hit = hit_point(new_x,new_y,0);
	else if (button == GLUT_RIGHT_BUTTON) {
		mouse_hit = -1;

		ctrlPts0[0][0] = -30;
		ctrlPts0[0][1] = 20;
		ctrlPts0[0][2] = 20;
		ctrlPts0[1][0] = -10;
		ctrlPts0[1][1] = 40;
		ctrlPts0[1][2] = 20;
		ctrlPts0[2][0] = 10;
		ctrlPts0[2][1] = 0;
		ctrlPts0[2][2] = 20;
		ctrlPts0[3][0] = 30;
		ctrlPts0[3][1] = 20.0;
		ctrlPts0[3][2] = 20.0;
		
		ctrlPts1[0][0] = -30.0;
		ctrlPts1[0][1] = -20.0;
		ctrlPts1[0][2] = -20.0;
		ctrlPts1[1][0] = -10.0;
		ctrlPts1[1][1] = 0.0;
		ctrlPts1[1][2] = -20.0;
		ctrlPts1[2][0] = 10.0;
		ctrlPts1[2][1] = -40.0;
		ctrlPts1[2][2] = -20.0;
		ctrlPts1[3][0] = 30.0;
		ctrlPts1[3][1] = -20.0;
		ctrlPts1[3][2] = -20.0;

		//For part 3
		for (int i = 0; i<=1024; i++) {
			CurvePts0[i][0] = cx0(float(i) /1024.f);
			CurvePts0[i][1] = cy0(float(i) /1024.f);
			CurvePts0[i][2] = cz0(float(i) /1024.f);


			CurvePts1[i][0] = cx1(float(i) /1024.f);
			CurvePts1[i][1] = cy1(float(i) /1024.f);
			CurvePts1[i][2] = cz1(float(i) /1024.f);
		}

		for (int j = 0; j<12; j++) {
			Epsil[0][j] = epsilon0(j);
			Epsil[1][j] = epsilon1(j);
		}

		volume_queue.push(Sample(0.f,0,1024,0,1024));

		shortest_line();
	}
	else mouse_hit = -1;
	redisplay_all();
}

void xy_motion(int x, int y) {
	y = sub_height - y;
	GLfloat new_x = (x - (GLfloat)sub_width/2.f) * 0.3f;
	GLfloat new_y = (y - (GLfloat)sub_height/2.f) * 0.3f;

	if (mouse_hit != -1) {
		if(mouse_hit < 4) {
			ctrlPts0[mouse_hit][0] = new_x;
			ctrlPts0[mouse_hit][1] = new_y;
		} else {
			ctrlPts1[mouse_hit-4][0] = new_x;
			ctrlPts1[mouse_hit-4][1] = new_y;
		}
	}

	//For part 3
	for (int i = 0; i<=1024; i++) {
		CurvePts0[i][0] = cx0(float(i) /1024.f);
		CurvePts0[i][1] = cy0(float(i) /1024.f);
		CurvePts0[i][2] = cz0(float(i) /1024.f);


		CurvePts1[i][0] = cx1(float(i) /1024.f);
		CurvePts1[i][1] = cy1(float(i) /1024.f);
		CurvePts1[i][2] = cz1(float(i) /1024.f);
	}

	for (int j = 0; j<12; j++) {
		Epsil[0][j] = epsilon0(j);
		Epsil[1][j] = epsilon1(j);
	}

	volume_queue.push(Sample(0.f,0,1024,0,1024));

	shortest_line();

	redisplay_all();
}


void yz_display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glBegin(GL_LINES);
		glColor3f (0.0, 1.0, 0.0);
		glVertex3f(0.f,0.f, 0.f);
		glVertex3f(0.f,33.f, 0.f);

		glColor3f (1.0, 0.0, 0.0);
		glVertex3f(0.f,0.f, 0.f);
		glVertex3f(33.f,0.f, 0.f);
		
		glColor3f (0.0, 0.0, 1.0);
		glVertex3f(0.f,0.f, 0.f);
		glVertex3f(0.f,0.f, 33.f);
	glEnd();

		//Curve0
	glMap1f(GL_MAP1_VERTEX_3, 0.0, 1.0, 3, 4, *ctrlPts0);
	glEnable(GL_MAP1_VERTEX_3);

	GLint k;

	glColor3f(0.0, 0.0, 0.0);
	glBegin (GL_LINE_STRIP);
		for (k = 0; k<=50; k++)
			glEvalCoord1f(GLfloat(k)/50.f);
	glEnd();

	glColor3f (1.0, 0.0, 0.0);
	glPointSize(7.0);
	glBegin(GL_POINTS);
		for(k=0; k<4; k++) 
			glVertex3fv(&ctrlPts0[k][0]);
	glEnd();

	// Curve1
	glMap1f(GL_MAP1_VERTEX_3, 0.0, 1.0, 3, 4, *ctrlPts1);
	glEnable(GL_MAP1_VERTEX_3);

	glColor3f(0.0, 0.0, 0.0);
	glBegin (GL_LINE_STRIP);
		for (k = 0; k<=50; k++)
			glEvalCoord1f(GLfloat(k)/50.f);
	glEnd();

	glColor3f (1.0, 0.5, 0.0);
	glPointSize(7.0);
	glBegin(GL_POINTS);
		for(k=0; k<4; k++) 
			glVertex3fv(&ctrlPts1[k][0]);
	glEnd();

	glColor3f(0.0, 1.0, 0.0);
	glBegin (GL_LINES);
		glVertex3f(ftemp[0],ftemp[1],ftemp[2]);
		glVertex3f(ftemp[3],ftemp[4],ftemp[5]);
	glEnd();

	 glutSwapBuffers();
}

void yz_reshape(int subWidth, int subHeight)
{
	//set the color of display window
		glClearColor (1, 224.f/255.f, 224.f/255.f, 0.f);

	glViewport(0, 0, subWidth, subHeight);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    glOrtho(-subWidth * 0.15 ,subWidth * 0.15, -subHeight * 0.15 ,subHeight * 0.15, -1000, 1000);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
	gluLookAt(0, 0, 0, -1, 0, 0, 0, 1, 0);
}


void yz_mouse(int button, int action, int z, int y)
{
	y = sub_height - y;

	GLfloat new_z = ((GLfloat)sub_width/2.f - z) * 0.3f;
	GLfloat new_y = (y - (GLfloat)sub_height/2.f) * 0.3f;

	if (button == GLUT_LEFT_BUTTON && action == GLUT_DOWN)		//마우스 왼쪽 키를 누르는 이벤트가 발생할 경우
		mouse_hit = hit_point(new_y,new_z,1);
	else if (button == GLUT_RIGHT_BUTTON) {
		mouse_hit = -1;

		ctrlPts0[0][0] = -30;
		ctrlPts0[0][1] = 20;
		ctrlPts0[0][2] = 20;
		ctrlPts0[1][0] = -10;
		ctrlPts0[1][1] = 40;
		ctrlPts0[1][2] = 20;
		ctrlPts0[2][0] = 10;
		ctrlPts0[2][1] = 0;
		ctrlPts0[2][2] = 20;
		ctrlPts0[3][0] = 30;
		ctrlPts0[3][1] = 20.0;
		ctrlPts0[3][2] = 20.0;
		
		ctrlPts1[0][0] = -30.0;
		ctrlPts1[0][1] = -20.0;
		ctrlPts1[0][2] = -20.0;
		ctrlPts1[1][0] = -10.0;
		ctrlPts1[1][1] = 0.0;
		ctrlPts1[1][2] = -20.0;
		ctrlPts1[2][0] = 10.0;
		ctrlPts1[2][1] = -40.0;
		ctrlPts1[2][2] = -20.0;
		ctrlPts1[3][0] = 30.0;
		ctrlPts1[3][1] = -20.0;
		ctrlPts1[3][2] = -20.0;

		//For part 3
		for (int i = 0; i<=1024; i++) {
			CurvePts0[i][0] = cx0(float(i) /1024.f);
			CurvePts0[i][1] = cy0(float(i) /1024.f);
			CurvePts0[i][2] = cz0(float(i) /1024.f);


			CurvePts1[i][0] = cx1(float(i) /1024.f);
			CurvePts1[i][1] = cy1(float(i) /1024.f);
			CurvePts1[i][2] = cz1(float(i) /1024.f);
		}

		for (int j = 0; j<12; j++) {
			Epsil[0][j] = epsilon0(j);
			Epsil[1][j] = epsilon1(j);
		}

		volume_queue.push(Sample(0.f,0,1024,0,1024));

		shortest_line();
	}
	else mouse_hit = -1;
	redisplay_all();
}

void yz_motion(int z, int y) {
	y = sub_height - y;

	GLfloat new_z = ((GLfloat)sub_width/2.f - z) * 0.3f;
	GLfloat new_y = (y - (GLfloat)sub_height/2.f) * 0.3f;

	if (mouse_hit != -1) {
		if(mouse_hit < 4) {
			ctrlPts0[mouse_hit][1] = new_y;
			ctrlPts0[mouse_hit][2] = new_z;
		} else {
			ctrlPts1[mouse_hit-4][1] = new_y;
			ctrlPts1[mouse_hit-4][2] = new_z;
		}
	}

	//For part 3
	for (int i = 0; i<=1024; i++) {
		CurvePts0[i][0] = cx0(float(i) /1024.f);
		CurvePts0[i][1] = cy0(float(i) /1024.f);
		CurvePts0[i][2] = cz0(float(i) /1024.f);


		CurvePts1[i][0] = cx1(float(i) /1024.f);
		CurvePts1[i][1] = cy1(float(i) /1024.f);
		CurvePts1[i][2] = cz1(float(i) /1024.f);
	}

	for (int j = 0; j<12; j++) {
		Epsil[0][j] = epsilon0(j);
		Epsil[1][j] = epsilon1(j);
	}

	volume_queue.push(Sample(0.f,0,1024,0,1024));

	shortest_line();

	redisplay_all();
}



void main_reshape(int newWidth,  int newHeight) 
{
	
	winWidth = newWidth;
	winHeight = newHeight;

    glViewport(0, 0, winWidth, winHeight);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, winWidth, winHeight, 0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    sub_width = winWidth/2;
    sub_height = winHeight/2;
    
	glutSetWindow(xyzwin);
	glutPositionWindow(sub_width, 0);
    glutReshapeWindow(sub_width, sub_height);

	glutSetWindow(zxwin);
	glutPositionWindow(0, 0);
    glutReshapeWindow(sub_width, sub_height);

	glutSetWindow(xywin);
	glutPositionWindow(0,sub_height);
    glutReshapeWindow(sub_width, sub_height);

	glutSetWindow(yzwin);
	glutPositionWindow(sub_width, sub_height);
    glutReshapeWindow(sub_width, sub_height);
}

void main_display(void)
{
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glutSwapBuffers();
}

void keyboard(unsigned char key, int x, int y) {
	switch(key) {
	case '1':
		speed = 0.005f;
		break;
	case '2':
		speed = 0.01f;
		break;
	case '3':
		speed = 0.02f;
		break;
	case '4':
		speed = 0.03f;
		break;
	case '5':
		speed = 0.06f;
		break;
	case 32:
		rotate_hit ^= true;
	default: break;
	}	
}

void mouse(int button, int action, int x, int y)
{
	if (button == GLUT_RIGHT_BUTTON) {
		ctrlPts0[0][0] = -30;
		ctrlPts0[0][1] = 20;
		ctrlPts0[0][2] = 20;
		ctrlPts0[1][0] = -10;
		ctrlPts0[1][1] = 40;
		ctrlPts0[1][2] = 20;
		ctrlPts0[2][0] = 10;
		ctrlPts0[2][1] = 0;
		ctrlPts0[2][2] = 20;
		ctrlPts0[3][0] = 30;
		ctrlPts0[3][1] = 20.0;
		ctrlPts0[3][2] = 20.0;
		
		ctrlPts1[0][0] = -30.0;
		ctrlPts1[0][1] = -20.0;
		ctrlPts1[0][2] = -20.0;
		ctrlPts1[1][0] = -10.0;
		ctrlPts1[1][1] = 0.0;
		ctrlPts1[1][2] = -20.0;
		ctrlPts1[2][0] = 10.0;
		ctrlPts1[2][1] = -40.0;
		ctrlPts1[2][2] = -20.0;
		ctrlPts1[3][0] = 30.0;
		ctrlPts1[3][1] = -20.0;
		ctrlPts1[3][2] = -20.0;

		//For part 3
		for (int i = 0; i<=1024; i++) {
			CurvePts0[i][0] = cx0(float(i) /1024.f);
			CurvePts0[i][1] = cy0(float(i) /1024.f);
			CurvePts0[i][2] = cz0(float(i) /1024.f);


			CurvePts1[i][0] = cx1(float(i) /1024.f);
			CurvePts1[i][1] = cy1(float(i) /1024.f);
			CurvePts1[i][2] = cz1(float(i) /1024.f);
		}

		for (int j = 0; j<12; j++) {
			Epsil[0][j] = epsilon0(j);
			Epsil[1][j] = epsilon1(j);
		}

		volume_queue.push(Sample(0.f,0,1024,0,1024));

		shortest_line();
	}
	redisplay_all();
}



void redisplay_all(void)
{
	glutSetWindow(xywin);
	xy_reshape(sub_width, sub_height);
    glutPostRedisplay();
	glutSetWindow(yzwin);
	yz_reshape(sub_width, sub_height);
    glutPostRedisplay();
	glutSetWindow(zxwin);
	zx_reshape(sub_width,sub_height);
	glutPostRedisplay();
}



int main(int argc, char** argv) 
{
	glutInit(&argc,argv);
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize (winWidth, winHeight);
	glutInitWindowPosition(50, 50);

	init();

	window = glutCreateWindow("Bezier Curve in 3D");
	glutReshapeFunc(main_reshape);
    glutDisplayFunc(main_display);
	glutKeyboardFunc(keyboard);
	glutMouseFunc(mouse);

	xyzwin = glutCreateSubWindow(window, winWidth/2 , 0, winWidth/2, winHeight/2);
	glutReshapeFunc(xyz_reshape);
	glutDisplayFunc(xyz_display);
	glutKeyboardFunc(keyboard);
	glutMouseFunc(mouse);

	zxwin = glutCreateSubWindow(window, 0, 0, winWidth/2, winHeight/2);
	glutReshapeFunc(zx_reshape);
	glutDisplayFunc(zx_display);
	glutKeyboardFunc(keyboard);
	glutMouseFunc(zx_mouse);
	glutMotionFunc(zx_motion);

	xywin = glutCreateSubWindow(window, 0, winHeight/2, winWidth/2, winHeight/2);
	glutReshapeFunc(xy_reshape);
	glutDisplayFunc(xy_display);
	glutKeyboardFunc(keyboard);
	glutMouseFunc(xy_mouse);
	glutMotionFunc(xy_motion);

	yzwin = glutCreateSubWindow(window, winWidth/2 ,winHeight/2 , winWidth/2, winHeight/2);
	glutReshapeFunc(yz_reshape);
	glutDisplayFunc(yz_display);
	glutKeyboardFunc(keyboard);
	glutMouseFunc(yz_mouse);
	glutMotionFunc(yz_motion);

	glutMainLoop();
	return 0;
}
