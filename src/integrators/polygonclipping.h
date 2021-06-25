#ifndef PBRT_INTEGRATORS_GHPolygonClipping_H
#define PBRT_INTEGRATORS_GHPolygonClipping_H

#include "polygon.h"

namespace PolyClip
{
	enum {MarkIntersection, MarkUnion, MarkDifferentiate};
	class PloygonOpration
	{
	public:
		/************************* Phase 1 **************************/
		// Detect intersection points between two polygons.  
		static void DetectIntersection(Polygon& clip_poly, Polygon& sub_poly);


		/************************* Phase 2 **************************/
		// Mark points for intersection
		// true: normally intersection
		// false: no intersection or one polygon is inside the other.
		static bool Mark(Polygon& clip_poly, Polygon& sub_poly, /*out parameter*/ std::vector<std::vector<Point2d>>& possible_result, int markType,
							bool *didIntersect, int *indicator);


		/************************* Phase 3 **************************/
		// extract result polygons (intersection, union or differentiate)
		static std::vector<std::vector<Point2d>> ExtractIntersectionResults(Polygon& clip_poly);
		static std::vector<std::vector<Point2d>> ExtractUnionResults(Polygon& clip_poly);
		static std::vector<std::vector<Point2d>> ExtractDifferentiateResults(Polygon& clip_poly);

		// print
		static  void print(Polygon& polygon)
		{
			for (auto iter = polygon.begin(); iter != polygon.end(); ++iter)
				std::cout << (*iter).x_ << " " << (*iter).y_ << "\n";
		}

	private:
		// utils function.
		// detect intersection proportion between two line segments.
		// Here we use reference, so that modify original data in case boundary condition
		static bool LineSegmentIntersection(Point2d& p1, Point2d& p2,
			Point2d& q1, Point2d& q2,
			float& alphaP, float& alphaQ);

		static Vertex CreateVertex(const Point2d& p1, const Point2d& p2, float alpha);

		// a point is inside a polygon?
		// See https://wrf.ecse.rpi.edu//Research/Short_Notes/pnpoly.html
		// true: inside polygon
		// false: outside
		static bool Pnpoly(const Point2d& p, Polygon& polygon);

	
	};


}

#endif