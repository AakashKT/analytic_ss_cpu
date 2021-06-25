#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_GHPolygon_H
#define PBRT_INTEGRATORS_GHPolygon_H

// integrators/path.h*
#include "pbrt.h"
#include "core/geometry.h"
#include "integrators/polygonutils.h"

namespace PolyClip {

class Polygon
{
public:
	Vertex * start_;
	Vertex * end_;

	typedef PolygonVertexIterator iterator;
	int vertexNum_;

	Polygon(std::vector<Point2d> vertices);
	~Polygon();

	/* IMPORTANT:
	* Copy constructor and assignment constructor. All of them are depth copy implementation.
	* Note: not all of vertices properties cannot be copied because some properties will be modified by polygon clipping operation.
	* So we just implement to copy basic connection relationship and vertex position property.
	* We also kick out intersection points which were added by polygon clipping operation.
	* If you use these function to copy, REMEMBER do polygon clipping operation from scratch.
	*/
	Polygon(Polygon& poly);
	Polygon& operator=(Polygon& poly);


	// insert vertex to pos
	void Insert(Vertex* new_v, Vertex * pos);
	void Insert(const Point2d& vertex, Vertex * pos);
	void Insert(float x, float y, Vertex* pos);
	// delete vertex on specific pos
	void Remove(Vertex* pos);

	// generate iterator
	iterator begin()
	{
		return PolygonVertexIterator(start_, 0);
	}

	// generate end iterator for visiting all of vertices.
	iterator end()
	{
		return PolygonVertexIterator(end_, vertexNum_);
	}
	// generate end iterator. Note: let iterator more visit 1 element to use start point as end point. 
	iterator endEdge()
	{
		return PolygonVertexIterator(end_, vertexNum_ + 1);
	}

protected:
	friend class PolygonVertexIterator;
};

}

#endif