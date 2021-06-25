#ifndef PBRT_INTEGRATORS_GHPolygonUtils_H
#define PBRT_INTEGRATORS_GHPolygonUtils_H

#include "pbrt.h"

/**
 * All of necessary data structures are defined here
 */


namespace PolyClip
{
	struct Point2d
	{
		float x_, y_, z_;
		bool hasZCoord;

		Point2d():x_(0.0),y_(0.0), hasZCoord(false) {}
		Point2d(float x,float y):x_(x),y_(y), hasZCoord(false) {}

		void setZCoord(float z) {
			hasZCoord = true;
			z_ = z;
		}
	};


	struct Vertex
	{
		float x_, y_, z_;
		bool hasZCoord;
		Vertex * next_;
		Vertex * prev_;
		bool intersect_;
		bool entryExit_;
		Vertex * neighbour_;
		float alpha_;
		bool processd;

		Vertex()
			:x_(0.0),y_(0.0), alpha_(0.0),next_(nullptr), prev_(nullptr),neighbour_(nullptr), intersect_(false), entryExit_(false),processd(false), hasZCoord(false)
		{}
		Vertex(const Point2d& p)
			:x_(p.x_),y_(p.y_), alpha_(0.0), next_(nullptr), prev_(nullptr), neighbour_(nullptr), intersect_(false), entryExit_(false), processd(false), hasZCoord(false)
		{}
		Vertex(float x, float y)
			:x_(x), y_(y), alpha_(0.0), next_(nullptr), prev_(nullptr), neighbour_(nullptr), intersect_(false), entryExit_(false), processd(false), hasZCoord(false)
		{}
	};

	// This object is used for sort insection vertices so that we can insert these at correct order
	struct VertexPtrDistance
	{
		Vertex * ptr;
		float distance;
		VertexPtrDistance(Vertex* vPtr, float dis):ptr(vPtr),distance(dis){}
	};
	struct SortVertexPtrDistance
	{
		bool operator()(const VertexPtrDistance& v1, const VertexPtrDistance& v2)
		{
			return v1.distance < v2.distance;
		}
	};



	class PolygonVertexIterator
	{
	public:
		PolygonVertexIterator(Vertex * cv, int visitedCount) :currentVertex_(cv), visitedCount_(visitedCount) {};

		PolygonVertexIterator operator++()
		{
			visitedCount_++;
			currentVertex_ = currentVertex_->next_;
			return *this;
		}
		PolygonVertexIterator operator--()
		{
			visitedCount_--;
			currentVertex_ = currentVertex_->prev_;
			return *this;
		}
		Vertex& operator*()
		{
			return *currentVertex_;
		}

		bool operator!=(const PolygonVertexIterator& p)
		{
			return std::abs(visitedCount_) != std::abs(p.visitedCount_);
		}
		Vertex* eval() { return currentVertex_; }

		PolygonVertexIterator next()
		{
			PolygonVertexIterator nextIter = *this;
			return ++nextIter;
		}

	private:
		Vertex * currentVertex_;
		int visitedCount_;
	};


	class Utils
	{
	public:
		// Shoelace formula. 
		// The area formula is valid for any non-self-intersecting (simple) polygon, which can be convex or concave.
		// See https://en.wikipedia.org/wiki/Shoelace_formula
		// Input: polygon vertices (ordered), it should equal to the number of edges
		// Output: the area.
		static float CalculatePolygonArea(const std::vector<Point2d>& polygon)
		{
			float a = 0.0, b = 0.0;
			for (int i=0;i<polygon.size();++i)
			{
				if (i == polygon.size() - 1)
				{
					a += polygon[i].x_*polygon[0].y_;
					b += polygon[0].x_*polygon[i].y_;
				}
				else
				{
					a += polygon[i].x_*polygon[i + 1].y_;
					b += polygon[i + 1].x_*polygon[i].y_;
				}
			}
			return std::abs((a - b) / 2);
		}
	};

	
}

#endif