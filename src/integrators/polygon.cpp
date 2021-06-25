#include "integrators/polygon.h"

namespace PolyClip {


PolyClip::Polygon::Polygon(std::vector<Point2d> vertices)
{
	vertexNum_ = 0;
	for (int i = 0; i<vertices.size(); ++i)
	{
		if (i == 0)
		{
			Vertex* new_v = new Vertex(vertices[0]);
			if(vertices[0].hasZCoord) {
				new_v->z_ = vertices[0].z_;
				new_v->hasZCoord = true;
			}

			new_v->next_ = new_v;
			new_v->prev_ = new_v;

			start_ = new_v;
			end_ = new_v;
			vertexNum_++;
		}
		else
		{
			Insert(vertices[i], end_);
			end_ = end_->next_;
		}
	}
}

PolyClip::Polygon::Polygon(Polygon& poly)
{
	int i = 0;
	for (auto iter = poly.begin(); iter != poly.end(); ++iter, ++i)
	{
		if ((*iter).intersect_)		// kick out intersection points which were added by polygon clipping operation.
			continue;
		if (i == 0)
		{
			Vertex* new_v = new Vertex(Point2d((*iter).x_, (*iter).y_));
			/* Note: these member variables cannot be validly copy since they are modified by polygon clipping operation.
			*new_v->alpha_ = (*iter).alpha_;
			new_v->entryExit_ = (*iter).entryExit_;
			new_v->intersect_ = (*iter).intersect_;
			new_v->neighbour_ = (*iter).neighbour_;*/

			new_v->next_ = new_v;
			new_v->prev_ = new_v;

			start_ = new_v;
			end_ = new_v;
			vertexNum_++;
		}
		else
		{
			Vertex* new_v = new Vertex(Point2d((*iter).x_, (*iter).y_));
			/*new_v->alpha_ = (*iter).alpha_;
			new_v->entryExit_ = (*iter).entryExit_;
			new_v->intersect_ = (*iter).intersect_;
			new_v->neighbour_ = (*iter).neighbour_;*/

			Insert(new_v, end_);
			end_ = end_->next_;
		}
	}
}

PolyClip::Polygon& PolyClip::Polygon::operator=(Polygon& poly)
{
	if (this == &poly) {
		return *this;
	}
	// first, deallocate memory that 'this' used to hold
	while (start_)
		Remove(start_);

	// second, copy from poly by using copy constructor
	Polygon* copy_poly = new Polygon(poly);

	// update 
	vertexNum_ = copy_poly->vertexNum_;
	start_ = copy_poly->start_;
	end_ = copy_poly->end_;

	return *this;
}

PolyClip::Polygon::~Polygon()
{
	while (start_)
		Remove(start_);
}

void PolyClip::Polygon::Insert(Vertex* new_v, Vertex* pos)
{
	auto next_v = pos->next_;

	// modify new vertex
	new_v->next_ = next_v;
	new_v->prev_ = pos;

	// modify current vertex
	pos->next_ = new_v;

	// modify next vertex
	next_v->prev_ = new_v;

	vertexNum_++;
}

void PolyClip::Polygon::Insert(const Point2d& vertex, Vertex* pos)
{
	Vertex * new_v = new Vertex(vertex);
	if(vertex.hasZCoord) {
		new_v->hasZCoord = true;
		new_v->z_ = vertex.z_;
	}

	auto next_v = pos->next_;

	// modify new vertex
	new_v->next_ = next_v;
	new_v->prev_ = pos;

	// modify current vertex
	pos->next_ = new_v;

	// modify next vertex
	next_v->prev_ = new_v;

	vertexNum_++;
}

void PolyClip::Polygon::Insert(float x, float y, Vertex* pos)
{
	Vertex * new_v = new Vertex(x, y);
	new_v->hasZCoord = false;
	new_v->z_ = 0.f;

	auto next_v = pos->next_;

	// modify new vertex
	new_v->next_ = next_v;
	new_v->prev_ = pos;

	// modify current vertex
	pos->next_ = new_v;

	// modify next vertex
	next_v->prev_ = new_v;

	vertexNum_++;
}

void PolyClip::Polygon::Remove(Vertex* pos)
{
	if (vertexNum_ == 1)
	{
		start_ = nullptr;
		end_ = nullptr;
		auto neighbor_v = pos->neighbour_;
		if (neighbor_v)
			neighbor_v->neighbour_ = nullptr;
		delete pos;
	}
	else
	{
		auto pre_v = pos->prev_;
		auto next_v = pos->next_;
		auto neighbor_v = pos->neighbour_;

		// modify previous vertex
		if (pre_v)
			pre_v->next_ = next_v;
		// modify next vertex
		if (next_v)
			next_v->prev_ = pre_v;
		// modify neighbor vertex
		if (neighbor_v)
			neighbor_v->neighbour_ = nullptr;
		if (pos == start_)
			start_ = next_v;
		if (pos == end_)
			end_ = pre_v;
		delete pos;
	}
	vertexNum_--;
}

}