#include "integrators/polygonclipping.h"

namespace PolyClip {

void PolyClip::PloygonOpration::DetectIntersection(Polygon& clip_poly, Polygon& sub_poly)
{
	int intersection_point_count = 0;
	// record current correct loop status, because insert operation will modify iterator. We do not use iterator to traversal.
	int loop1_count = 0;
	int loop1_total = clip_poly.vertexNum_;
	
	for (auto iter1 = clip_poly.begin(); loop1_count<loop1_total; ++loop1_count)
	{
		auto next_check = iter1.next();
		int loop2_count = 0;
		int loop2_total = sub_poly.vertexNum_;

		std::vector<VertexPtrDistance> vertices_insert_poly1;

		for (auto iter2 = sub_poly.begin(); loop2_count<loop2_total; ++loop2_count)
		{
			auto next_check2 = iter2.next();

			float alphaP = -1.0, alphaQ = -1.0;
			Point2d p1((*iter1).x_, (*iter1).y_);
			Point2d p2((*iter1).next_->x_, (*iter1).next_->y_);
			Point2d q1((*iter2).x_, (*iter2).y_);
			Point2d q2((*iter2).next_->x_, (*iter2).next_->y_);
			if(LineSegmentIntersection(p1,p2,q1,q2,alphaP,alphaQ))			// here may modify p,q original (x,y) position
			{
				intersection_point_count++;
				// update original data
				(*iter1).x_ = p1.x_;	(*iter1).y_ = p1.y_;
				(*iter1).next_->x_ = p2.x_;	(*iter1).next_->y_ = p2.y_;
				(*iter2).x_ = q1.x_;	(*iter2).y_ = q1.y_;
				(*iter2).next_->x_ = q2.x_;	(*iter2).next_->y_ = q2.y_;

				Vertex i1 = CreateVertex(p1, p2, alphaP);
				Vertex i2 = CreateVertex(q1, q2, alphaQ);
				//printf("Line1 (%f, %f)--(%f, %f) and Line2 (%f, %f)--(%f, %f) intersect at (%f, %f) on Line1 and (%f,%f) on Line2\n\n",
				//	p1.x_, p1.y_, p2.x_, p2.y_, q1.x_, q1.y_, q2.x_, q2.y_, i1.x_, i1.y_, i2.x_, i2.y_);
				Vertex * I1 = new Vertex(i1);
				Vertex * I2 = new Vertex(i2);
				I1->intersect_ = true;
				I2->intersect_ = true;
				I1->neighbour_ = I2;
				I2->neighbour_ = I1;
				sub_poly.Insert(I2, iter2.eval());

				// Note: simply inserting the intersection vertex into father list will destory correct order. 
				//clip_poly.Insert(I1, iter1.eval());
				float distance = (I1->x_ - (*iter1).x_)*(I1->x_ - (*iter1).x_) + (I1->y_ - (*iter1).y_)*(I1->y_ - (*iter1).y_);		// calculate the distance between intersection point and edge start point 
				VertexPtrDistance vpt(I1, distance);
				vertices_insert_poly1.push_back(vpt);
			}
			iter2 = next_check2;
		}

		// sort objects which are waiting for insert
		std::sort(vertices_insert_poly1.begin(), vertices_insert_poly1.end(), SortVertexPtrDistance());
		auto temp_cur = iter1.eval();
		for (auto item : vertices_insert_poly1)
		{
			clip_poly.Insert(item.ptr, temp_cur);
			temp_cur = temp_cur->next_;
		}
		iter1 = next_check;
	}
	assert(intersection_point_count % 2 == 0 && "The intersection points should be even!");
}

bool PolyClip::PloygonOpration::Mark(Polygon& clip_poly, Polygon& sub_poly, /*out parameter*/ std::vector<std::vector<Point2d>>& possible_result, int markType,
										bool *didIntersect, int *indicator)
{
	bool noIntersection = true;
	int innerIndicator = 0;

	//---------------------- polygon 1 --------------------------
	bool status;		// false: exit, true: entry
	auto p0 = *(clip_poly.begin());
	if (Pnpoly(Point2d(p0.x_, p0.y_), sub_poly))
	{
		status = false;
		innerIndicator = 1;
	}
	else
		status = true;

	for(auto iter = clip_poly.begin(); iter != clip_poly.end(); ++iter)
	{
		if((*iter).intersect_)
		{
			(*iter).entryExit_ = status;
			status = !status;
			if (noIntersection)	noIntersection = false;		// have intersection point
		}
	}

	//---------------------- polygon 2 --------------------------
	p0 = *(sub_poly.begin());
	if (Pnpoly(Point2d(p0.x_, p0.y_), clip_poly))
	{
		status = false;
		innerIndicator = 2;
	}
	else
		status = true;

	for (auto iter = sub_poly.begin(); iter != sub_poly.end(); ++iter)
	{
		if ((*iter).intersect_)
		{
			(*iter).entryExit_ = status;
			status = !status;
			if (noIntersection)	noIntersection = false;		// have intersection point
		}
	}

	*didIntersect = !noIntersection;
	*indicator = innerIndicator;

	// no intersection
	if (noIntersection && innerIndicator == 0)
	{
		if (markType == MarkIntersection)	
			return false;
		if(markType == MarkUnion)
		{
			std::vector<Point2d> poly;
			for (auto iter = clip_poly.begin(); iter != clip_poly.end(); ++iter)
				poly.push_back(Point2d((*iter).x_, (*iter).y_));
			possible_result.push_back(poly);
			poly.clear();
			for (auto iter = sub_poly.begin(); iter != sub_poly.end(); ++iter)
				poly.push_back(Point2d((*iter).x_, (*iter).y_));
			possible_result.push_back(poly);
			return false;
		}
		if(markType == MarkDifferentiate)
		{
			std::vector<Point2d> poly;
			for (auto iter = clip_poly.begin(); iter != clip_poly.end(); ++iter)
				poly.push_back(Point2d((*iter).x_, (*iter).y_));
			possible_result.push_back(poly);
			return false;
		}
	}
	// polygon 1 is inside polygon 2
	if (noIntersection && innerIndicator == 1)
	{
		if(markType == MarkIntersection)
		{
			std::vector<Point2d> poly;
			for (auto iter = clip_poly.begin(); iter != clip_poly.end(); ++iter)
				poly.push_back(Point2d((*iter).x_, (*iter).y_));
			possible_result.push_back(poly);
			return false;
		}
		if(markType == MarkUnion)
		{
			std::vector<Point2d> poly;
			for (auto iter = sub_poly.begin(); iter != sub_poly.end(); ++iter)
				poly.push_back(Point2d((*iter).x_, (*iter).y_));
			possible_result.push_back(poly);
			return false;
		}
		if(markType == MarkDifferentiate)
		{
			std::vector<Point2d> poly;
			for (auto iter = sub_poly.begin(); iter != sub_poly.end(); ++iter)
				poly.push_back(Point2d((*iter).x_, (*iter).y_));
			possible_result.push_back(poly);
			poly.clear();
			for (auto iter = clip_poly.begin(); iter != clip_poly.end(); ++iter)
				poly.push_back(Point2d((*iter).x_, (*iter).y_));
			possible_result.push_back(poly);
			return false;
		}
	}
	// polygon 2 is inside polygon 1
	if (noIntersection && innerIndicator == 2)
	{
		if (markType == MarkIntersection)
		{
			std::vector<Point2d> poly;
			for (auto iter = sub_poly.begin(); iter != sub_poly.end(); ++iter)
				poly.push_back(Point2d((*iter).x_, (*iter).y_));
			possible_result.push_back(poly);
			return false;
		}
		if (markType == MarkUnion)
		{
			std::vector<Point2d> poly;
			for (auto iter = clip_poly.begin(); iter != clip_poly.end(); ++iter)
				poly.push_back(Point2d((*iter).x_, (*iter).y_));
			possible_result.push_back(poly);
			return false;
		}
		if (markType == MarkDifferentiate)
		{
			std::vector<Point2d> poly;
			for (auto iter = clip_poly.begin(); iter != clip_poly.end(); ++iter)
				poly.push_back(Point2d((*iter).x_, (*iter).y_));
			possible_result.push_back(poly);
			poly.clear();
			for (auto iter = sub_poly.begin(); iter != sub_poly.end(); ++iter)
				poly.push_back(Point2d((*iter).x_, (*iter).y_));
			possible_result.push_back(poly);
			return false;
		}
	}
	return true;
}

std::vector<std::vector<PolyClip::Point2d>> PolyClip::PloygonOpration::ExtractIntersectionResults(Polygon& clip_poly)
{
	// first, detect all of unprocessed intersection points
	std::vector<Vertex*> unprocessed_intersection_points;
	for(auto iter = clip_poly.begin();iter!=clip_poly.end();++iter)
	{
		if ((*iter).intersect_ && !(*iter).processd)
			unprocessed_intersection_points.push_back(iter.eval());
	}

	// extract polygons
	std::vector<std::vector<PolyClip::Point2d>> results;
	for (auto ptr : unprocessed_intersection_points)
	{
		if (ptr->processd) continue;

		ptr->processd = true;
		auto st = ptr;
		auto current = ptr;
		std::vector<PolyClip::Point2d> poly;

		auto ver = PolyClip::Point2d(current->x_, current->y_);
		ver.hasZCoord = true;
		ver.z_ = current->prev_->z_;
		poly.push_back(ver);

		do
		{
			if (current->entryExit_)
			{
				do
				{
					float temp = current->z_;

					current = current->next_;
					auto ver = PolyClip::Point2d(current->x_, current->y_);
					ver.hasZCoord = true;
					ver.z_ = temp;

					poly.push_back(ver);
				} while (!current->intersect_);
			}
			else
			{
				do
				{
					float temp = current->z_;

					current = current->prev_;
					auto ver = PolyClip::Point2d(current->x_, current->y_);
					ver.hasZCoord = true;
					ver.z_ = temp;
					
					poly.push_back(ver);
				} while (!current->intersect_);
			}
			current->processd = true;

			float temp = current->z_;
			current = current->neighbour_;
			current->z_ = temp;
			current->hasZCoord = true;

			current->processd = true;
			
		} while (current != st);

		results.push_back(poly);

	}

	return results;
}

std::vector<std::vector<PolyClip::Point2d>> PolyClip::PloygonOpration::ExtractUnionResults(Polygon& clip_poly)
{
	// first, detect all of unprocessed intersection points
	std::vector<Vertex*> unprocessed_intersection_points;
	for (auto iter = clip_poly.begin(); iter != clip_poly.end(); ++iter)
	{
		if ((*iter).intersect_ && !(*iter).processd)
			unprocessed_intersection_points.push_back(iter.eval());
	}

	// extract polygons
	std::vector<std::vector<PolyClip::Point2d>> results;
	for (auto ptr : unprocessed_intersection_points)
	{
		if (ptr->processd) continue;

		ptr->processd = true;
		auto st = ptr;
		auto current = ptr;
		std::vector<PolyClip::Point2d> poly;

		auto ver = PolyClip::Point2d(current->x_, current->y_);
		ver.hasZCoord = true;
		ver.z_ = current->prev_->z_;
		poly.push_back(ver);
		int idx = 0;
		do
		{
			if (current->entryExit_)
			{
				do
				{
					float temp = current->z_;

					current = current->prev_;
					auto ver = PolyClip::Point2d(current->x_, current->y_);
					ver.hasZCoord = true;
					ver.z_ = temp;

					poly.push_back(ver);

				} while (!current->intersect_);
			}
			else
			{
				do
				{
					float temp = current->z_;

					current = current->next_;
					auto ver = PolyClip::Point2d(current->x_, current->y_);
					ver.hasZCoord = true;
					ver.z_ = temp;

					poly.push_back(ver);
				} while (!current->intersect_);
			}
			current->processd = true;
			current = current->neighbour_;
			current->processd = true;
			idx++;
		} while (current != st && idx <= clip_poly.vertexNum_);

		results.push_back(poly);
	}

	return results;
}

std::vector<std::vector<PolyClip::Point2d>> PolyClip::PloygonOpration::ExtractDifferentiateResults(Polygon& clip_poly)
{
	// first, detect all of unprocessed intersection points
	std::vector<Vertex*> unprocessed_intersection_points;
	for (auto iter = clip_poly.begin(); iter != clip_poly.end(); ++iter)
	{
		if ((*iter).intersect_ && !(*iter).processd)
			unprocessed_intersection_points.push_back(iter.eval());
	}

	// extract polygons
	std::vector<std::vector<PolyClip::Point2d>> results;
	for (auto ptr : unprocessed_intersection_points)
	{
		if (ptr->processd) continue;

		bool self = true;
		ptr->processd = true;
		auto st = ptr;
		auto current = ptr;
		std::vector<PolyClip::Point2d> poly;
		poly.push_back(PolyClip::Point2d(current->x_, current->y_));
		int idx = 0;
		do
		{
			if (current->entryExit_)
			{
				do
				{
					if (self)
						current = current->prev_;
					else
						current = current->next_;
					poly.push_back(PolyClip::Point2d(current->x_, current->y_));
				} while (!current->intersect_);
			}
			else
			{
				do
				{
					if(self)
						current = current->next_;
					else current = current->prev_;
					poly.push_back(PolyClip::Point2d(current->x_, current->y_));
				} while (!current->intersect_);
			}
			self = !self;
			current->processd = true;
			current = current->neighbour_;
			current->processd = true;

			idx++;
		} while (current != st && idx <= clip_poly.vertexNum_*10);

		results.push_back(poly);
	}

	return results;
}


bool PolyClip::PloygonOpration::LineSegmentIntersection(
	Point2d& p1, Point2d& p2,
	Point2d& q1, Point2d& q2,
	float& alphaP, float& alphaQ)
{
	const float perturbation = static_cast<float>(1.001);
	// float perturbation = rand()/float(RAND_MAX/2.0) - 1;
	// perturbation = 1.f + perturbation * 1e-3;
	// std::cout << std::to_string(perturbation) << std::endl;

	Point2d vec_p1_q1(p1.x_ - q1.x_, p1.y_ - q1.y_);
	Point2d vec_p2_q1(p2.x_ - q1.x_, p2.y_ - q1.y_);
	Point2d vec_q2_q1(q2.x_ - q1.x_, q2.y_ - q1.y_);
	// vec_q2_q1 rotation 90
	Point2d vec_q2_q1_90(-vec_q2_q1.y_, vec_q2_q1.x_);
	// dot product
	float WEC_P1 = vec_p1_q1.x_*vec_q2_q1_90.x_ + vec_p1_q1.y_*vec_q2_q1_90.y_;
	float WEC_P2 = vec_p2_q1.x_*vec_q2_q1_90.x_ + vec_p2_q1.y_*vec_q2_q1_90.y_;
	
	if(WEC_P1 == 0.0)
	{
		// add perturbation
		p1.x_ = (p1.x_ - p2.x_)*perturbation + p2.x_;
		p1.y_ = (p1.y_ - p2.y_)*perturbation + p2.y_;
		vec_p1_q1.x_ = p1.x_ - q1.x_;	vec_p1_q1.y_ = p1.y_ - q1.y_;
		WEC_P1 = vec_p1_q1.x_*vec_q2_q1_90.x_ + vec_p1_q1.y_*vec_q2_q1_90.y_;
		WEC_P2 = vec_p2_q1.x_*vec_q2_q1_90.x_ + vec_p2_q1.y_*vec_q2_q1_90.y_;
	}
	if(WEC_P2 == 0.0)
	{
		// add perturbation
		p2.x_ = (p2.x_ - p1.x_)*perturbation + p1.x_;
		p2.y_ = (p2.y_ - p1.y_)*perturbation + p1.y_;
		vec_p2_q1.x_ = p2.x_ - q1.x_;	vec_p2_q1.y_ = p2.y_ - q1.y_;
		WEC_P1 = vec_p1_q1.x_*vec_q2_q1_90.x_ + vec_p1_q1.y_*vec_q2_q1_90.y_;
		WEC_P2 = vec_p2_q1.x_*vec_q2_q1_90.x_ + vec_p2_q1.y_*vec_q2_q1_90.y_;
	}
	if(WEC_P1*WEC_P2<0)
	{
		Point2d vec_q1_p1(q1.x_ - p1.x_, q1.y_ - p1.y_);
		Point2d vec_q2_p1(q2.x_ - p1.x_, q2.y_ - p1.y_);
		Point2d vec_p2_p1(p2.x_ - p1.x_, p2.y_ - p1.y_);
		// vec_p2_p1 rotation 90
		Point2d vec_p2_p1_90(-vec_p2_p1.y_, vec_p2_p1.x_);
		// dot product
		float WEC_Q1 = vec_q1_p1.x_*vec_p2_p1_90.x_ + vec_q1_p1.y_*vec_p2_p1_90.y_;
		float WEC_Q2 = vec_q2_p1.x_*vec_p2_p1_90.x_ + vec_q2_p1.y_*vec_p2_p1_90.y_;
		
		if(WEC_Q1 == 0)
		{
			// add perturbation
			q1.x_ = (q1.x_ - q2.x_)*perturbation + q2.x_;
			q1.y_ = (q1.y_ - q2.y_)*perturbation + q2.y_;
			vec_q1_p1.x_ = q1.x_ - p1.x_;	vec_q1_p1.y_ = q1.y_ - p1.y_;
			WEC_Q1 = vec_q1_p1.x_*vec_p2_p1_90.x_ + vec_q1_p1.y_*vec_p2_p1_90.y_;
			WEC_Q2 = vec_q2_p1.x_*vec_p2_p1_90.x_ + vec_q2_p1.y_*vec_p2_p1_90.y_;
		}
		if(WEC_Q2 == 0)
		{
			// add perturbation
			q2.x_ = (q2.x_ - q1.x_)*perturbation + q1.x_;
			q2.y_ = (q2.y_ - q1.y_)*perturbation + q1.y_;
			vec_q2_p1.x_ = q2.x_ - p1.x_;	vec_q2_p1.y_ = q2.y_ - p1.y_;
			WEC_Q1 = vec_q1_p1.x_*vec_p2_p1_90.x_ + vec_q1_p1.y_*vec_p2_p1_90.y_;
			WEC_Q2 = vec_q2_p1.x_*vec_p2_p1_90.x_ + vec_q2_p1.y_*vec_p2_p1_90.y_;
		}
		if(WEC_Q1*WEC_Q2<=0)
		{
			alphaP = WEC_P1 / (WEC_P1 - WEC_P2);
			alphaQ = WEC_Q1 / (WEC_Q1 - WEC_Q2);
			return true;
		}
	}

	return false;
}


//bool PolyClip::PloygonOpration::LineSegmentIntersection(
//	Point2d& p1, Point2d& p2,
//	Point2d& q1, Point2d& q2,
//	float& alphaP, float& alphaQ)
//{
//	// one side
//	Point2d vec_p1_q1(p1.x_ - q1.x_, p1.y_ - q1.y_);
//	Point2d vec_p2_q1(p2.x_ - q1.x_, p2.y_ - q1.y_);
//	Point2d vec_q2_q1(q2.x_ - q1.x_, q2.y_ - q1.y_);
//	// vec_q2_q1 rotation 90
//	Point2d vec_q2_q1_90(-vec_q2_q1.y_, vec_q2_q1.x_);
//	// dot product
//	float WEC_P1 = vec_p1_q1.x_*vec_q2_q1_90.x_ + vec_p1_q1.y_*vec_q2_q1_90.y_;
//	float WEC_P2 = vec_p2_q1.x_*vec_q2_q1_90.x_ + vec_p2_q1.y_*vec_q2_q1_90.y_;
//
//	// other side
//	Point2d vec_q1_p1(q1.x_ - p1.x_, q1.y_ - p1.y_);
//	Point2d vec_q2_p1(q2.x_ - p1.x_, q2.y_ - p1.y_);
//	Point2d vec_p2_p1(p2.x_ - p1.x_, p2.y_ - p1.y_);
//	// vec_p2_p1 rotation 90
//	Point2d vec_p2_p1_90(-vec_p2_p1.y_, vec_p2_p1.x_);
//	// dot product
//	float WEC_Q1 = vec_q1_p1.x_*vec_p2_p1_90.x_ + vec_q1_p1.y_*vec_p2_p1_90.y_;
//	float WEC_Q2 = vec_q2_p1.x_*vec_p2_p1_90.x_ + vec_q2_p1.y_*vec_p2_p1_90.y_;
//	
//	if(WEC_P1 == 0.0 && WEC_P2 == 0)
//	{
//		return false;	// line p1-p2 and line q1-q2 are overlap
//	}
//	else if((WEC_P1 == 0.0 && WEC_P2 !=0) || (WEC_P2 == 0.0 && WEC_P1 != 0))
//	{
//		if (WEC_Q1*WEC_Q2 < 0)
//		{
//			alphaP = WEC_P1 / (WEC_P1 - WEC_P2);
//			alphaQ = WEC_Q1 / (WEC_Q1 - WEC_Q2);
//			return true;
//		}
//	}
//	else if(WEC_P1*WEC_P2<0)
//	{
//		assert(!(WEC_Q1 == 0 && WEC_Q2 == 0));
//		if (WEC_Q1*WEC_Q2 <= 0)
//		{
//			alphaP = WEC_P1 / (WEC_P1 - WEC_P2);
//			alphaQ = WEC_Q1 / (WEC_Q1 - WEC_Q2);
//			return true;
//		}
//	}
//
//	return false;
//}

PolyClip::Vertex PolyClip::PloygonOpration::CreateVertex(const Point2d& p1, const Point2d& p2, float alpha)
{
	float x = alpha*(p2.x_ - p1.x_) + p1.x_;
	float y = alpha*(p2.y_ - p1.y_) + p1.y_;
	return Vertex(x, y);
}

bool PolyClip::PloygonOpration::Pnpoly(const Point2d& p, Polygon& polygon)
{
	bool c = false;
	
	for (auto iter = polygon.begin(); iter != polygon.end(); ++iter)
	{
		Point2d p1((*iter).x_, (*iter).y_);
		Point2d p2((*iter).next_->x_, (*iter).next_->y_);

		if (((p2.y_ > p.y_) != (p1.y_ > p.y_)) &&
			(p.x_ < (p1.x_ - p2.x_)*(p.y_ - p2.y_) / (p1.y_ - p2.y_) + p2.x_))
			c = !c;
	}

	return c;
	
}

}