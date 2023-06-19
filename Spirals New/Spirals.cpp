#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"
#include <cmath>
#include <vector>
#include <string>
#include <complex>
#include <algorithm>


const float PI = 3.1415926535897932384626433832795f;

const float RED_TO_LUMINOSITY_RATIO = 0.299f;
const float GREEN_TO_LUMINOSITY_RATIO = 0.587f;
const float BLUE_TO_LUMINOSITY_RATIO = 0.114f;

const int32_t WIDTH = 1280;
const int32_t HEIGHT = 720;
const int32_t PIXEL_WIDTH = 1;
const int32_t PIXEL_HEIGHT = 1;
const bool FULL_SCREEN = false;
const bool VSYNC = false;

using Vector2 = olc::vf2d;
using complexf = std::complex<float>;
using complexi = std::complex<int32_t>;

int32_t iroundf(float f)
{
	return (int32_t)lroundf(f);
}
int32_t iround(double d)
{
	return (int32_t)round(d);
}

template<class T>
olc::v2d_generic<T> to_v2d(std::complex<T> c)
{
	return olc::v2d_generic<T>(c.real(), c.imag());
}
olc::vi2d to_vi2d(std::complex<float> c)
{
	return olc::vi2d(iroundf(c.real()), iroundf(c.imag()));
}
olc::vi2d to_vi2d(std::complex<double> c)
{
	return olc::vi2d(iround(c.real()), iround(c.imag()));
}

complexf to_complexf(olc::vi2d p)
{
	return complexf(p.x, p.y);
}

bool point_vs_rectangle(complexf p, complexf rect1, complexf rect2)
{
	float min_x = min(rect1.real(), rect2.real());
	float max_x = max(rect1.real(), rect2.real());
	float min_y = min(rect1.imag(), rect2.imag());
	float max_y = max(rect1.imag(), rect2.imag());

	return min_x < p.real() && p.real() < max_x &&
		   min_y < p.imag() && p.imag() < max_y;
}

bool segment_vs_segment(complexf line1_p1, complexf line1_p2, complexf line2_p1, complexf line2_p2)
{
	return true;

	if (line1_p2.real() < line1_p1.real())
	{
		swap(line1_p1, line1_p2);
	}
	if (line2_p2.real() < line2_p1.real())
	{
		swap(line2_p1, line2_p2);
	}

	// Check if any line is parallel to one of the axes

	const float EPSILON = 1e-6;

	bool parallel = false;
	float dx1 = line1_p2.real() - line1_p1.real();
	float dy1 = line1_p2.imag() - line1_p1.imag();
	float dx2 = line2_p2.real() - line2_p1.real();
	float dy2 = line2_p2.imag() - line2_p1.imag();
	for (float difference : { dx1, dy1, dx2, dy2 })
	{
		if (abs(difference) < EPSILON)
		{
			parallel = true;
		}
	}

	// TODO: Write a proper check if any line is parallel to any axis
	//if (parallel)
	//{
	//	return true;
	//}

	// If any line is parallel to one of the axes, make a rotation rotate them around a point.
	// This will not change the relative position of the lines,
	// so the rotated lines will perserve the intersection from the original ones,
	// but it will help us not to deal with a lot of edge cases
	
	// We should choose the angle such that after the rotation all lines are not parallel to any axis
	const complexf rotation_angle = complexf(cosf(PI / 6), sinf(PI / 6));
	std::vector<complexf*> points;
	points.push_back(&line1_p1);
	points.push_back(&line1_p2);
	points.push_back(&line2_p1);
	points.push_back(&line2_p2);
	
	for (complexf* point : points)
	{
		*point = (*point) - line1_p1;
		if (parallel)
		{
			(*point) *= rotation_angle;
		}
	}

	float line1_x1 = line1_p1.real();
	float line1_x2 = line1_p2.real();
	float line1_y1 = line1_p1.imag();
	float line1_y2 = line1_p2.imag();
	float line2_x1 = line2_p1.real();
	float line2_x2 = line2_p2.real();
	float line2_y1 = line2_p1.imag();
	float line2_y2 = line2_p2.imag();

	float k1 = (line1_y2 - line1_y1) / (line1_x2 - line1_x1);
	float k2 = (line2_y2 - line2_y1) / (line2_x2 - line2_x1);
	float c1 = line1_y1 - k1 * line1_x1;
	float c2 = line2_y1 - k2 * line2_x1;

	bool two_segments_on_one_line_intesect = false;
	struct CoordinateInfo
	{
		float x;
		int i;
	};
	std::vector<CoordinateInfo> coordinates = { {line1_x1, 0}, {line1_x2, 0}, {line2_x1, 1}, {line2_x2, 1} };
	std::sort(coordinates.begin(), coordinates.end(), [](CoordinateInfo c1, CoordinateInfo c2){
		return c1.x < c2.x;
	});
	two_segments_on_one_line_intesect = coordinates[0].i != coordinates[1].i;
	
	if (abs(k2 - k1) < EPSILON && abs(c2 - c1) < EPSILON && two_segments_on_one_line_intesect)
	{
		return true;
	}
	
	float x = (c2 - c1) / (k1 - k2);
	//float y = k1 * x + c1;

	if (min(line1_x1, line1_x2) <= x && x <= max(line1_x1, line1_x2) &&
		min(line2_x1, line2_x2) <= x && x <= max(line2_x1, line2_x2))
	{
		return true;
	}

	return false;
}

bool segment_vs_rectangle(complexf line1, complexf line2, complexf rect1, complexf rect2)
{
	if (point_vs_rectangle(line1, rect1, rect2) || point_vs_rectangle(line2, rect1, rect2))
	{
		return true;
	}

	float min_x = min(rect1.real(), rect2.real());
	float max_x = max(rect1.real(), rect2.real());
	float min_y = min(rect1.imag(), rect2.imag());
	float max_y = max(rect1.imag(), rect2.imag());

	std::vector<std::pair<complexf, complexf>> lines = {
		{{min_x, min_y}, {min_x, max_y}},
		{{min_x, max_y}, {max_x, max_y}},
		{{max_x, max_y}, {max_x, min_y}},
		{{max_x, min_y}, {min_x, min_y}}
	};

	for (std::pair<complexf, complexf> line : lines)
	{
		if (segment_vs_segment(line1, line2, line.first, line.second))
		{
			return true;
		}
	}

	return false;
}

bool triangle_vs_rectangle(complexf tri1, complexf tri2, complexf tri3, complexf rect1, complexf rect2)
{
	// TODO: Actually wrong, a triangle can intersect a rectangle, if the latter is contained inside the former
	return segment_vs_rectangle(tri1, tri2, rect1, rect2) ||
		   segment_vs_rectangle(tri2, tri3, rect1, rect2) ||
		   segment_vs_rectangle(tri3, tri1, rect1, rect2);
}

olc::Pixel gradient(const std::vector<olc::Pixel>& key_points, float time, float time_of_full_cycle)
{
	if (time < 0)
	{
		throw "Time (time) cannot be negative";
	}
	if (time_of_full_cycle <= 0)
	{
		throw "time_of_full_cycle cannot be negative or equal to 0";
	}

	time = fmodf(time, time_of_full_cycle);
	float part = time / time_of_full_cycle;

	int i1 = part * key_points.size();
	int i2 = (i1 + 1) % key_points.size();
	olc::Pixel color1 = key_points[i1];
	olc::Pixel color2 = key_points[i2];
	int32_t dr = color2.r - color1.r;
	int32_t dg = color2.g - color1.g;
	int32_t db = color2.b - color1.b;

	int32_t r = color1.r, g = color1.g, b = color1.b;
	float k = fmodf(part * key_points.size(), 1);
	r += dr * k;
	g += dg * k;
	b += db * k;

	return olc::Pixel(r, g, b);
}

// Custom Types
namespace ct
{
	struct Polygon
	{
		complexf center;
		std::vector<complexf> points;
		olc::Pixel color;
		
		Polygon()
		{
			center = complexf(0, 0);
			points = {};
			color = olc::Pixel(0, 0, 0);
		}

		Polygon(complexf center, std::vector<complexf> points, olc::Pixel color)
		{
			this->center = center;
			this->points = points;
			this->color = color;
		}

		int number_of_points()
		{
			return points.size();
		}

		void Rotate(float degrees)
		{
			float radians = degrees * PI / 180;
			complexf rotation = complexf(cosf(radians), sinf(radians));
			for (complexf& point : points)
			{
				point *= rotation;
			}
		}

		std::vector<complexf>* reconstruct_points(int number_of_points, float degrees = 0)
		{
			return &points;
		}
	};

	struct NGon : Polygon
	{
		float radius;

		NGon(int number_of_points, complexf center, float radius, olc::Pixel color, float degrees = 0) : Polygon(center, {}, color)
		{
			this->radius = radius;

			reconstruct_points(number_of_points, radius, degrees);
		}

		std::vector<complexf>* reconstruct_points(int number_of_points, float radius, float start_rotation_degrees = 0)
		{
			this->points = {};

			float alpha = 2 * PI / number_of_points;
			//float radius = side / sinf(alpha / 2);
			this->radius = radius;
			//complexf point = complexf(radius, 0);
			//if (number_of_points % 2 != 0)
			//{
			//	point *= complexf(cosf(PI / 2 + degrees), sinf(PI / 2 + degrees));
			//}
			complexf point = complexf(0, radius);
			float start_rotation_radians = PI * start_rotation_degrees / 180.0f;
			point *= complexf(cosf(start_rotation_radians), sinf(start_rotation_radians));
			complexf rotation = complexf(cosf(alpha), sinf(alpha));
			this->points.push_back(point);
			for (int i = 0; i < number_of_points - 1; ++i)
			{
				point *= rotation;
				this->points.push_back(point);
			}

			return &(this->points);
		}
	};
	struct Triangle : NGon
	{
		Triangle(complexf center, float side, olc::Pixel color, float degrees = 0) : NGon(3, center, side, color, degrees)
		{
			this->Rotate(90);
		}
	};
	struct Square : NGon
	{
		Square(complexf center, float side, olc::Pixel color, float degrees = 0) : NGon(4, center, side, color, degrees)
		{
			this->Rotate(45);
		}
	};
	struct Hexagon : NGon
	{
		Hexagon(complexf center, float side, olc::Pixel color, float degrees = 0) : NGon(6, center, side, color, degrees)
		{

		}
	};
}


class Engine : public olc::PixelGameEngine
{
private:
	float time_since_start = 0;

	olc::vi2d previousMousePos;

	olc::vi2d center_of_screen;
	complexf origin = complexf(0, 0);
	float default_zoom = 4.0f;
	float zoom;
	float max_zoom = 0.01f;
	float min_zoom = 100.0f;
	float zoom_multiplier = 1.25f;
	float zoom_per_turn = 0.1f;
	complexf center_in_world_space = complexf(0, 0);
	complexf point_in_center;
	
	complexf polygons_center = center_in_world_space;
	uint32_t polygons_count = 128;
	uint32_t vertices_count = 4;
	float rotation_of_spiral = 45.0f;
	float rotation_of_spiral_per_press = 0.25f;
	float rotation_degrees = 0.0f;
	float rotation_per_press = 1.0f;
	float rotation_per_polygon;
	std::vector<ct::NGon> polygons;

	uint32_t text_scale = 8;

	std::vector<olc::Pixel> rainbow_rgb = {
		olc::Pixel(255, 0, 0), olc::Pixel(255, 255, 0), olc::Pixel(0, 255, 0),
		olc::Pixel(0, 255, 255), olc::Pixel(0, 0, 255), olc::Pixel(255, 0, 255)
	};
	std::vector<olc::Pixel> rainbow_with_white = {
		olc::Pixel(255, 0, 0), olc::Pixel(255, 255, 0), olc::Pixel(0, 255, 0),
		olc::Pixel(0, 255, 255), olc::Pixel(0, 0, 255), olc::Pixel(255, 0, 255),
		olc::Pixel(255, 255, 255)
	};
	std::vector<olc::Pixel> rainbow = {
		olc::Pixel(255, 0, 0), olc::Pixel(255, 102, 0), olc::Pixel(255, 255, 0),
		olc::Pixel(0, 255, 0), olc::Pixel(0, 191, 255), olc::Pixel(0, 0, 255),
		olc::Pixel(139, 0, 255)
	};
	std::vector<olc::Pixel> rainbow_3 = {
		olc::Pixel(255, 0, 0), olc::Pixel(0, 255, 0), olc::Pixel(0, 0, 255)
	};

	std::vector<olc::Pixel> pretty_gradient = {
		olc::Pixel(236, 68, 176), olc::Pixel(17, 104, 226), olc::Pixel(17, 228, 181), olc::Pixel(17, 104, 226)
	};

	std::vector<olc::Pixel> red_blue = {
		olc::Pixel(255, 0, 0), olc::Pixel(0, 0, 255)
	};

	std::vector<olc::Pixel> green = {
		olc::Pixel(0, 255, 0)
	};

	std::vector<olc::Pixel> gradient_points = rainbow_rgb;
	std::vector<float> times;
	float time_of_full_cycle = 18;

	
	olc::vi2d to_screen_space(complexf point)
	{
		olc::vi2d from_center_of_screen = to_vi2d(zoom * (point - point_in_center));
		from_center_of_screen.y = -from_center_of_screen.y;

		return center_of_screen + from_center_of_screen;
	}

	complexf to_world_space(olc::vi2d point)
	{
		olc::vi2d from_center_of_screen = point - center_of_screen;
		from_center_of_screen.y = -from_center_of_screen.y;
		
		return complexf(to_complexf(from_center_of_screen) / zoom + point_in_center);
	}

	bool is_in_screen(olc::vi2d& point)
	{
		return 0 <= point.x && point.x < ScreenWidth() &&
			   0 <= point.y && point.y < ScreenHeight();
	}

	bool is_out_of_screen(olc::vi2d& point)
	{
		return !is_in_screen(point);
	}

	//void FillTriangle(complexf p1, complexf p2, complexf p3, olc::Pixel color = olc::WHITE)
	//{
	//	FillTriangle(, color);
	//}

	void DrawPolygon(ct::Polygon& polygon, int thickness)
	{
		float min_offset = (thickness % 2 == 0) ? -thickness / 2 + 1 : -thickness / 2;
		float max_offset = thickness / 2;
		
		for (size_t index = 0; index < polygon.points.size(); ++index)
		{
			complexf p1 = polygon.points[index];
			complexf p2 = polygon.points[(index + 1) % polygon.points.size()];
			float length1 = sqrtf(p1.real() * p1.real() + p1.imag() * p1.imag());
			float length2 = sqrtf(p2.real() * p2.real() + p2.imag() * p2.imag());
			complexf delta1_towards = p1 / length1 / zoom * min_offset;
			complexf delta1_outwards = p1 / length1 / zoom * max_offset;
			complexf delta2_towards = p2 / length1 / zoom * min_offset;
			complexf delta2_outwards = p2 / length1 / zoom * max_offset;

			complexf p1_inwards = p1 + delta1_towards;
			complexf p1_outwards = p1 + delta1_outwards;
			complexf p2_inwards = p2 + delta2_towards;
			complexf p2_outwards = p2 + delta2_outwards;

			p1_inwards += polygon.center;
			p1_outwards += polygon.center;
			p2_inwards += polygon.center;
			p2_outwards += polygon.center;

			// pss_i - point in screen space inwards
			olc::vi2d pss1_i = to_screen_space(p1_inwards);
			olc::vi2d pss1_o = to_screen_space(p1_outwards);
			olc::vi2d pss2_i = to_screen_space(p2_inwards);
			olc::vi2d pss2_o = to_screen_space(p2_outwards);

			if (is_out_of_screen(pss1_i) && is_out_of_screen(pss1_o) &&
				is_out_of_screen(pss2_i) && is_out_of_screen(pss2_o))
			{
				continue;
			}

			this->FillTriangle(pss1_i, pss1_o, pss2_i, polygon.color);
			//this->FillTriangle(pss1_i, pss2_o, pss2_o, polygon.color);
			//this->FillTriangle(pss2_i, pss2_o, pss1_i, polygon.color);
			this->FillTriangle(pss2_i, pss2_o, pss1_o, polygon.color);
		}
	}

	void FillNGon(ct::NGon& polygon)
	{
		for (int i = 0; i < polygon.points.size(); i++)
		{
			this->FillTriangle(to_screen_space(polygon.center), to_screen_space(polygon.points[i]), to_screen_space(polygon.points[(i+1) % polygon.points.size()]), polygon.color);
		}
	}

	void ReconstructPolygons(int vertices_count, float rotation_per_polygon)
	{
		for (int i = 0; i < polygons.size(); i++)
		{
			polygons[i].reconstruct_points(vertices_count, polygons[i].radius, rotation_of_spiral + rotation_per_polygon * i);
		}
	}
	
	bool OnUserCreate()
	{
		center_of_screen = olc::vi2d(ScreenWidth() / 2, ScreenHeight() / 2);
		point_in_center = center_in_world_space;
		zoom = default_zoom;

		rotation_per_polygon = rotation_degrees / polygons_count;

		for (int i = 0; i < polygons_count; ++i)
		{
			polygons.push_back(ct::NGon(vertices_count, polygons_center, (i + 1), olc::Pixel(0, 0, 0), rotation_of_spiral));
			polygons[i].Rotate(i * rotation_per_polygon);
			//times.push_back(float((sinf(float(i) * PI / polygons_count)+1)/2));
			times.push_back(float(i) / gradient_points.size() * time_of_full_cycle / polygons_count);
		}

		return true;
	}

	bool OnUserUpdate(float fElapsedTime)
	{
		FillRect({ 0, 0 }, { ScreenWidth(), ScreenHeight() }, olc::BLACK);
		time_since_start += fElapsedTime;
		
		int32_t scroll = GetMouseWheel();
		if (scroll > 0)
		{
			zoom = min(min_zoom, zoom * zoom_multiplier);
		}
		else if (scroll < 0)
		{
			zoom = max(max_zoom, zoom / zoom_multiplier);
		}

		if (GetMouse(0).bHeld)
		{
			olc::vi2d drag = olc::vi2d(GetMouseX(), GetMouseY()) - previousMousePos;
			drag.y = -drag.y;
			complexf move = (1 / zoom) * complexf(drag.x, drag.y);
			point_in_center -= move;
		}
		if (GetMouse(1).bHeld)
		{
			olc::vi2d drag = olc::vi2d(GetMouseX(), GetMouseY()) - previousMousePos;
		}
		previousMousePos = olc::vi2d(GetMouseX(), GetMouseY());

		if (GetKey(olc::C).bPressed)
		{
			point_in_center = polygons_center;
		}
		if (GetKey(olc::R).bPressed)
		{
			zoom = default_zoom;
		}


		bool increasedPolygons = GetKey(olc::W).bPressed;
		bool decreasedPolygons = GetKey(olc::S).bPressed;
		bool changedPolygons = false;
		if (increasedPolygons || decreasedPolygons)
		{
			int verticesIncrease = increasedPolygons - decreasedPolygons;
			vertices_count = max(3u, vertices_count + verticesIncrease);

			changedPolygons = true;
			
		}

		bool increasedRotation = GetKey(olc::UP).bHeld;
		bool decreasedRotation = GetKey(olc::DOWN).bHeld;
		bool changedRotation = false;
		if (increasedRotation || decreasedRotation)
		{
			float rotationIncrease = rotation_per_press * ((float)increasedRotation - (float)decreasedRotation);
			rotation_degrees += rotationIncrease;
			rotation_per_polygon = rotation_degrees / polygons_count;

			changedRotation = true;
		}

		bool increasedRotationSpiral = GetKey(olc::LEFT).bHeld;
		bool decreasedRotationSpiral = GetKey(olc::RIGHT).bHeld;
		bool changedRotationSpiral = false;
		if (increasedRotationSpiral || decreasedRotationSpiral)
		{
			float rotationIncrease = rotation_of_spiral_per_press * ((float)increasedRotationSpiral - (float)decreasedRotationSpiral);
			rotation_of_spiral += rotationIncrease;

			changedRotationSpiral = true;
		}

		if (changedPolygons || changedRotation || changedRotationSpiral)
		{
			ReconstructPolygons(vertices_count, rotation_per_polygon);
		}


		for (int i = 0; i < polygons.size(); i++)
		{
			polygons[i].color = gradient(gradient_points, times[i], time_of_full_cycle);

			float luminosity_factor = float(polygons.size() - i) / polygons.size();
			//float luminosity_factor = float(i + 1) / polygons.size();
			//polygons[i].color.r -= roundf(luminosity_factor * polygons[i].color.r);
			//polygons[i].color.g -= roundf(luminosity_factor * polygons[i].color.g);
			//polygons[i].color.b -= roundf(luminosity_factor * polygons[i].color.b);

			times[i] += fElapsedTime;
		}


		for (int i = polygons.size() - 1; i >= 1; i--)
		{
			ct::NGon polygon1 = polygons[i];
			ct::NGon polygon2 = polygons[i-1];
			for (int j = 0; j < polygons[i].points.size(); j++)
			{
				int j1 = j;
				int j2 = (j + 1) % polygons[i].points.size();

				complexf tri1_p1 = polygon1.points[j1];
				complexf tri1_p2 = polygon1.points[j2];
				complexf tri1_p3 = polygon2.points[j1];

				complexf tri2_p1 = polygon1.points[j2];
				complexf tri2_p2 = polygon2.points[j1];
				complexf tri2_p3 = polygon2.points[j2];
		
				
				if (triangle_vs_rectangle(tri1_p1, tri1_p2, tri1_p3, to_world_space(olc::vi2d(0, 0)), to_world_space(olc::vi2d(WIDTH, HEIGHT))))
				{
					this->FillTriangle(to_screen_space(tri1_p1), to_screen_space(tri1_p2), to_screen_space(tri1_p3), polygon1.color);
				}
				if (triangle_vs_rectangle(tri2_p1, tri2_p2, tri2_p3, to_world_space(olc::vi2d(0, 0)), to_world_space(olc::vi2d(WIDTH, HEIGHT))))
				{
					this->FillTriangle(to_screen_space(tri2_p1), to_screen_space(tri2_p2), to_screen_space(tri2_p3), polygon1.color);
				}
			}
			
			//DrawPolygon(polygons[i]);
		}
		FillNGon(polygons[0]);


		if (GetKey(olc::I).bHeld)
		{
			olc::Pixel c = polygons[0].color;
			olc::Pixel text_color = olc::Pixel(255 - c.r, 255 - c.g, 255 - c.b);
			DrawString(olc::vi2d(0, 0), std::to_string(vertices_count), text_color, text_scale);
		}

		//if (GetKey(olc::))

		//for (int i = 0; i < polygons.size(); i++)
		//{
		//	DrawPolygon(polygons[i], 3);
		//}

		return true;
	}
};

int main()
{
	Engine engine;
	if (engine.Construct(WIDTH, HEIGHT, PIXEL_WIDTH, PIXEL_HEIGHT, FULL_SCREEN, VSYNC))
	{
		engine.Start();
	}

	return 0;
}