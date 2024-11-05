// simple_raytracer.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <glm.hpp>
#include <gtc/matrix_transform.hpp>
#include <gtc/type_ptr.hpp>
#include "CImg.h"
#include <chrono>

using namespace cimg_library;
#define TINYOBJLOADER_IMPLEMENTATION // define this in only *one* .cc
// Optional. define TINYOBJLOADER_USE_MAPBOX_EARCUT gives robust trinagulation. Requires C++11
#define TINYOBJLOADER_USE_MAPBOX_EARCUT
#include "tiny_obj_loader.h"

struct Ray {
	glm::vec3 origin;
	glm::vec3 direction;
};

struct Triangle {
	glm::vec3 point_one;
	glm::vec3 point_two;
	glm::vec3 point_three;

};
// Intersection of a ray with a triangle
// This implementation uses the Möller–Trumbore intersection algorithm


inline float rayTriangleIntersection(const Ray* ray, const Triangle* triangle) {
	glm::vec3 p1p2 = triangle->point_two - triangle->point_one;
	glm::vec3 p1p3 = triangle->point_three - triangle->point_one;
	glm::vec3 pvec = glm::cross(ray->direction, p1p3);
	float det = glm::dot(p1p2, pvec);

	if (det < 0.000001f) return -INFINITY;

	float invDet = 1.0f / det;
	glm::vec3 tvec = ray->origin - triangle->point_one;
	float u = glm::dot(tvec, pvec) * invDet;
	if (u < 0.0f || u > 1.0f) return -INFINITY;

	glm::vec3 qvec = glm::cross(tvec, p1p2);
	float v = glm::dot(ray->direction, qvec) * invDet;
	if (v < 0.0f || u + v > 1.0f) return -INFINITY;

	return glm::dot(p1p3, qvec) * invDet;
}


std::vector<Triangle> triangleobjloader(std::string objfilename){
	std::string inputfile = objfilename;
	tinyobj::ObjReader reader;

	if (!reader.ParseFromFile(inputfile)) {
		if (!reader.Error().empty()) {
			std::cerr << "TinyObjReader: " << reader.Error();
		}
	}

	if (!reader.Warning().empty()) {
		std::cout << "TinyObjReader: " << reader.Warning();
	}

	auto& attrib = reader.GetAttrib();
	auto& shapes = reader.GetShapes();
	auto& materials = reader.GetMaterials();

	std::vector<Triangle> triangles;

	for (const auto& shape : shapes) {
		for (const auto& index : shape.mesh.indices) {
			glm::vec3 vertex(
				attrib.vertices[3 * index.vertex_index + 0],
				attrib.vertices[3 * index.vertex_index + 1],
				attrib.vertices[3 * index.vertex_index + 2]
			);

			static std::vector<glm::vec3> vertices;
			vertices.push_back(vertex);

			if (vertices.size() == 3) {
				Triangle triangle;
				triangle.point_one = vertices[0];
				triangle.point_two = vertices[1];
				triangle.point_three = vertices[2];
				triangle.point_one *= 25.0f;
				triangle.point_two *= 25.0f;
				triangle.point_three *= 25.0f;
				//triangle.point_two.z += 2.6;
				//triangle.point_one.z += 2.6;
				//triangle.point_three.z += 2.6;

				triangles.push_back(triangle);
				vertices.clear();
			}
		}
	}

	// Output the triangles
	for (const auto& triangle : triangles) {
		std::cout << "Triangle: \n";
		std::cout << "  Point One: (" << triangle.point_one.x << ", " << triangle.point_one.y << ", " << triangle.point_one.z << ")\n";
		std::cout << "  Point Two: (" << triangle.point_two.x << ", " << triangle.point_two.y << ", " << triangle.point_two.z << ")\n";
		std::cout << "  Point Three: (" << triangle.point_three.x << ", " << triangle.point_three.y << ", " << triangle.point_three.z << ")\n";
	}
	std::cout << "Size:\n";
	std::cout << triangles.size();
	return triangles;
}

std::pair<glm::vec2, glm::vec3> rayIntersection(Ray ray,std::vector<Triangle> triangles, int point_x, int point_y){
	float distance_comparison = INFINITY;
	glm::vec3 color_point(0,0,0);
	for (int k = 0; k < triangles.size(); k++) {
		float f_distance = rayTriangleIntersection(&ray, &triangles[k]);
		if (f_distance != -INFINITY) {
			if (f_distance < distance_comparison) {
				distance_comparison = f_distance;
				int i_distance = int(f_distance * 70);
				glm::vec3 temp_color(i_distance * 3, i_distance * 3, i_distance * 3);
				color_point = temp_color;
			}
		}
	}
	glm::vec2 image_point(point_x, point_y);
	return { image_point, color_point };
}
int main()
{
	// create a triangle
	glm::vec3 triangleP1(600.0f, 100.0f, 3.0f);
	glm::vec3 triangleP2(200.0f, 100.0f, 3.0f);
	glm::vec3 triangleP3(300.0f, 1000.0f, 3.0f);

	Triangle triangle;
	triangle.point_one = triangleP1;
	triangle.point_two = triangleP2;
	triangle.point_three = triangleP3;


	// create a ray starting at z = -100 so you can see objects 
	// which are centered at 0 0 0
	glm::vec3 vecOrig(0.0f, 0.0f, -100.0f);
	glm::vec3 vecDir(2.0f, 2.0f, 100.0f);

	Ray ray;
	ray.direction = vecDir;
	ray.origin = vecOrig;


	// create image 600x400
	int image_width = 300;
	int image_height = 250;
	CImg<unsigned char> img(image_width, image_height, 1, 3);

	// Set pixel values to 0 (color : black)
	img.fill(0);

	// define color for rays
	unsigned char color[] = { 255,128,64 };
	// go through each pixel in image and check if there is an intersection with triangle
	
	std::vector<Triangle> circle_triangles = triangleobjloader("sphere.obj");


	std::vector<glm::vec2> image_points;
	std::vector<glm::vec3> image_colors;
	for (int i = -image_width / 2; i < image_width / 2; ++i)
	{
		for (int j = -image_height / 2; j < image_height / 2; ++j)
		{
			ray.direction.x = i;
			ray.direction.y = j;
			/*if (rayTriangleIntersection(ray, triangle) != -INFINITY) {
				img.draw_point(i, j, color);
			} */
			std::pair<glm::vec2, glm::vec3> points = rayIntersection(ray, circle_triangles, i + image_width / 2, j + image_height / 2);
			image_points.push_back(points.first);
			image_colors.push_back(points.second);
		}
	}
	std::cout << "Size of image points" << image_points.size();
	// draw pixels found in circle
	for (int i = 0; i < image_points.size(); i++) {
		color[0] = image_colors[i].x;
		color[1] = image_colors[i].y;
		color[2] = image_colors[i].z;

		img.draw_point(image_points[i].x, image_points[i].y, color);
	}
	img.display("Simple Raytracer by Leon Lang");
}
