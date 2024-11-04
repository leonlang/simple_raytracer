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
float rayTriangleIntersection(Ray ray, Triangle triangle) {
	// Step 1: Calculate vectors for two edges sharing triangle.point_one
	glm::vec3 p1p2 = triangle.point_two - triangle.point_one;
	glm::vec3 p1p3 = triangle.point_three - triangle.point_one;

	// Step 2: Calculate the determinant
	glm::vec3 pvec = glm::cross(ray.direction, p1p3);
	float det = glm::dot(p1p2, pvec);

	// Step 3: If the determinant is near zero, the ray lies in the plane of the triangle
	if (det < 0.000001)
		return -INFINITY;

	// Step 4: Calculate the inverse of the determinant
	float invDet = 1.0 / det;

	// Step 5: Calculate the distance from triangle.point_one to the ray origin
	glm::vec3 tvec = ray.origin - triangle.point_one;

	// Step 6: Calculate the u parameter and test bounds
	float u = glm::dot(tvec, pvec) * invDet;
	if (u < 0 || u > 1)
		return -INFINITY;

	// Step 7: Prepare to test the v parameter
	glm::vec3 qvec = glm::cross(tvec, p1p2);

	// Step 8: Calculate the v parameter and test bounds
	float v = glm::dot(ray.direction, qvec) * invDet;
	if (v < 0 || u + v > 1)
		return -INFINITY;

	// Step 9: Calculate t, the distance along the ray to the intersection point
	return glm::dot(p1p3, qvec) * invDet;
}

std::vector<Triangle> triangleobjloader(){
	std::string inputfile = "simple_sphere.obj";
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

	// ray at origin aber ich verändere die Strahlen,
	// so dass sie in alle richtungen gehen, nicht nur nach links und rechts
	// create a ray starting at origin 0 0 0 
	glm::vec3 vecOrig(0.0f, 0.0f, -100.0f);
	glm::vec3 vecDir(2.0f, 2.0f, 1.0f);

	Ray ray;
	ray.direction = vecDir;
	ray.origin = vecOrig;


	// create image 600x400
	int image_width = 600;
	int image_height = 400;
	CImg<unsigned char> img(image_width, image_height, 1, 3);

	// Set pixel values to 0 (color : black)
	img.fill(0);

	// define color for rays
	unsigned char color[] = { 255,128,64 };
	// go through each pixel in image and check if there is an intersection with triangle
	
	std::vector<Triangle> triangleto = triangleobjloader();

	// start timer to measure time
	auto start = std::chrono::high_resolution_clock::now();
	for (int i = -image_width/2; i < image_width/2; ++i)
	{
		for (int j = -image_height/2; j < image_height/2; ++j)
		{
			glm::vec3 vecNewDir(i, j, 100.0f);
			ray.direction = vecNewDir;
			/*if (rayTriangleIntersection(ray, triangle) != -INFINITY) {
				img.draw_point(i, j, color);
			} */
			for (int k = 0; k < triangleto.size(); k++){
				float f_distance = rayTriangleIntersection(ray, triangleto[k]);
				if (f_distance != -INFINITY) {
					int i_distance = int(f_distance * 70);
					color[0] = i_distance*3;
					color[1] = i_distance*3;
					color[2] = i_distance*3;
					img.draw_point(i+image_width/2, j+image_height/2, color);
					k = triangleto.size();
				}
			}
		}
	}

	// Display the image.
	img.display("Simple Raytracer by Leon Lang");

	// unsigned char purple[] = { 255,0,255 };        // Define a purple color
	// img.draw_text(100, 100, "Hello World", purple); // Draw a purple "Hello world" at coordinates (100,100).
            
	/**
	std::string inputfile = "cornell_box.obj";
	tinyobj::ObjReaderConfig reader_config;
	reader_config.mtl_search_path = "./"; // Path to material files

	tinyobj::ObjReader reader;

	if (!reader.ParseFromFile(inputfile, reader_config)) {
		if (!reader.Error().empty()) {
			std::cerr << "TinyObjReader: " << reader.Error();
		}
		exit(1);
	}

	if (!reader.Warning().empty()) {
		std::cout << "TinyObjReader: " << reader.Warning();
	}

	auto& attrib = reader.GetAttrib();
	auto& shapes = reader.GetShapes();
	auto& materials = reader.GetMaterials();

	// Loop over shapes
	for (size_t s = 0; s < shapes.size(); s++) {
		// Loop over faces(polygon)
		size_t index_offset = 0;
		for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
			size_t fv = size_t(shapes[s].mesh.num_face_vertices[f]);
			std::cout << f;
			// Loop over vertices in the face.
			for (size_t v = 0; v < fv; v++) {
				// access to vertex
				tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
				tinyobj::real_t vx = attrib.vertices[3 * size_t(idx.vertex_index) + 0];
				tinyobj::real_t vy = attrib.vertices[3 * size_t(idx.vertex_index) + 1];
				tinyobj::real_t vz = attrib.vertices[3 * size_t(idx.vertex_index) + 2];

				// Check if `normal_index` is zero or positive. negative = no normal data
				if (idx.normal_index >= 0) {
					tinyobj::real_t nx = attrib.normals[3 * size_t(idx.normal_index) + 0];
					tinyobj::real_t ny = attrib.normals[3 * size_t(idx.normal_index) + 1];
					tinyobj::real_t nz = attrib.normals[3 * size_t(idx.normal_index) + 2];
				}

				// Check if `texcoord_index` is zero or positive. negative = no texcoord data
				if (idx.texcoord_index >= 0) {
					tinyobj::real_t tx = attrib.texcoords[2 * size_t(idx.texcoord_index) + 0];
					tinyobj::real_t ty = attrib.texcoords[2 * size_t(idx.texcoord_index) + 1];
				}

				// Optional: vertex colors
				// tinyobj::real_t red   = attrib.colors[3*size_t(idx.vertex_index)+0];
				// tinyobj::real_t green = attrib.colors[3*size_t(idx.vertex_index)+1];
				// tinyobj::real_t blue  = attrib.colors[3*size_t(idx.vertex_index)+2];
			}
			index_offset += fv;

			// per-face material
			shapes[s].mesh.material_ids[f];
		}
	}
	**/
}
