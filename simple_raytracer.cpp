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
	glm::vec3 normal_one;  // Normal at point_one
	glm::vec3 normal_two;  // Normal at point_two
	glm::vec3 normal_three; // Normal at point_three

};

std::vector<Triangle> triangleobjloader(std::string objfilename) {
	std::string inputfile = objfilename;
	tinyobj::ObjReader reader;

	if (!reader.ParseFromFile(inputfile)) {
		if (!reader.Error().empty()) {
			std::cerr << "TinyObjReader: " << reader.Error();
		}
	}
	/*
	if (!reader.Warning().empty()) {
		std::cout << "TinyObjReader: " << reader.Warning();
	} */

	auto& attrib = reader.GetAttrib();
	auto& shapes = reader.GetShapes();
	auto& materials = reader.GetMaterials();

	std::vector<Triangle> triangles;

	// store triangles and normals
	for (const auto& shape : shapes) {
		for (const auto& index : shape.mesh.indices) {
			glm::vec3 vertex(
				attrib.vertices[3 * index.vertex_index + 0],
				attrib.vertices[3 * index.vertex_index + 1],
				attrib.vertices[3 * index.vertex_index + 2]
			);
			glm::vec3 normal(
				attrib.normals[3 * index.normal_index + 0],
				attrib.normals[3 * index.normal_index + 1],
				attrib.normals[3 * index.normal_index + 2]

			);
			static std::vector<glm::vec3> vertices;
			vertices.push_back(vertex);
			static std::vector<glm::vec3> normals;
			normals.push_back(normal);
			if (vertices.size() == 3) {
				Triangle triangle;
				triangle.point_one = vertices[0];
				triangle.point_two = vertices[1];
				triangle.point_three = vertices[2];
				triangle.normal_one = normals[0];
				triangle.normal_two = normals[1];
				triangle.normal_three = normals[2];
				triangles.push_back(triangle);
				vertices.clear();
				normals.clear();
			}
		}
	}

	// Output the triangles
	/*
	for (const auto& triangle : triangles) {
		std::cout << "Triangle: \n";
		std::cout << "  Point One: (" << triangle.point_one.x << ", " << triangle.point_one.y << ", " << triangle.point_one.z << ")\n";
		std::cout << "  Point Two: (" << triangle.point_two.x << ", " << triangle.point_two.y << ", " << triangle.point_two.z << ")\n";
		std::cout << "  Point Three: (" << triangle.point_three.x << ", " << triangle.point_three.y << ", " << triangle.point_three.z << ")\n";
	}
	std::cout << "Size:\n";
	std::cout << triangles.size(); */
	return triangles;
}

inline float rayTriangleIntersection(const Ray* ray, const Triangle* triangle) {
	// Intersection of a ray with a triangle
	// This implementation uses the Möller–Trumbore intersection algorithm

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

glm::vec3 calculateBarycentricCoords(const Triangle& triangle, const glm::vec3& point) {
	// Calculate the vectors from point_one to the other two vertices of the triangle
	// These vectors represent the edges of the triangle and the vector from the
	// first vertex to the point of interest. They are used for determining the relative
	// position of the point within the triangle.
	glm::vec3 v0 = triangle.point_two - triangle.point_one;
	glm::vec3 v1 = triangle.point_three - triangle.point_one;
	glm::vec3 v2 = point - triangle.point_one;

	// The dot products are used to calculate the areas and angles between the vectors.
	// These values are crucial for the barycentric coordinates formula, which determines
	// how much of each vertex's influence is present at the point.
	float d00 = glm::dot(v0, v0); 
	float d01 = glm::dot(v0, v1); 
	float d11 = glm::dot(v1, v1);
	float d20 = glm::dot(v2, v0);
	float d21 = glm::dot(v2, v1);

	// Compute the denominator of the barycentric coordinates formula
	// It normalizes the coordinates, ensuring they sum to 1. 
	// This step ensures that the point lies within the triangle
	float denom = d00 * d11 - d01 * d01;

	// Compute the barycentric coordinates (v, w, u)
	// v and w are calculated using the dot products and the denominator
	float v = (d11 * d20 - d01 * d21) / denom;
	float w = (d00 * d21 - d01 * d20) / denom;
	// u is calculated to ensure the sum of the barycentric coordinates is 1
	float u = 1.0f - v - w;

	// Return the barycentric coordinates as a vec3
	return glm::vec3(u, v, w);
}


glm::vec3 interpolateNormal(const Triangle& triangle, const glm::vec3& barycentricCoords) {
	// By multiplying each vertex normal by its corresponding barycentric coordinate, 
	// I am weighting each normal by how much influence that vertex has at the point of interest.
	return glm::normalize(
		barycentricCoords.x * triangle.normal_one +
		barycentricCoords.y * triangle.normal_two +
		barycentricCoords.z * triangle.normal_three
	);
}

glm::vec3 phongIllumination(const Triangle& triangle, const Ray ray, const glm::vec3& lightPos, const glm::vec3& Llight, const glm::vec3& p, float ambientStrength, float s, float m, float distance) {
	// Phong illumination model

	// p = object color
	// Llight = Light Color
	// m = specular radius
	// s = specular Strength

	// rView is a constant factor for the reflection model, representing the light reflected to the view position
	// smaller number = less light reflected to view position = diffuse plays a smaller role
	constexpr float rView = 1.0f / glm::pi<float>();

	// Calculate the intersection point of the ray with the triangle
	glm::vec3 intersectionPoint = ray.origin + distance * ray.direction;

	// Calculate barycentric coordinates for the intersection point within the triangle
	glm::vec3 barycentricCoords = calculateBarycentricCoords(triangle, intersectionPoint);

	// Interpolate the normal at the intersection point using barycentric coordinates
	glm::vec3 n = interpolateNormal(triangle, barycentricCoords); // normal

	// Calculate the direction vector from the intersection point to the light source
	glm::vec3 l = glm::normalize(lightPos - intersectionPoint); // lightDirection

	// Calculate Diffuse Reflection
	// Diffuse reflection is based on Lambert's cosine law, which states that the intensity of light is proportional to the 
	// cosine of the angle between the light direction and the surface normal
	// Means: intensity of light is is higher if the angle is sharper
	// The dot product n * l represents this cosine value, and I use max to ensure it is non-negative
	// p * Llight represents the final color of object with light 
	glm::vec3 diffuse = rView * p * Llight * glm::max(glm::dot(n, l), 0.00f);

	// Calculate Ambient Reflection
	// Ambient reflection represents the constant illumination of the object by the environment
	// It is usually a small constant value added to ensure that objects are visible even when not directly lit
	// Higher AmbientStrenght = All of the Object brightens up more by the same amount
	glm::vec3 ambient = (1 / glm::pi<float>()) * ambientStrength * p * Llight;

	// Calculate Specular Reflection
	// Specular reflection represents the mirror-like reflection of light sources on shiny surfaces
	// It does not use the object color (p) because specular highlights are typically the color of the light source
	// Higher shininess (m) means a smaller specular highlight
	glm::vec3 v = glm::normalize(-ray.direction); // View Direction
	glm::vec3 r = glm::reflect(-l, n); // Reflect Direction

	// The specular term is calculated using the Phong reflection model
	// It is based on the dot product between the view direction and the reflection direction, raised to the power of the shininess factor (m)
	// s = specular strength. Smaller specular strength means less intensity
	glm::vec3 specular = Llight * s * glm::max(glm::dot(n, l), 0.00f) * glm::pow(glm::max(glm::dot(r, v), 0.0f), m);

	// Combine the three components (diffuse, specular, and ambient) to get the final color
	return diffuse + specular + ambient;
}


std::pair<glm::vec2, glm::vec3> rayIntersection(Ray ray,std::vector<Triangle> triangles, int point_x, int point_y){
	float distance_comparison = INFINITY;
	glm::vec3 color_point(0,0,0);
	for (int k = 0; k < triangles.size(); k++) {
		float f_distance = rayTriangleIntersection(&ray, &triangles[k]);
		if (f_distance != -INFINITY) {
			if (f_distance < distance_comparison) {
				distance_comparison = f_distance;

				// Initialize Phong Illumination with a Red Object
				glm::vec3 lightPos(300.0f, -300.0f, -1000.0f);
				glm::vec3 lightColor(1.0f, 1.0f, 1.0f); // White light
				glm::vec3 objectColor(1.0f, 0.0f, 0.0f); // Red object
				float ambientStrength = 0.2f;
				float specularStrength = 0.5f;
				float shininess = 15.0f;
				glm::vec3 color = phongIllumination(triangles[k], ray, lightPos, lightColor, objectColor, ambientStrength, specularStrength, shininess, f_distance);
				color = glm::clamp(color, 0.0f, 1.0f);
				color_point.x = int((color.x * 255));
				color_point.y = int((color.y * 255));
				color_point.z = int((color.z * 255));

				/*
				* Draw the Triangles based on the intersection distance
				int i_distance = int(f_distance * 40);
				glm::vec3 temp_color(i_distance * 8, i_distance *8, i_distance * 8);
				glm::vec3 clampcolortest = glm::vec3(f_distance, f_distance, f_distance);
				glm::vec3 clampedColor = glm::clamp(clampcolortest, 0.4f, 2.5f);
				clampedColor.x = int((clampedColor.x * 255));
				clampedColor.y = int((clampedColor.y * 255));
				clampedColor.z = int((clampedColor.z * 255));
				color_point = temp_color;
				color_point = clampedColor;
				*/
			}
		}
	}
	glm::vec2 image_point(point_x, point_y);
	return { image_point, color_point };
}

int main()
{
	// create a triangle
	Triangle triangle;
	triangle.point_one = glm::vec3(5.0f, 0.0f, 0.0f);
	triangle.point_two = glm::vec3(-5.0f, 0.0f, 0.0f);
	triangle.point_three = glm::vec3(0.0f, 8.0f, 0.0f);

	// create a ray starting at z = -100 so you can see objects 
	// which are centered at 0 0 0
	Ray ray;
	ray.direction = glm::vec3(0.0f, 0.0f,4000.0f);
	ray.origin = glm::vec3 (0.0f, 0.0f, -100.0f);

	// create image
	int image_width = 300;
	int image_height = 250;
	CImg<unsigned char> img(image_width, image_height, 1, 3);
	img.fill(0);

	// define default color for rays
	unsigned char color[] = { 255,128,64 };

	// load circle triangles from obj
	std::vector<Triangle> circle_triangles = triangleobjloader("sphere.obj");

	// load cube triangles from obj
	std::vector<Triangle> cube_triangles = triangleobjloader("cube.obj");
	//circle_triangles.insert(circle_triangles.end(), cube_triangles.begin(), cube_triangles.end());
	
	// add triangle
	//circle_triangles.push_back(triangle);

	// send rays out based from the center of the ray origin and intersect them with triangles
	std::vector<glm::vec2> image_points;
	std::vector<glm::vec3> image_colors;
	for (int i = -image_width / 2; i < image_width / 2; ++i)
	{
		for (int j = -image_height / 2; j < image_height / 2; ++j)
		{
			ray.direction.x = i;
			ray.direction.y = j;
			std::pair<glm::vec2, glm::vec3> points = rayIntersection(ray, circle_triangles, i + image_width / 2, j + image_height / 2);
			image_points.push_back(points.first);
			image_colors.push_back(points.second);
		}
	}

	// draw pixels found in ray intersection
	for (int i = 0; i < image_points.size(); i++) {
		color[0] = image_colors[i].x;
		color[1] = image_colors[i].y;
		color[2] = image_colors[i].z;;
		img.draw_point(image_points[i].x, image_points[i].y, color);
	}
	img.display("Simple Raytracer by Leon Lang");
}
