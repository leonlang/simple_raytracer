// simple_raytracer.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <glm.hpp>
#include <gtc/matrix_transform.hpp>
#include <gtc/type_ptr.hpp>
#include "CImg.h"
#include <chrono>
#include "Transformation.h"
#include "Object.h"


using namespace cimg_library;
#define TINYOBJLOADER_IMPLEMENTATION // define this in only *one* .cc
#define TINYOBJLOADER_USE_MAPBOX_EARCUT
#include "tiny_obj_loader.h"
/*
struct Ray {
	glm::vec3 origin;
	glm::vec3 direction;
};

struct Triangle {
	glm::vec4 pointOne;
	glm::vec4 pointTwo;
	glm::vec4 pointThree;
	glm::vec3 normalOne;  // Normal at pointOne
	glm::vec3 normalTwo;  // Normal at pointTwo
	glm::vec3 normalThree; // Normal at pointThree
};
*/
glm::vec3 calculateTriangleNormal(Triangle triangle) {
	glm::vec3 v1 = triangle.pointTwo - triangle.pointOne;
	glm::vec3 v2 = triangle.pointThree - triangle.pointOne;
	glm::vec3 normal = glm::cross(v1, v2);
	return glm::normalize(normal);
}

std::vector<Triangle> triangleObjLoader(std::string objFilename) {
	std::string inputFile = objFilename;
	tinyobj::ObjReader reader;

	if (!reader.ParseFromFile(inputFile)) {
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
				// Convert them to triangle Points and add 1 for homogeneous notation
				triangle.pointOne = glm::vec4(vertices[0], 1);
				triangle.pointTwo = glm::vec4(vertices[1], 1);
				triangle.pointThree = glm::vec4(vertices[2], 1);
				triangle.normalOne = normals[0];
				triangle.normalTwo = normals[1];
				triangle.normalThree = normals[2];
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
		std::cout << "  Point One: (" << triangle.pointOne.x << ", " << triangle.pointOne.y << ", " << triangle.pointOne.z << ")\n";
		std::cout << "  Point Two: (" << triangle.pointTwo.x << ", " << triangle.pointTwo.y << ", " << triangle.pointTwo.z << ")\n";
		std::cout << "  Point Three: (" << triangle.pointThree.x << ", " << triangle.pointThree.y << ", " << triangle.pointThree.z << ")\n";
	}
	std::cout << "Size:\n";
	std::cout << triangles.size(); */
	return triangles;
}

std::vector<Triangle> objTransform(std::vector<Triangle> triangles, glm::mat4 matrix) {
	// Triangle Matrix Multiplication
	std::vector<Triangle> obj;
	for (int k = 0; k < triangles.size(); k++) {
		Triangle t = triangles[k];
		t.pointOne = matrix * triangles[k].pointOne;
		t.pointTwo = matrix * triangles[k].pointTwo;
		t.pointThree = matrix * triangles[k].pointThree;
		obj.push_back(t);
	}
	return obj;
}


inline float rayTriangleIntersection(const Ray* ray, const Triangle* triangle) {
	// Intersection of a ray with a triangle
	// This implementation uses the Möller–Trumbore intersection algorithm

	// Triangle Point1 in Cartesian Form
	glm::vec3 tP1Cartesian = glm::vec3(triangle->pointOne) / triangle->pointOne.w;
	glm::vec3 tP2Cartesian = glm::vec3(triangle->pointTwo) / triangle->pointTwo.w;
	glm::vec3 tP3Cartesian = glm::vec3(triangle->pointThree) / triangle->pointThree.w;

	glm::vec3 p1p2 = tP2Cartesian - tP1Cartesian;
	glm::vec3 p1p3 = tP3Cartesian - tP1Cartesian;
	glm::vec3 pvec = glm::cross(ray->direction, p1p3);
	float det = glm::dot(p1p2, pvec);

	if (det < 0.000001f) return -INFINITY;

	float invDet = 1.0f / det;
	glm::vec3 tvec = ray->origin - tP1Cartesian;
	float u = glm::dot(tvec, pvec) * invDet;
	if (u < 0.0f || u > 1.0f) return -INFINITY;

	glm::vec3 qvec = glm::cross(tvec, p1p2);
	float v = glm::dot(ray->direction, qvec) * invDet;
	if (v < 0.0f || u + v > 1.0f) return -INFINITY;

	float t = glm::dot(p1p3, qvec) * invDet;

	// Check if triangle is behind ray
	if (t < 0.0f) return -INFINITY;
	return t;
}

glm::vec3 calculateBarycentricCoords(const Triangle& triangle, const glm::vec3& point) {
	// Calculate the vectors from pointOne to the other two vertices of the triangle
	// These vectors represent the edges of the triangle and the vector from the
	// first vertex to the point of interest. They are used for determining the relative
	// position of the point within the triangle.
	glm::vec3 tP1Cartesian = glm::vec3(triangle.pointOne) / triangle.pointOne.w;
	glm::vec3 tP2Cartesian = glm::vec3(triangle.pointTwo) / triangle.pointTwo.w;
	glm::vec3 tP3Cartesian = glm::vec3(triangle.pointThree) / triangle.pointThree.w;

	glm::vec3 v0 = tP2Cartesian - tP1Cartesian;
	glm::vec3 v1 = tP3Cartesian - tP1Cartesian;
	glm::vec3 v2 = point - tP1Cartesian;

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
		barycentricCoords.x * triangle.normalOne +
		barycentricCoords.y * triangle.normalTwo +
		barycentricCoords.z * triangle.normalThree
	);
}

glm::vec3 phongIllumination(const Triangle& triangle, const Ray ray, const glm::vec3& lightPos, const glm::vec3& lightColor, const glm::vec3& objectColor, float ambientStrength, float specularStrength, float shininess, float distance) {
	// Phong illumination model

	// objectColor = object color
	// lightColor = Light Color
	// shininess = specular radius
	// specularStrength = specular Strength

	// rView is a constant factor for the reflection model, representing the light reflected to the view position
	// smaller number = less light reflected to view position = diffuse plays a smaller role
	constexpr float rView = 1.0f / glm::pi<float>();

	// Calculate the intersection point of the ray with the triangle
	glm::vec3 intersectionPoint = ray.origin + distance * ray.direction;

	// Calculate barycentric coordinates for the intersection point within the triangle
	glm::vec3 barycentricCoords = calculateBarycentricCoords(triangle, intersectionPoint);

	// Interpolate the normal at the intersection point using barycentric coordinates
	// glm::vec3 n = interpolateNormal(triangle, barycentricCoords); // normal
	// test with normal
	glm::vec3 n = calculateTriangleNormal(triangle);
	// Calculate the direction vector from the intersection point to the light source
	glm::vec3 l = glm::normalize(lightPos - intersectionPoint); // lightDirection

	// Calculate Diffuse Reflection
	// Diffuse reflection is based on Lambert's cosine law, which states that the intensity of light is proportional to the 
	// cosine of the angle between the light direction and the surface normal
	// Means: intensity of light is is higher if the angle is sharper
	// The dot product n * l represents this cosine value, and I use max to ensure it is non-negative
	// objectColor * lightColor represents the final color of object with light 
	float dotProduct = glm::dot(n, l);
	if (dotProduct < 0.0f) {
		dotProduct = -dotProduct;
	}
	glm::vec3 diffuse = rView * objectColor * lightColor * glm::max(dotProduct, 0.00f);

	// Calculate Ambient Reflection
	// Ambient reflection represents the constant illumination of the object by the environment
	// It is usually a small constant value added to ensure that objects are visible even when not directly lit
	// Higher AmbientStrenght = All of the Object brightens up more by the same amount
	glm::vec3 ambient = (1 / glm::pi<float>()) * ambientStrength * objectColor * lightColor;

	// Calculate Specular Reflection
	// Specular reflection represents the mirror-like reflection of light sources on shiny surfaces
	// It does not use the object color (objectColor) because specular highlights are typically the color of the light source
	// Higher shininess means a smaller specular highlight
	glm::vec3 v = glm::normalize(-ray.direction); // View Direction
	glm::vec3 r = glm::reflect(-l, n); // Reflect Direction

	// The specular term is calculated using the Phong reflection model
	// It is based on the dot product between the view direction and the reflection direction, raised to the power of the shininess factor (shininess)
	// specularStrength = specular strength. Smaller specular strength means less intensity
	glm::vec3 specular = lightColor * specularStrength * glm::max(dotProduct, 0.00f) * glm::pow(glm::max(glm::dot(r, v), 0.0f), shininess);

	// Combine the three components (diffuse, specular, and ambient) to get the final color
	return diffuse + specular + ambient;
}

std::pair<glm::vec2, glm::vec3> rayIntersection(Ray ray, std::vector<Triangle> triangles, int pointX, int pointY, glm::vec3 lightPos) {
	float distanceComparison = INFINITY;
	glm::vec3 colorPoint(0, 0, 0);
	for (int k = 0; k < triangles.size(); k++) {
		float fDistance = rayTriangleIntersection(&ray, &triangles[k]);
		if (fDistance != -INFINITY) {
			if (fDistance < distanceComparison) {
				distanceComparison = fDistance;

				// Initialize Phong Illumination with a Red Object
				glm::vec3 lightColor(1.0f, 1.0f, 1.0f); // White light
				glm::vec3 objectColor(1.0f, 0.0f, 0.0f); // Red object
				float ambientStrength = 0.2f;
				float specularStrength = 0.5f;
				float shininess = 15.0f;
				glm::vec3 color = phongIllumination(triangles[k], ray, lightPos, lightColor, objectColor, ambientStrength, specularStrength, shininess, fDistance);
				color = glm::clamp(color, 0.0f, 1.0f);
				colorPoint.x = int((color.x * 255));
				colorPoint.y = int((color.y * 255));
				colorPoint.z = int((color.z * 255));

				/*
				* Draw the Triangles based on the intersection distance
				int iDistance = int(fDistance * 40);
				glm::vec3 tempColor(iDistance * 8, iDistance *8, iDistance * 8);
				glm::vec3 clampColorTest = glm::vec3(fDistance, fDistance, fDistance);
				glm::vec3 clampedColor = glm::clamp(clampColorTest, 0.4f, 2.5f);
				clampedColor.x = int((clampedColor.x * 255));
				clampedColor.y = int((clampedColor.y * 255));
				clampedColor.z = int((clampedColor.z * 255));
				colorPoint = tempColor;
				colorPoint = clampedColor;
				*/
			}
		}
	}
	glm::vec2 imagePoint(pointX, pointY);
	return { imagePoint, colorPoint };
}

int main()
{
	for (float angleDegree = 40; angleDegree < 360; angleDegree = angleDegree + 1000) {

		// Create a circle and get the x and z coordinates for a specific degree
		// in its radius. With this, we can spin the camera around the center.
		float radius = 100.0f; // Radius of the circle on which the camera moves
		float radians = glm::radians(angleDegree); // Convert angle from degrees to radians
		float circleX = radius * std::cos(radians); // Calculate x coordinate on the circle
		float circleZ = radius * std::sin(radians); // Calculate z coordinate on the circle

		// create a triangle
		Triangle triangle;
		triangle.pointOne = glm::vec4(5.0f, 0.0f, 0.0f, 1.0f);
		triangle.pointTwo = glm::vec4(-5.0f, 0.0f, 0.0f, 1.0f);
		triangle.pointThree = glm::vec4(0.0f, 8.0f, 0.0f, 1.0f);

		// create a ray 
		Ray ray(glm::vec3(0.0f, 0.0f, 400.0f));
		// ray.direction = glm::vec3(0.0f, 0.0f, 400.0f);
		// ray.origin = glm::vec3(0.0f, 0.0f, 0.0f);

		// define default color for rays
		unsigned char color[] = { 255,128,64 };

		// create image
		int imageWidth = 300;
		int imageHeight = 250;
		CImg<unsigned char> img(imageWidth, imageHeight, 1, 3);
		img.fill(0);
		for (int i = 0; i < imageWidth; i++) {
			for (int j = 0; j < imageHeight; j++) {
				img.draw_point(i, j, color);
			}
		}

		// create viewMatrix
		glm::mat4 viewMatrix = Transformation::createViewMatrix(glm::vec3(circleX, -100.f, circleZ), glm::vec3(glm::radians(50.f), glm::radians(angleDegree + 90.f), 0.f));

		// load circle triangles from obj
		std::vector<Triangle> circleTriangles = triangleObjLoader("sphere.obj");
		// load cube triangles from obj
		std::vector<Triangle> cubeTriangles = triangleObjLoader("cube.obj");

		// Transform Cube Triangles
		cubeTriangles = objTransform(cubeTriangles, Transformation::scaleObj(10.0f, 10.0f, 10.0f));
		cubeTriangles = objTransform(cubeTriangles, glm::inverse(viewMatrix));

		// combine objects
		circleTriangles.insert(circleTriangles.end(), cubeTriangles.begin(), cubeTriangles.end());

		// add triangle
		//circleTriangles.push_back(triangle);

		// send rays out based from the center of the ray origin and intersect them with triangles
		std::vector<glm::vec2> imagePoints;
		std::vector<glm::vec3> imageColors;
		glm::vec4 lightPos(200.0f, -300.0f, -1000.4f, 1.0f);
		glm::vec2 rayXY = glm::vec2(ray.direction.x, ray.direction.y);
		// lightPos = viewSpaceTransformation(angleDegree) * lightPos;
		for (int i = -imageWidth / 2; i < imageWidth / 2; ++i)
		{
			for (int j = -imageHeight / 2; j < imageHeight / 2; ++j)
			{
				ray.direction.x = i + rayXY.x;
				ray.direction.y = j + rayXY.y;

				std::pair<glm::vec2, glm::vec3> points = rayIntersection(ray, cubeTriangles, i + imageWidth / 2, j + imageHeight / 2, lightPos);
				imagePoints.push_back(points.first);
				imageColors.push_back(points.second);
			}
		}

		// draw pixels found in ray intersection
		for (int i = 0; i < imagePoints.size(); i++) {
			color[0] = imageColors[i].x;
			color[1] = imageColors[i].y;
			color[2] = imageColors[i].z;
			img.draw_point(imagePoints[i].x, imagePoints[i].y, color);
		}

		unsigned char lightBlue[] = { 173, 216, 230 }; // RGB values for light blue

		// Iterate through all the pixels
		cimg_forXY(img, x, y) {
			// Check if the pixel is black (all channels are 0)
			if (img(x, y, 0, 0) == 0 && img(x, y, 0, 1) == 0 && img(x, y, 0, 2) == 0) {
				// Change the pixel to light blue
				img(x, y, 0, 0) = lightBlue[0]; // Red channel
				img(x, y, 0, 1) = lightBlue[1]; // Green channel
				img(x, y, 0, 2) = lightBlue[2]; // Blue channel
			}
		}
		std::string imgName = "output";
		imgName += std::to_string(static_cast<int>(angleDegree));
		imgName += ".bmp";
		img.save_bmp(imgName.c_str()); // Use c_str() to get a const char* from std::string
		img.display("Simple Raytracer by Leon Lang");

	}
}
