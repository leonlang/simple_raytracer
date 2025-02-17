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

// Notes:
// Each Function presents an introduction. For example: Triangle Intersection corresponds to the same section in the report

struct ImageData { std::vector<glm::vec2> imagePoints; std::vector<glm::vec3> imageColors; };

// Normal Interpolation
// Concept for the Algorithm: https://www.khronos.org/opengl/wiki/Calculating_a_Surface_Normal
glm::vec3 calculateTriangleNormal(const Triangle* triangle) {
	glm::vec3 v1 = triangle->pointTwo - triangle->pointOne;
	glm::vec3 v2 = triangle->pointThree - triangle->pointOne;
	glm::vec3 normal = glm::cross(v1, v2);
	return glm::normalize(normal);
}

// Triangle Intersection
// This implementation uses the Möller–Trumbore intersection algorithm
// Concept for the Algorithm: https://www.graphics.cornell.edu/pubs/1997/MT97.pdf
inline float rayTriangleIntersection(const Ray* ray, const Triangle* triangle) {
	// Intersection of a ray with a triangle
	// Triangle Point1 in Cartesian Form
	glm::vec3 tP1Cartesian = glm::vec3(triangle->pointOne) / triangle->pointOne.w;
	glm::vec3 tP2Cartesian = glm::vec3(triangle->pointTwo) / triangle->pointTwo.w;
	glm::vec3 tP3Cartesian = glm::vec3(triangle->pointThree) / triangle->pointThree.w;

	// edge vectors
	glm::vec3 p1p2 = tP2Cartesian - tP1Cartesian;
	glm::vec3 p1p3 = tP3Cartesian - tP1Cartesian;

	// determinant
	glm::vec3 pvec = glm::cross(ray->direction, p1p3);
	float det = glm::dot(p1p2, pvec);

	if (fabs(det) < 1e-12f) return -INFINITY;

	// inverse Determinant
	float invDet = 1.0f / det;
	// distance vector
	glm::vec3 tvec = ray->origin - tP1Cartesian;
	// u and v paramter
	float u = glm::dot(tvec, pvec) * invDet;
	if (u < 0.0f || u > 1.0f) return -INFINITY;
	glm::vec3 qvec = glm::cross(tvec, p1p2);
	float v = glm::dot(ray->direction, qvec) * invDet;
	if (v < 0.0f || u + v > 1.0f) return -INFINITY;
	// intersection distance
	float t = glm::dot(p1p3, qvec) * invDet;

	// Check if triangle is behind ray
	if (t < 0.0f) return -INFINITY;
	return t;
}

// Phong Shading
// Concept for the Algorithm: https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/barycentric-coordinates.html
glm::vec3 calculateBarycentricCoords(const Triangle* triangle, const glm::vec3& point) {
	// Triangle Point1 in Cartesian Form
	glm::vec3 tP1Cartesian = glm::vec3(triangle->pointOne) / triangle->pointOne.w;
	glm::vec3 tP2Cartesian = glm::vec3(triangle->pointTwo) / triangle->pointTwo.w;
	glm::vec3 tP3Cartesian = glm::vec3(triangle->pointThree) / triangle->pointThree.w;

	// Calculate the vectors from pointOne to the other two vertices of the triangle
	// These vectors represent the edges of the triangle and the vector from the
	// first vertex to the point of interest. They are used for determining the relative
	// position of the point within the triangle.
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
	// std::cout << "1:" << u << "2:" << v << "3:" << w << std::endl;
	return glm::vec3(u, v, w);
}
// Function to compute the interpolated texture coordinate
glm::vec2 getTextureCoordinate(const glm::vec3& barycentricCoords, const glm::vec2& texCoordA, const glm::vec2& texCoordB, const glm::vec2& texCoordC) {
	// Interpolate the texture coordinates using the barycentric coordinates
	glm::vec2 texCoord = barycentricCoords.x * texCoordA
		+ barycentricCoords.y * texCoordB
		+ barycentricCoords.z * texCoordC;

	return texCoord;
}

// Phong Shading
glm::vec3 interpolateNormal(const Triangle& triangle, const glm::vec3& barycentricCoords) {
	// By multiplying each vertex normal by its corresponding barycentric coordinate, 
	// I am weighting each normal by how much influence that vertex has at the point of interest.
	return glm::normalize(
		barycentricCoords.x * triangle.normalOne +
		barycentricCoords.y * triangle.normalTwo +
		barycentricCoords.z * triangle.normalThree
	);
}

// Lighting with Phong Illumination Model
// Concept is found here: https://cg.informatik.uni-freiburg.de/course_notes/graphics_02_shading.pdf
glm::vec3 phongIllumination(const Triangle* triangle, const Ray* ray, const glm::vec3& lightPos, const glm::vec3& lightColor, const glm::vec3& objectColor, const float& ambientStrength, const float& specularStrength, const float& shininess, const float& distance) {
	// Phong illumination model
	// objectColor = object color
	// lightColor = Light Color
	// shininess = specular radius
	// specularStrength = specular Strength

	// rView is a constant factor for the reflection model, representing the light reflected to the view position
	// smaller number = less light reflected to view position = diffuse plays a smaller role
	constexpr float rView = 1.0f / glm::pi<float>();

	// Calculate the intersection point of the ray with the triangle
	glm::vec3 intersectionPoint = ray->origin + distance * ray->direction;

	// Calculate barycentric coordinates for the intersection point within the triangle
	glm::vec3 barycentricCoords = calculateBarycentricCoords(triangle, intersectionPoint);

	// Interpolate the normal at the intersection point using barycentric coordinates
	// glm::vec3 n = interpolateNormal(*triangle, barycentricCoords); // normal
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
	glm::vec3 v = glm::normalize(-ray->direction); // View Direction
	glm::vec3 r = glm::reflect(-l, n); // Reflect Direction

	// The specular term is calculated using the Phong reflection model
	// It is based on the dot product between the view direction and the reflection direction, raised to the power of the shininess factor (shininess)
	// specularStrength = specular strength. Smaller specular strength means less intensity
	glm::vec3 specular = lightColor * specularStrength * glm::max(dotProduct, 0.00f) * glm::pow(glm::max(glm::dot(r, v), 0.0f), shininess);

	// Combine the three components (diffuse, specular, and ambient) to get the final color
	return diffuse + specular + ambient;
}


bool intersectRayAabb(const glm::vec3& direction, const glm::vec3& minBox, const glm::vec3& maxBox) {
	// Calculate the t values for each pair of planes

	// origin is always at 0 because we are at view Space so we don't need to include it in 
	// the calculations

	// x-slab intersection t
	float minXT = minBox.x / direction.x;
	float maxXT = maxBox.x / direction.x;
	if (minXT > maxXT) {
		std::swap(minXT, maxXT);
	}

	float minYT = minBox.y / direction.y;
	float maxYT = maxBox.y / direction.y;
	if (minYT > maxYT) {
		std::swap(minYT, maxYT);
	}

	// Check if the intervals overlap
	if (maxXT < minYT || maxYT < minXT) {
		return false;
	}

	// Now do it for z axis
	if (minYT > minXT) {
		minXT = minYT;
	}
	if (maxYT < maxXT) {
		maxXT = maxYT;
	}

	float minZT = minBox.z / direction.z;
	float maxZT = maxBox.z / direction.z;

	if (minZT > maxZT) {
		std::swap(minZT, maxZT);
	}

	if ((minXT > maxZT) || (minZT > maxXT)) {
		return false;
	}

	return true;
}

bool intersectRayAabbNoOrigin(const Ray& ray, const glm::vec3& minBox, const glm::vec3& maxBox) {
	// Calculate the t values for each pair of planes

	// origin is always at 0 because we are at view Space so we don't need to include it in 
	// the calculations

	// x-slab intersection t
	float minXT = (minBox.x - ray.origin.x) / ray.direction.x;
	float maxXT = (maxBox.x - ray.origin.x) / ray.direction.x;
	if (minXT > maxXT) {
		std::swap(minXT, maxXT);
	}

	float minYT = (minBox.y - ray.origin.y) / ray.direction.y;
	float maxYT = (maxBox.y - ray.origin.y) / ray.direction.y;
	if (minYT > maxYT) {
		std::swap(minYT, maxYT);
	}

	// Check if the intervals overlap
	if (maxXT < minYT || maxYT < minXT) {
		return false;
	}

	// Now do it for z axis
	if (minYT > minXT) {
		minXT = minYT;
	}
	if (maxYT < maxXT) {
		maxXT = maxYT;
	}

	float minZT = (minBox.z - ray.origin.z) / ray.direction.z;
	float maxZT = (maxBox.z - ray.origin.z) / ray.direction.z;

	if (minZT > maxZT) {
		std::swap(minZT, maxZT);
	}

	if ((minXT > maxZT) || (minZT > maxXT)) {
		return false;
	}

	return true;
}
std::vector<Triangle> boundingBoxIntersection(Node* node, const Ray& ray) {
	// Only Collect Triangles if Ray and Box have intersection
	if (!intersectRayAabbNoOrigin(ray, node->minBox, node->maxBox)) {
		return {};
	}

	// If leaf node, just return its triangles.
	else if (!node->left && !node->right) {
		return node->triangles;
	}
	else {
		// Otherwise, collect all triangles from left and right.
		std::vector<Triangle> leftTriangles = boundingBoxIntersection(node->left, ray);
		std::vector<Triangle> rightTriangles = boundingBoxIntersection(node->right, ray);

		// Merge them.
		leftTriangles.insert(leftTriangles.end(),
			rightTriangles.begin(),
			rightTriangles.end());
		return leftTriangles;
	}
}

bool shadowIntersection(ObjectManager* objManager, const std::string& currentObjFilename, const glm::vec3& lightPos, const float& fDistance, const Ray& ray) {
	for (const auto& pairShadow : objManager->objTriangles) {
		const std::string& shadowObjFilename = pairShadow.first;

		/*glm::vec3 P = ray.origin + ray.direction * fDistance;
		glm::vec3 lightDir = glm::normalize(lightPos - P);

		Ray shadowRay(lightDir);
		shadowRay.origin = P; // + eps * normalAtP; // offset a bit above the surface
		// shadowRay.direction = lightDir; */

		// const std::vector<Triangle>& trianglesBox = pairShadow.second;
		Ray shadowRay(lightPos - ray.direction * fDistance);
		shadowRay.origin = ray.direction * fDistance;
		// shadowRay.origin += ray.direction * 0.001f; // add small value to prevent shadowAcne
		const std::vector<Triangle>& trianglesBox = boundingBoxIntersection(objManager->boundingVolumeHierarchy[shadowObjFilename], shadowRay);

		// if (intersectRayAabbNoOrigin(shadowRay, objManager.minBox[shadowObjFilename], objManager.maxBox[shadowObjFilename])) {
			// Prevents intersection between same object
			if (shadowObjFilename != currentObjFilename) {
				for (int j = 0; j < trianglesBox.size(); j++) {

					float shadowDistance = rayTriangleIntersection(&shadowRay, &trianglesBox[j]);
					if (shadowDistance != -INFINITY) {
						return true;
					// }
				}
			}
		}
	}
	return false;
}

// Generate multiple light Sources in near distance to each other
// And use tone mapping
glm::vec3 softShadow(int lightAmount,ObjectManager* objManager, const std::string& objFilename,const Triangle* triangle, const Ray* ray, const glm::vec3& lightPos, const glm::vec3& lightColor, const glm::vec3& objColor, const float& distance) {
	glm::vec3 testColor = objColor;
	if (!triangle->textureName.empty()) {
		glm::vec3 intersectionPoint = ray->origin + distance * ray->direction;
		glm::vec2 interpolatedTexCoordinate = getTextureCoordinate(calculateBarycentricCoords(triangle, intersectionPoint), triangle->colorOneCoordinate, triangle->colorTwoCoordinate, triangle->colorThreeCoordinate);
		// std::cout << "Interpolated X:" << interpolatedTexCoordinate.x << "Xone:" << triangle->colorOneCoordinate.x << "Xtwo:" << triangle->colorTwoCoordinate.x << "XThree:" << triangle->colorThreeCoordinate.x << std::endl;
		unsigned char* texData = objManager->textureData[triangle->textureName];
		glm::ivec2 texDim = objManager->textureDimensions[triangle->textureName];

		size_t texIndex = (static_cast<int>(interpolatedTexCoordinate.y) *texDim.x + static_cast<int>(interpolatedTexCoordinate.x)) * 3;
		testColor.x = texData[texIndex + 0] / 255.0f;
		testColor.y = texData[texIndex + 1] / 255.0f;
		testColor.z = texData[texIndex + 2] / 255.0f;
	}
	glm::vec3 color(0.0f,0.0f,0.0f);
	glm::vec3 lightPosChanged = lightPos;
	// Generate the lightPos based on the last lightPos and change always one coordinate
	for (int i = 0; i < lightAmount; i++) {
		bool isShadow = shadowIntersection(objManager, objFilename, lightPosChanged, distance, *ray);
		glm::vec3 colorPhong = phongIllumination(triangle, ray, lightPosChanged, lightColor, testColor, objManager->objProperties[objFilename][0], objManager->objProperties[objFilename][1], objManager->objProperties[objFilename][2], distance);
		if (isShadow) { colorPhong /= 5; };
		color += colorPhong;
		// Before it was 2
		switch (i % 3) {
		case 0:
			lightPosChanged.x += 3.0f;
			break;
		case 1:
			lightPosChanged.y += 3.0f; 
			break;
		case 2:
			lightPosChanged.z += 3.0f;
			break;
		}

	}

	// Reinhardt Tone Mapping 
	// color = color / (color + 1.0f);
	color = color / (color + 0.1f);
	// color = color / (color + 4.0f);

	// Add Gamma
	glm::vec3 gamma (2.2f, 2.2f, 2.2f);
	// color = glm::pow(color, gamma);

	return color;
}

std::pair<glm::vec2, glm::vec3> rayIntersection(const Ray& ray, ObjectManager* objManager, const int& pointX,const int& pointY, const glm::vec3& lightPos) {

	glm::vec3 colorPoint(0, 0, 0);
	float distanceComparison = INFINITY;
	for (const auto& pair : objManager->objTriangles) {
		const std::string& objFilename = pair.first;
		// std::cout << objFilename << pair.second.size();

		// const std::vector<Triangle>& triangles = pair.second;
		/*
		// Old Intersection for Speed Testing Purposes with OneBox and goes through all Triangles in an object
		if (intersectRayAabb(ray.direction, objManager.minBox[objFilename], objManager.maxBox[objFilename])) {

			for (int k = 0; k < triangles.size(); k++) {
			*/
		// std::cout << objFilename << std::endl;
		const std::vector<Triangle>& trianglesBox = boundingBoxIntersection(objManager->boundingVolumeHierarchy[objFilename], ray);

		for (int k = 0; k < trianglesBox.size(); k++) {

			float fDistance = rayTriangleIntersection(&ray, &trianglesBox[k]);

			if (fDistance != -INFINITY) {
				if (fDistance < distanceComparison) {

					distanceComparison = fDistance;

					glm::vec3 lightColor(1.0f, 1.0f, 1.0f); // White light

					// Code for many light Sources

					glm::vec3 objColor(1,0,0);
					if (trianglesBox[k].textureName.empty()) {
						objColor = objManager->objColors[objFilename];
					}
					else {
						objColor = trianglesBox[k].color;
					}
					// 36 Shadows are a good value
					glm::vec3 color = softShadow(1,objManager,objFilename, &trianglesBox[k], &ray,lightPos,lightColor, objColor,fDistance);
					// Convert 0...1 color values to 1...255 color Values
					colorPoint.x = int((color.x * 255));
					colorPoint.y = int((color.y * 255));
					colorPoint.z = int((color.z * 255));
					// }
				}
			}
		}
	}
	glm::vec2 imagePoint(pointX, pointY);
	return { imagePoint, colorPoint };
}

void drawImage(const glm::vec2& imgSize, const std::vector<glm::vec2>& imagePoints, const std::vector<glm::vec3>& imageColors, const int& angleDegree,const bool& saveImage, const bool& displayImage) {
	// create image
	CImg<unsigned char> img(imgSize.x, imgSize.y, 1, 3);
	img.fill(0);


	// draw pixels found in ray intersection
	for (int i = 0; i < imagePoints.size(); i++) {
		unsigned char color[] = { 0,0,0 };
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
	std::string imgName = "images/generation/output";
	imgName += std::to_string(static_cast<int>(angleDegree));
	imgName += ".bmp";

	if (saveImage == true) {
		img.save_bmp(imgName.c_str()); // Use c_str() to get a const char* from std::string
	}
	if (displayImage == true) {
		img.display("Simple Raytracer by Leon Lang");
	}
}



// Sending out Rays
// Sends out Rays and returns the corresponding color for each pixel
ImageData sendRaysAndIntersectPointsColors(const glm::vec2& imageSize, const glm::vec4& lightPos, ObjectManager* objManager) {
	Ray ray(glm::vec3(0.0f, 0.0f, 400.0f));
	glm::vec2 rayXY = glm::vec2(ray.direction.x, ray.direction.y);
	ImageData imageData;

	// std::cout << "Intersection: " << intersectRayAabb(glm::vec2(0, 0), glm::vec2(6, 3), glm::vec2(300, 300), glm::vec2(400, 400)) << " Intersects";
	for (int i = -imageSize.x / 2; i < imageSize.x / 2; ++i) {

		for (int j = -imageSize.y / 2; j < imageSize.y / 2; ++j) {
			ray.direction.x = i + rayXY.x;
			ray.direction.y = j + rayXY.y;

			std::pair<glm::vec2, glm::vec3> points = rayIntersection(ray, objManager, i + imageSize.x / 2, j + imageSize.y / 2, lightPos);
			if (points.second != glm::vec3(0, 0, 0)) {
				imageData.imagePoints.push_back(points.first);
				imageData.imageColors.push_back(points.second);
			}
		}	
	}
	return imageData;
}




int main()
{
	/*
	std::string filePath = "uyhh44.png"; // Image file path
	int x = 50; // X-coordinate
	int y = 100; // Y-coordinate

	getColorAtCoordinates(filePath, x, y);
	*/
	// save images at different degrees based on camera
	for (float angleDegree = 0; angleDegree < 360; angleDegree = angleDegree + 10) {

		// Create an ObjectManager instance 
		ObjectManager objManager;


		/*

		// create a complex scene with trees, bunnys and cats
		// Uses many features of the simple raytracer but takes time to output
		// Obj.z = left and right
		// Obj.y = up and down
		// Obj.x = 
		
		
		// Initiate Scene with Camera and View Matrix
		float radius = 50.0f; // Radius of the circle on which the camera moves
		float radians = glm::radians(angleDegree); // Convert angle from degrees to radians
		float circleX = radius * std::cos(radians); // Calculate x coordinate on the circle
		float circleZ = radius * std::sin(radians); // Calculate z coordinate on the circle
		glm::mat4 viewMatrix = Transformation::createViewMatrix(glm::vec3(circleX, -50.f, circleZ), glm::vec3(glm::radians(30.f), glm::radians(angleDegree + 90), glm::radians(0.f)));

		// 1 Cube Ground
		objManager.loadObjFile("./obj/cube.obj");
		objManager.setColor("./obj/cube.obj", glm::vec3(0.f, 1.f, 0.f));
		objManager.transformTriangles("./obj/cube.obj", Transformation::scaleObj(35.0f, 35.0f, 35.0f));
		objManager.transformTriangles("./obj/cube.obj", Transformation::changeObjPosition(glm::vec3(0.f, 10.f, 0.f)));
		objManager.transformTriangles("./obj/cube.obj", glm::inverse(viewMatrix));
		objManager.createBoundingHierarchy("./obj/cube.obj");
		
		
		// 2 Cats
		objManager.loadObjFile("./obj/cat/cat.obj");
		objManager.objProperties["./obj/cat/cat.obj"].y = 0.0f;
		objManager.objTriangles["./obj/cat/cat.obj1"] = objManager.getTriangles("./obj/cat/cat.obj");
		objManager.objProperties["./obj/cat/cat.obj1"] = objManager.objProperties["./obj/cat/cat.obj"];
		// Cat 0
		objManager.transformTriangles("./obj/cat/cat.obj", Transformation::scaleObj(0.35f, 0.35f, 0.35f));
		objManager.transformTriangles("./obj/cat/cat.obj", Transformation::rotateObjX(glm::radians(-90.f)));
		objManager.transformTriangles("./obj/cat/cat.obj", Transformation::rotateObjY(glm::radians(125.f)));
		objManager.transformTriangles("./obj/cat/cat.obj", Transformation::changeObjPosition(glm::vec3(25.f, -25.f, -14.f)));
		objManager.transformTriangles("./obj/cat/cat.obj", glm::inverse(viewMatrix));
		objManager.createBoundingHierarchy("./obj/cat/cat.obj");
		//Cat 1
		objManager.transformTriangles("./obj/cat/cat.obj1", Transformation::scaleObj(0.35f, 0.35f, 0.35f));
		objManager.transformTriangles("./obj/cat/cat.obj1", Transformation::rotateObjX(glm::radians(-90.f)));
		objManager.transformTriangles("./obj/cat/cat.obj1", Transformation::rotateObjY(glm::radians(70.f)));
		objManager.transformTriangles("./obj/cat/cat.obj1", Transformation::changeObjPosition(glm::vec3(25.f, -25.f, 8.f)));
		objManager.transformTriangles("./obj/cat/cat.obj1", glm::inverse(viewMatrix));
		objManager.createBoundingHierarchy("./obj/cat/cat.obj1");
		
		// 1 Stanford Bunny
		objManager.loadObjFile("./obj/stanford-bunny.obj");
		objManager.transformTriangles("./obj/stanford-bunny.obj", Transformation::scaleObj(50.f, 50.0f, 50.0f));
		objManager.transformTriangles("./obj/stanford-bunny.obj", Transformation::rotateObjX(glm::radians(181.f)));
		objManager.transformTriangles("./obj/stanford-bunny.obj", Transformation::rotateObjY(glm::radians(90.f)));

		objManager.transformTriangles("./obj/stanford-bunny.obj", Transformation::changeObjPosition(glm::vec3(25.f, -23.f, 0.f)));
		objManager.transformTriangles("./obj/stanford-bunny.obj", glm::inverse(viewMatrix));
		objManager.createBoundingHierarchy("./obj/stanford-bunny.obj");
		
		// 3 Trees
		objManager.loadObjFile("./obj/tree/tree.obj");
		objManager.objProperties["./obj/tree/tree.obj"].y = 0.0f;
		// Copy Tree onto Tree1 and 2
		objManager.objTriangles["./obj/tree/tree.obj1"] = objManager.getTriangles("./obj/tree/tree.obj");
		objManager.objProperties["./obj/tree/tree.obj1"] = objManager.objProperties["./obj/tree/tree.obj"];
		objManager.objTriangles["./obj/tree/tree.obj2"] = objManager.getTriangles("./obj/tree/tree.obj");
		objManager.objProperties["./obj/tree/tree.obj2"] = objManager.objProperties["./obj/tree/tree.obj"];
		// Place Tree 0
		objManager.transformTriangles("./obj/tree/tree.obj", Transformation::scaleObj(0.03f, 0.03f, 0.03f));
		objManager.transformTriangles("./obj/tree/tree.obj", Transformation::rotateObjX(glm::radians(-90.f)));
		objManager.transformTriangles("./obj/tree/tree.obj", Transformation::changeObjPosition(glm::vec3(-6.f, -25.f, -25.f)));
		objManager.transformTriangles("./obj/tree/tree.obj", glm::inverse(viewMatrix));
		objManager.createBoundingHierarchy("./obj/tree/tree.obj");
		// Place Tree 1
		objManager.transformTriangles("./obj/tree/tree.obj1", Transformation::scaleObj(0.035f, 0.035f, 0.035f));
		objManager.transformTriangles("./obj/tree/tree.obj1", Transformation::rotateObjX(glm::radians(-90.f)));
		objManager.transformTriangles("./obj/tree/tree.obj1", Transformation::changeObjPosition(glm::vec3(-6.f, -25.f, 0.f)));
		objManager.transformTriangles("./obj/tree/tree.obj1", glm::inverse(viewMatrix));
		objManager.createBoundingHierarchy("./obj/tree/tree.obj1");
		// Place Tree 2
		objManager.transformTriangles("./obj/tree/tree.obj2", Transformation::scaleObj(0.03f, 0.03f, 0.03f));
		objManager.transformTriangles("./obj/tree/tree.obj2", Transformation::rotateObjX(glm::radians(-90.f)));
		objManager.transformTriangles("./obj/tree/tree.obj2", Transformation::changeObjPosition(glm::vec3(-6.f, -25.f, 25.f)));
		objManager.transformTriangles("./obj/tree/tree.obj2", glm::inverse(viewMatrix));
		objManager.createBoundingHierarchy("./obj/tree/tree.obj2");
		*/


		/*
		// Scene: 6 Sphere Triangles transformed
		 
		// Inititate Scene
		float radius = 100.0f; // Radius of the circle on which the camera moves
		float radians = glm::radians(angleDegree); // Convert angle from degrees to radians
		float circleX = radius * std::cos(radians); // Calculate x coordinate on the circle
		float circleZ = radius * std::sin(radians); // Calculate z coordinate on the circle

		// If you want to have interpolated spheres:
		// Please uncomment this line in Phong Illumination function: 
		// glm::vec3 n = interpolateNormal(*triangle, barycentricCoords); // normal
		// And comment this line:
		// glm::vec3 n = calculateTriangleNormal(triangle);

		objManager.loadObjFile("sphere.obj");
		objManager.transformTriangles("sphere.obj", Transformation::changeObjPosition(glm::vec3(0.f, 6.f, 30.f)));

		// Load 5 more Spheres for time Testing
		objManager.objTriangles["sphere1.obj"] = objManager.getTriangles("sphere.obj");
		objManager.objColors["sphere1.obj"] = glm::vec3(1.f, 0.f, 0.f);
		objManager.transformTriangles("sphere1.obj", Transformation::changeObjPosition(glm::vec3(6.f, 0.f, 0.f)));

		// Load 8 more Spheres for time Testing
		objManager.objTriangles["sphere2.obj"] = objManager.getTriangles("sphere.obj");
		objManager.objColors["sphere2.obj"] = glm::vec3(1.f, 0.f, 0.f);
		objManager.transformTriangles("sphere2.obj", Transformation::changeObjPosition(glm::vec3(-6.f, 0.f, 0.f)));

		// Load 8 more Spheres for time Testing
		objManager.objTriangles["sphere3.obj"] = objManager.getTriangles("sphere.obj");
		objManager.objColors["sphere3.obj"] = glm::vec3(1.f, 0.f, 0.f);
		objManager.transformTriangles("sphere3.obj", Transformation::changeObjPosition(glm::vec3(0.f, -12.f, 0.f)));

		// Load 8 more Spheres for time Testing
		objManager.objTriangles["sphere4.obj"] = objManager.getTriangles("sphere.obj");
		objManager.objColors["sphere4.obj"] = glm::vec3(1.f, 0.f, 0.f);
		objManager.transformTriangles("sphere4.obj", Transformation::changeObjPosition(glm::vec3(6.f, -12.f, 0.f)));

		// Load 8 more Spheres for time Testing
		objManager.objTriangles["sphere5.obj"] = objManager.getTriangles("sphere.obj");
		objManager.objColors["sphere5.obj"] = glm::vec3(1.f, 0.f, 0.f);
		objManager.transformTriangles("sphere5.obj", Transformation::changeObjPosition(glm::vec3(-6.f, -12.f, 0.f)));
		objManager.createBoundingHierarchy("sphere.obj");
		objManager.createBoundingHierarchy("sphere1.obj");
		objManager.createBoundingHierarchy("sphere2.obj");
		objManager.createBoundingHierarchy("sphere3.obj");
		objManager.createBoundingHierarchy("sphere4.obj");
		objManager.createBoundingHierarchy("sphere5.obj");
		*/


		
		/*
		// Big Cat Sideways for testing purposes

		// Inititate Scene with view Matrix
		float radius = 50.0f; // Radius of the circle on which the camera moves
		float radians = glm::radians(angleDegree); // Convert angle from degrees to radians
		float circleX = radius * std::cos(radians); // Calculate x coordinate on the circle
		float circleZ = radius * std::sin(radians); // Calculate z coordinate on the circle
		glm::mat4 viewMatrix = Transformation::createViewMatrix(glm::vec3(circleX, -50.f, circleZ), glm::vec3(glm::radians(30.f), glm::radians(angleDegree + 90), glm::radians(0.f)));

		objManager.loadObjFile("./obj/cat/cat.obj");
		objManager.objProperties["./obj/cat/cat.obj"].y = 0.0f;
		objManager.transformTriangles("./obj/cat/cat.obj", Transformation::scaleObj(0.5f, 0.5f, 0.5f));
		objManager.transformTriangles("./obj/cat/cat.obj", Transformation::rotateObjX(glm::radians(-90.f)));
		// objManager.transformTriangles("./obj/cat/cat.obj", Transformation::rotateObjY(glm::radians(90.f)));
		// objManager.transformTriangles("chair.obj", Transformation::rotateObjZ(glm::radians(90.f)));
		objManager.transformTriangles("./obj/cat/cat.obj", Transformation::changeObjPosition(glm::vec3(14.f, -20.f, -8.f)));
		objManager.transformTriangles("./obj/cat/cat.obj", glm::inverse(viewMatrix));
		objManager.createBoundingHierarchy("./obj/cat/cat.obj");
		*/
		
		



		
		// One Sample Cube for Testing
		
		// Inititate Scene with view Matrix
		float radius = 100.0f; // Radius of the circle on which the camera moves
		float radians = glm::radians(angleDegree); // Convert angle from degrees to radians
		float circleX = radius * std::cos(radians); // Calculate x coordinate on the circle
		float circleZ = radius * std::sin(radians); // Calculate z coordinate on the circle
		glm::mat4 viewMatrix = Transformation::createViewMatrix(glm::vec3(circleX, 0.f, circleZ), glm::vec3(glm::radians(0.f), glm::radians(angleDegree + 90), glm::radians(0.f)));

		
		objManager.loadObjFile("cube.obj");
		// objManager.setColor("cube.obj", glm::vec3(1.f, 1.f, 0.f));
		objManager.transformTriangles("cube.obj", Transformation::scaleObj(20.0f, 20.0f, 20.0f));
		objManager.transformTriangles("cube.obj", Transformation::rotateObjY(glm::radians(25.f)));
		// objManager.transformTriangles("cube.obj", Transformation::mirrorObj(true, false, false));
		// objManager.transformTriangles("cube.obj", Transformation::changeObjPosition(glm::vec3(0.f, 10.f, 0.f)));
		objManager.transformTriangles("cube.obj", glm::inverse(viewMatrix));
		objManager.createBoundingHierarchy("cube.obj");
		



		/*
		// Scene with 4 Cubes in different colors

		// Inititate Scene with view Matrix
		float radius = 100.0f; // Radius of the circle on which the camera moves
		float radians = glm::radians(angleDegree); // Convert angle from degrees to radians
		float circleX = radius * std::cos(radians); // Calculate x coordinate on the circle
		float circleZ = radius * std::sin(radians); // Calculate z coordinate on the circle
		glm::mat4 viewMatrix = Transformation::createViewMatrix(glm::vec3(circleX, 0.f, circleZ), glm::vec3(glm::radians(0.f), glm::radians(angleDegree + 90), glm::radians(0.f)));

		// Load Cube Triangles and scale it
		objManager.loadObjFile("cube.obj");
		objManager.setColor("cube.obj", glm::vec3(1.f, 1.f, 0.f));
		objManager.transformTriangles("cube.obj", Transformation::scaleObj(10.0f, 10.0f, 10.0f));

		// Create Cube Triangles Clone with different Position and Colors
		objManager.objTriangles["cube1.obj"] = objManager.getTriangles("cube.obj");
		objManager.objColors["cube1.obj"] = glm::vec3(1.f,0.f,1.f);
		objManager.transformTriangles("cube1.obj", Transformation::changeObjPosition(glm::vec3(0.f, -15.f, -15.f)));

		objManager.objTriangles["cube2.obj"] = objManager.getTriangles("cube.obj");
		objManager.objColors["cube2.obj"] = glm::vec3(1.f, 0.f, 0.f);
		objManager.transformTriangles("cube2.obj", Transformation::changeObjPosition(glm::vec3(0.f, -15.f, 15.f)));

		objManager.objTriangles["cube3.obj"] = objManager.getTriangles("cube.obj");
		objManager.objColors["cube3.obj"] = glm::vec3(0.f, 1.f, 0.f);
		objManager.transformTriangles("cube3.obj", Transformation::changeObjPosition(glm::vec3(0.f, 15.f, 15.f)));

		// transform position of original cube
		objManager.transformTriangles("cube.obj", Transformation::changeObjPosition(glm::vec3(0.f, 15.f, -15.f)));
		
		// Bring Obj Models into ViewPosition
		objManager.transformTriangles("cube.obj", glm::inverse(viewMatrix));
		objManager.transformTriangles("cube1.obj", glm::inverse(viewMatrix));
		objManager.transformTriangles("cube2.obj", glm::inverse(viewMatrix));
		objManager.transformTriangles("cube3.obj", glm::inverse(viewMatrix));
		
		// Create Bounding Volume Hierarchy for Cubes
		objManager.createBoundingHierarchy("cube.obj");
		objManager.createBoundingHierarchy("cube1.obj");
		objManager.createBoundingHierarchy("cube2.obj");
		objManager.createBoundingHierarchy("cube3.obj");
		*/
		

		// Draw Image
		glm::vec2 imageSize(600, 400);
		//glm::vec4 lightPos(-200.0f, -300.0f, -1000.4f, 1.0f);
		
		glm::vec4 lightPos(500.0f, -300.0f, -200.f, 1.0f);
		
		lightPos = glm::inverse(viewMatrix) * lightPos;
		// lightPos.z = -lightPos.z;
		// lightPos = lightPos * 20001.f;
		// Start the timer 
		auto start = std::chrono::high_resolution_clock::now();

		ImageData points = sendRaysAndIntersectPointsColors(imageSize, lightPos, &objManager);


		// End the timer 
		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = end - start;
		// Print the time taken 
		std::cout << "Time taken for Intersection: " << elapsed.count() << " seconds " << std::endl;

		drawImage(imageSize, points.imagePoints, points.imageColors, angleDegree, true, true);

	}
}
