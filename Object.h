#ifndef OBJECT_H
#define OBJECT_H

#include <iostream>
#include <glm.hpp>
#include <gtc/matrix_transform.hpp>
#include <gtc/type_ptr.hpp>
#include <vector>
#include <string>
#include <unordered_map>
#include "tiny_obj_loader.h"

class Triangle {
public:
    glm::vec4 pointOne;
    glm::vec4 pointTwo;
    glm::vec4 pointThree;
    glm::vec3 normalOne;  // Normal at pointOne
    glm::vec3 normalTwo;  // Normal at pointTwo
    glm::vec3 normalThree; // Normal at pointThree

    // Default constructor
    Triangle();

    // Parameterized constructor 
    Triangle(glm::vec4 p1, glm::vec4 p2, glm::vec4 p3, glm::vec3 n1, glm::vec3 n2, glm::vec3 n3);

    // Method to calculate the normal of the triangle
    glm::vec3 calculateNormal() const;
};

class Ray {
public:
    glm::vec3 origin;
    glm::vec3 direction;

    // Constructor
    Ray(glm::vec3 d);
};

class ObjectManager {
public:
    std::unordered_map<std::string, std::vector<Triangle>> objTriangles;

    // Method to load triangles from an OBJ file
    void loadObjFile(const std::string& objFilename);

    // Method to get triangles for a specific OBJ file
    const std::vector<Triangle>& getTriangles(const std::string& objFilename) const;

    // Method to set triangles for a specific OBJ file 
    void setTriangles(const std::string& objFilename, const std::vector<Triangle>& triangles);

    // Method to transform triangles for a specific OBJ file 
    void transformTriangles(const std::string& objFilename, const glm::mat4& matrix);
};

#endif // OBJECT_H
