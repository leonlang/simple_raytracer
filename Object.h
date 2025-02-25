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
#include <algorithm>


class Triangle {
public:
    glm::vec4 pointOne;
    glm::vec4 pointTwo;
    glm::vec4 pointThree;
    glm::vec3 normalOne;  // Normal at pointOne
    glm::vec3 normalTwo;  // Normal at pointTwo
    glm::vec3 normalThree; // Normal at pointThree
    glm::vec2 colorOneCoordinate;
    glm::vec2 colorTwoCoordinate;
    glm::vec2 colorThreeCoordinate;
    glm::vec3 color;
    std::string textureName;
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

struct Node {
    glm::vec3 minBox;
    glm::vec3 maxBox;
    std::vector<Triangle> triangles;
    Node* left;
    Node* right;

    // Constructor to initialize the node with a value
    Node(const std::vector<Triangle>& triangleD, const glm::vec3& minBoxD, const glm::vec3& maxBoxD)
        : triangles(triangleD), minBox(minBoxD), maxBox(maxBoxD), left(nullptr), right(nullptr) {
    }
};

class ObjectManager {
public:
    std::unordered_map<std::string, glm::vec3> minBox;
    std::unordered_map<std::string, glm::vec3> maxBox;
    // ambientStrength, specularStrength and shininess
    std::unordered_map<std::string, glm::vec3> objProperties;
    std::unordered_map<std::string, std::vector<Triangle>> objTriangles;
    std::unordered_map<std::string, Node*> boundingVolumeHierarchy;
    // define object color
    std::unordered_map<std::string, glm::vec3> objColors;
    // loaded textures
    std::unordered_map<std::string, unsigned char*> textureData;
    std::unordered_map<std::string, glm::ivec2> textureDimensions;
    // Method to load triangles from an OBJ file
    void loadObjFile(const std::string& objFilename);

    // Method to get triangles for a specific OBJ file
    const std::vector<Triangle>& getTriangles(const std::string& objFilename) const;

    // Method to set triangles for a specific OBJ file 
    void setTriangles(const std::string& objFilename, const std::vector<Triangle>& triangles);

    // Method to transform triangles for a specific OBJ file 
    void transformTriangles(const std::string& objFilename, const glm::mat4& matrix);
    void splitTrianglesForBox(Node* root);
    void createBoundingHierarchy(const std::string& objFilename);
    // Method to set color for a specific OBJ file 
    void setColor(const std::string& objFilename, const glm::vec3& color);
    // Method to get color for a specific OBJ file 
    glm::vec3 getColor(const std::string& objFilename) const;
};

#endif // OBJECT_H
