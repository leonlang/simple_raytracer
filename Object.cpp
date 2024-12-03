#include "Object.h"
#include <iostream>

// Default constructor implementation 
Triangle::Triangle() 
    : pointOne(0.0f), pointTwo(0.0f), pointThree(0.0f), normalOne(0.0f), normalTwo(0.0f), normalThree(0.0f) {}

// Parameterized constructor implementation
Triangle::Triangle(glm::vec4 p1, glm::vec4 p2, glm::vec4 p3, glm::vec3 n1, glm::vec3 n2, glm::vec3 n3)
    : pointOne(p1), pointTwo(p2), pointThree(p3), normalOne(n1), normalTwo(n2), normalThree(n3) {}

glm::vec3 Triangle::calculateNormal() const {
    glm::vec3 v1 = glm::vec3(pointTwo) - glm::vec3(pointOne);
    glm::vec3 v2 = glm::vec3(pointThree) - glm::vec3(pointOne);
    glm::vec3 normal = glm::cross(v1, v2);
    return glm::normalize(normal);
}

// Ray class implementation 
Ray::Ray(glm::vec3 d) 
    : origin(0.0f, 0.0f, 0.0f), direction(d) {}


// ObjectManager class implementation
void ObjectManager::loadObjFile(const std::string& objFilename) {
    std::string inputFile = objFilename;
    tinyobj::ObjReader reader;

    objColors[objFilename] = glm::vec3(1.f,0.f,1.f);

    if (!reader.ParseFromFile(inputFile)) {
        if (!reader.Error().empty()) {
            std::cerr << "TinyObjReader: " << reader.Error();
        }
    }

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
                Triangle triangle(
                    glm::vec4(vertices[0], 1),
                    glm::vec4(vertices[1], 1),
                    glm::vec4(vertices[2], 1),
                    normals[0],
                    normals[1],
                    normals[2]
                );
                triangles.push_back(triangle);
                vertices.clear();
                normals.clear();
            }
        }
    }

    objTriangles[objFilename] = triangles;
}

const std::vector<Triangle>& ObjectManager::getTriangles(const std::string& objFilename) const {
    return objTriangles.at(objFilename);
}
void ObjectManager::setTriangles(const std::string& objFilename, const std::vector<Triangle>& triangles) {
    objTriangles[objFilename] = triangles; 
}

void ObjectManager::transformTriangles(const std::string& objFilename, const glm::mat4& matrix) { 
    std::vector<Triangle>& triangles = objTriangles[objFilename]; 
    for (Triangle& t : triangles) { 
        t.pointOne = matrix * t.pointOne;
        t.pointTwo = matrix * t.pointTwo;
        t.pointThree = matrix * t.pointThree; 
    } 
}

void ObjectManager::setColor(const std::string& objFilename, const glm::vec3& color) {
    objColors[objFilename] = color; 
} 
glm::vec3 ObjectManager::getColor(const std::string& objFilename) const {
    return objColors.at(objFilename); 
}