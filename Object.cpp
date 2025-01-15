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

    objColors[objFilename] = glm::vec3(1.f,0.f,0.f);

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
            /*
            glm::vec3 normal(
                attrib.normals[3 * index.normal_index + 0],
                attrib.normals[3 * index.normal_index + 1],
                attrib.normals[3 * index.normal_index + 2]
            ); */
            static std::vector<glm::vec3> vertices;
            vertices.push_back(vertex);
            static std::vector<glm::vec3> normals;
            // normals.push_back(normal);
            if (vertices.size() == 3) {
                /*
                Triangle triangle(
                    glm::vec4(vertices[0], 1),
                    glm::vec4(vertices[1], 1),
                    glm::vec4(vertices[2], 1),
                    normals[0],
                    normals[1],
                    normals[2]
                );*/
                Triangle triangle(
                    glm::vec4(vertices[0], 1),
                    glm::vec4(vertices[1], 1),
                    glm::vec4(vertices[2], 1),
                    glm::vec3(0,0,0),
                    glm::vec3(0, 0, 0),
                    glm::vec3(0, 0, 0)
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
bool compareXPointsOfTriangle(const Triangle& a, const Triangle& b) {
    return a.pointOne.x < b.pointOne.x;
}
bool compareYPointsOfTriangle(const Triangle& a, const Triangle& b) {
    return a.pointOne.y < b.pointOne.y;
}
bool compareZPointsOfTriangle(const Triangle& a, const Triangle& b) {
    return a.pointOne.z < b.pointOne.z;
}

std::pair<glm::vec3,glm::vec3> calculateBoundingBoxes(const std::vector<Triangle>& triangles) {

    glm::vec3 minBox = glm::vec3(FLT_MAX, FLT_MAX, FLT_MAX);;
    glm::vec3 maxBox = glm::vec3(-FLT_MAX, -FLT_MAX, -FLT_MAX);

    for (const Triangle& t : triangles) {
        minBox = glm::min(minBox, glm::vec3(t.pointOne));
        minBox = glm::min(minBox, glm::vec3(t.pointTwo));
        minBox = glm::min(minBox, glm::vec3(t.pointThree));

        maxBox = glm::max(maxBox, glm::vec3(t.pointOne));
        maxBox = glm::max(maxBox, glm::vec3(t.pointTwo));
        maxBox = glm::max(maxBox, glm::vec3(t.pointThree));
    }

    return { minBox, maxBox };
}


void ObjectManager::splitTrianglesForBox(Node* root) {
    // Sort triangles: The largest Side of the Box is sorted so I can cut this side of the box in Half
    // (Triangles get put into right or left side of Tree)
    float boxXSize = fabs(root->maxBox.x - root->minBox.x);
    float boxYSize = fabs(root->maxBox.y - root->minBox.y);
    float boxZSize = fabs(root->maxBox.z - root->minBox.z);

    glm::vec3 minBoxLeft;
    glm::vec3 maxBoxLeft;
    glm::vec3 minBoxRight;
    glm::vec3 maxBoxRight;
    std::vector<Triangle> trianglesLeftSide;
    std::vector<Triangle> trianglesRightSide;
    std::size_t triangleSize = root->triangles.size();

    if (boxXSize > boxYSize && boxXSize > boxZSize) {
        std::sort(root->triangles.begin(), root->triangles.end(), compareXPointsOfTriangle);
    }
    else if (boxYSize > boxXSize && boxYSize > boxZSize) {
        std::sort(root->triangles.begin(), root->triangles.end(), compareYPointsOfTriangle);
    }
    else{
        std::sort(root->triangles.begin(), root->triangles.end(), compareZPointsOfTriangle);
    }

    trianglesLeftSide.insert(trianglesLeftSide.end(), root->triangles.begin(), root->triangles.begin() + triangleSize/2);
    trianglesRightSide.insert(trianglesRightSide.end(), root->triangles.begin() + triangleSize / 2, root->triangles.end());
    std::tie(minBoxLeft, maxBoxLeft) = calculateBoundingBoxes(trianglesLeftSide);
    std::tie(minBoxRight, maxBoxRight) = calculateBoundingBoxes(trianglesRightSide);

    root->left = new Node(trianglesLeftSide, minBoxLeft, maxBoxLeft);
    root->right = new Node(trianglesRightSide, minBoxRight, maxBoxRight);
    // Build Binary Tree Recursively until size is <= triangleSizeStop
    int triangleSizeStop = 8;
    if (trianglesLeftSide.size() > triangleSizeStop) {
        splitTrianglesForBox(root->left);
    }
    if (trianglesRightSide.size() > triangleSizeStop) {
        splitTrianglesForBox(root->right);
    }


    // std::cout << "Triangles" << triangleSize << "Left" << trianglesLeftSide.size() << "Right" << trianglesRightSide.size() << std::endl;

}
void  ObjectManager::createBoundingHierarchy(const std::string& objFilename) {

    std::vector<Triangle>& triangles = objTriangles[objFilename];

    // Assign First Min and Max Box to Triangles
    std::tie(minBox[objFilename], maxBox[objFilename]) = calculateBoundingBoxes(triangles);
    Node* root = new Node(triangles, minBox[objFilename], maxBox[objFilename]);
    splitTrianglesForBox(root);
    boundingVolumeHierarchy[objFilename] = root;
}

void ObjectManager::setColor(const std::string& objFilename, const glm::vec3& color) {
    objColors[objFilename] = color; 
} 
glm::vec3 ObjectManager::getColor(const std::string& objFilename) const {
    return objColors.at(objFilename); 
}