#include "Transformation.h"

// Homogeneous Notation and ModelView Transform
// Concept for the Transformations and for ModelView: https://cg.informatik.uni-freiburg.de/course_notes/graphics_03_homogeneousNotation.pdf

glm::mat4 Transformation::scaleObj(float sx, float sy, float sz) {
    glm::mat4 matrix = glm::mat4(0.0f);
    matrix[0][0] = sx;
    matrix[1][1] = sy;
    matrix[2][2] = sz;
    matrix[3][3] = 1.0f;
    return matrix;
}

glm::mat4 Transformation::rotateObjX(float degree) {
    glm::mat4 matrix = glm::mat4(0.0f);
    matrix[0][0] = 1;
    matrix[1][1] = std::cos(degree);
    matrix[1][2] = -std::sin(degree);
    matrix[2][1] = std::sin(degree);
    matrix[2][1] = std::sin(degree);
    matrix[2][2] = std::cos(degree);
    matrix[3][3] = 1.0f;
    return matrix;
}

glm::mat4 Transformation::rotateObjY(float degree) {
    glm::mat4 matrix = glm::mat4(0.0f);
    matrix[0][0] = std::cos(degree);
    matrix[0][2] = std::sin(degree);
    matrix[1][1] = 1.0f;
    matrix[2][0] = -std::sin(degree);
    matrix[2][2] = std::cos(degree);
    matrix[3][3] = 1.0f;
    return matrix;
}

glm::mat4 Transformation::rotateObjZ(float degree) {
    glm::mat4 matrix = glm::mat4(0.0f);
    matrix[0][0] = std::cos(degree);
    matrix[0][1] = -std::sin(degree);
    matrix[1][0] = std::sin(degree);
    matrix[1][1] = std::cos(degree);
    matrix[2][2] = 1.0f;
    matrix[3][3] = 1.0f;
    return matrix;
}

glm::mat4 Transformation::mirrorObj(bool mirrorX, bool mirrorY, bool mirrorZ) {
    glm::mat4 matrix = glm::mat4(1.0f);

    if (mirrorX) {
        matrix[0][0] = -1.0f;
    }
    if (mirrorY) {
        matrix[1][1] = -1.0f;
    }
    if (mirrorZ) {
        matrix[2][2] = -1.0f;
    }

    return matrix;
}

glm::mat4 Transformation::shearObj(float shearXY, float shearXZ, float shearYX, float shearYZ, float shearZX, float shearZY) {
    glm::mat4 matrix = glm::mat4(1.0f);

    matrix[1][0] = shearXY;
    matrix[2][0] = shearXZ;
    matrix[0][1] = shearYX;
    matrix[2][1] = shearYZ;
    matrix[0][2] = shearZX;
    matrix[1][2] = shearZY;

    return matrix;
}

glm::mat4 Transformation::changeObjPosition(glm::vec3 position) {
    glm::mat4 translationMatrix = glm::mat4(1.0f);
    translationMatrix[3] = glm::vec4(position, 1.0f);
    return translationMatrix;
}
// Model View Transform
glm::mat4 Transformation::createViewMatrix(glm::vec3 position, glm::vec3 rotation) {
    glm::mat4 matrix = changeObjPosition(position);
    matrix *= rotateObjZ(rotation.z);
    matrix *= rotateObjY(rotation.y);
    matrix *= rotateObjX(rotation.x);
    return matrix;
}
