#include "Transformation.h"

glm::mat4 Transformation::viewSpaceTransformation(float angleDegree) {
    float radius = 80.0f;
    float radians = glm::radians(angleDegree);

    float circleX = radius * std::cos(radians);
    float circleZ = radius * std::sin(radians);

    glm::mat4 viewMatrix(1.0f);
    glm::vec3 cameraPosition(0.f, 0.f, 100.f);
    glm::vec3 targetPosition(0.0f, 0.f, 0.0f);
    glm::vec3 upVector(0.0f, 1.0f, 0.0f);

    glm::vec3 zAxis = glm::normalize(cameraPosition - targetPosition);
    glm::vec3 xAxis = glm::normalize(glm::cross(upVector, zAxis));
    glm::vec3 yAxis = glm::cross(zAxis, xAxis);

    viewMatrix[0][0] = xAxis.x;
    viewMatrix[1][0] = xAxis.y;
    viewMatrix[2][0] = xAxis.z;
    viewMatrix[0][1] = yAxis.x;
    viewMatrix[1][1] = yAxis.y;
    viewMatrix[2][1] = yAxis.z;
    viewMatrix[0][2] = zAxis.x;
    viewMatrix[1][2] = zAxis.y;
    viewMatrix[2][2] = zAxis.z;

    viewMatrix[3][0] = -glm::dot(xAxis, cameraPosition);
    viewMatrix[3][1] = -glm::dot(yAxis, cameraPosition);
    viewMatrix[3][2] = -glm::dot(zAxis, cameraPosition);

    return viewMatrix;
}

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

glm::mat4 Transformation::createViewMatrix(glm::vec3 position, glm::vec3 rotation) {
    glm::mat4 matrix = changeObjPosition(position);
    matrix *= rotateObjZ(rotation.z);
    matrix *= rotateObjY(rotation.y);
    matrix *= rotateObjX(rotation.x);
    return matrix;
}