#ifndef TRANSFORMATION_H
#define TRANSFORMATION_H

#include <glm.hpp>
#include <gtc/matrix_transform.hpp>
#include <gtc/type_ptr.hpp>
#include <iostream>


class Transformation {
public:
    static glm::mat4 scaleObj(float sx, float sy, float sz);
    static glm::mat4 rotateObjX(float degree);
    static glm::mat4 rotateObjY(float degree);
    static glm::mat4 rotateObjZ(float degree);
    static glm::mat4 mirrorObj(bool mirrorX = false, bool mirrorY = false, bool mirrorZ = false);
    static glm::mat4 shearObj(float shearXY = 0.0f, float shearXZ = 0.0f, float shearYX = 0.0f, float shearYZ = 0.0f, float shearZX = 0.0f, float shearZY = 0.0f);
    static glm::mat4 changeObjPosition(glm::vec3 position);
    static glm::mat4 createViewMatrix(glm::vec3 position, glm::vec3 rotation);
};

#endif // TRANSFORMATION_H
