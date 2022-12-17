#include <iostream>
#include <glm/glm.hpp>              

using Point = glm::vec3;          
using Point2D = glm::vec2;     

Point homogeneous_coords(const Point2D &p)
{
    return Point(p, 1.0);
}

Point2D affinize(const Point &p)
{
    return Point2D(glm::round(p[0] / p[2]), glm::round(p[1] / p[2]));
}

Point average(const Point &P1, const Point &P2, const Point &P3)
{
    Point T;
    for (unsigned i = 0; i < 3; ++i)
        T[i] = (P1[i] + P2[i] + P3[i]) / 3.0;
    return T;
}

Point2D invisible_point(const Point2D &P1, const Point2D &P2, const Point2D &P3, const Point2D &P5, const Point2D &P6, const Point2D &P7, const Point2D &P8)
{
    Point T1 = homogeneous_coords(P1);
    Point T2 = homogeneous_coords(P2);
    Point T3 = homogeneous_coords(P3);
    Point T5 = homogeneous_coords(P5);
    Point T6 = homogeneous_coords(P6);
    Point T7 = homogeneous_coords(P7);
    Point T8 = homogeneous_coords(P8);

    // Xinf
    Point T2T6 = glm::cross(T2, T6);
    Point T1T5 = glm::cross(T1, T5);
    Point T3T7 = glm::cross(T3, T7);

    Point X1_inf = glm::cross(T2T6, T1T5);
    Point X2_inf = glm::cross(T2T6, T3T7);
    Point X3_inf = glm::cross(T1T5, T3T7);
    Point X_inf  = average(X1_inf, X2_inf, X3_inf);
    
    // Yinf
    Point T5T6 = glm::cross(T5, T6);
    Point T7T8 = glm::cross(T7, T8);
    Point T1T2 = glm::cross(T1, T2);

    Point Y1_inf = glm::cross(T5T6, T7T8);
    Point Y2_inf = glm::cross(T1T2, T7T8);
    Point Y3_inf = glm::cross(T1T2, T5T6);
    Point Y_inf  = average(Y1_inf, Y2_inf, Y3_inf);
    
    Point P = glm::cross(T8, X_inf);
    Point Q = glm::cross(T3, Y_inf);

    Point P4 = glm::cross(P, Q);
    return affinize(P4);
}

int main ()
{
    // Point2D P1 = Point2D(595, 301);
    // Point2D P2 = Point2D(292, 517);
    // Point2D P3 = Point2D(157, 379);
    // Point2D P5 = Point2D(665, 116);
    // Point2D P6 = Point2D(304, 295);
    // Point2D P7 = Point2D(135, 163);
    // Point2D P8 = Point2D(509, 43);

    Point2D P1 = Point2D(870, 333);
    Point2D P2 = Point2D(305, 629);
    Point2D P3 = Point2D(111, 325);
    Point2D P5 = Point2D(892, 272);
    Point2D P6 = Point2D(300, 551);
    Point2D P7 = Point2D(97, 267);
    Point2D P8 = Point2D(610, 148);    

    Point2D P4 = invisible_point(P1, P2, P3, P5, P6, P7, P8);
    std::cout << "Invisible point coordinates: " << "(" << P4[0] << ", " << P4[1] << ")" << std::endl;

    return 0;
}
