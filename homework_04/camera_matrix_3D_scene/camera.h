#ifndef CAMERA_CAMERA_H
#define CAMERA_CAMERA_H

#include <iostream>

#include <Eigen/Eigen>
#include <Eigen/QR>

class Camera {
public:
    Camera();

    Eigen::Matrix<double, 3, 4> get_camera_matrix();

    void set_camera_matrix(Eigen::Matrix<double, 3, 4> new_camera_matrix);

    void camera_parameters(Eigen::Matrix<double, 3, 4> &T);
    Eigen::Matrix<double, 3, 4> camera_DLT(Eigen::Matrix<double, 6, 4> &original_points, Eigen::Matrix<double, 6, 3> &projected_points);
    Eigen::Matrix<double, 12, 12> matrix_for_DLT(Eigen::Matrix<double, 6, 4> &original_points, Eigen::Matrix<double, 6, 3> &projected_points);

    void print_matrix();
    void normalize_matrix();
private:
    Eigen::Matrix<double, 3, 4> m_camera_matrix;
    Eigen::Matrix3d m_K, m_A;
    Eigen::Vector3d m_C;
};

void normalize_matrix_snd(Eigen::Matrix3d A);

#endif //CAMERA_CAMERA_H
