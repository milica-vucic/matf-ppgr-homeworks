#include "camera.h"

#include <utility>

Camera::Camera()
{}

Eigen::Matrix<double, 3, 4> Camera::get_camera_matrix()
{
    return m_camera_matrix;
}

void Camera::set_camera_matrix(Eigen::Matrix<double, 3, 4> new_camera_matrix)
{
    m_camera_matrix = std::move(new_camera_matrix);
}

void Camera::camera_parameters(Eigen::Matrix<double, 3, 4> &T)
{
    Eigen::Matrix3d T0 = T.block<3, 3>(0, 0);
    if (T0.determinant() < 0)
        T = -1 * T;

    Eigen::Matrix3d tmp;
    double c4;

    tmp << T.col(1), T.col(2), T.col(3);
    m_C(0) = (tmp.transpose()).determinant();

    tmp << T.col(0), T.col(2), T.col(3);
    m_C(1) = -(tmp.transpose()).determinant();

    tmp << T.col(0), T.col(1), T.col(3);
    m_C(2) = (tmp.transpose()).determinant();

    c4 = -T0.determinant();

    m_C(0) = m_C(0) / c4;
    m_C(1) = m_C(1) / c4;
    m_C(2) = m_C(2) / c4;

    Eigen::HouseholderQR<Eigen::Matrix3d> qr(T0.inverse());
    Eigen::Matrix3d Q = qr.householderQ() * Eigen::Matrix3d::Identity();
    Eigen::Matrix3d R = qr.matrixQR().triangularView<Eigen::Upper>();

    if (R(0, 0) < 0) {
        R.row(0) = -1 * R.row(0);
        Q.col(0) = -1 * Q.col(0);
    }

    if (R(1, 1) < 0) {
        R.row(1) = -1 * R.row(1);
        Q.col(1) = -1 * Q.col(1);
    }

    if (R(2, 2) < 0) {
        R.row(2) = -1 * R.row(2);
        Q.col(2) = -1 * Q.col(2);
    }

    m_A = Q.transpose();
    m_K = R.inverse();

    ::normalize_matrix_snd(m_K);
}

Eigen::Matrix<double, 3, 4> Camera::camera_DLT(Eigen::Matrix<double, 6, 4> &original_points, Eigen::Matrix<double, 6, 3> &projected_points)
{
    Eigen::Matrix<double, 12, 12> matrix = matrix_for_DLT(original_points, projected_points);
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(matrix, Eigen::ComputeFullV);
    Eigen::MatrixXd V = svd.matrixV();
    unsigned k = V.cols();

    Eigen::MatrixXd result_matrix(3, 4);
    for (unsigned i = 0; i < 3; ++i) {
        for (unsigned j = 0; j < 4; ++j) {
            result_matrix(i, j) = V(4 * i + j, k - 1);
        }
    }

    return result_matrix;
}

Eigen::Matrix<double, 12, 12> Camera::matrix_for_DLT(Eigen::Matrix<double, 6, 4> &original_points, Eigen::Matrix<double, 6, 3> &projected_points)
{
    Eigen::MatrixXd matrix(12, 12);
    unsigned j = 0;
    for (unsigned i = 0; i < original_points.rows(); ++i) {
        matrix.row(j) << 0.0, 0.0, 0.0, 0.0, -projected_points(i, 2) * original_points(i, 0), -projected_points(i, 2) * original_points(i, 1),
                -projected_points(i, 2) * original_points(i, 2), -projected_points(i, 2) * original_points(i, 3),
                projected_points(i, 1) * original_points(i, 0), projected_points(i, 1) * original_points(i, 1),
                projected_points(i, 1) * original_points(i, 2), projected_points(i, 1) * original_points(i, 3);

        matrix.row(j + 1) << projected_points(i, 2) * original_points(i, 0), projected_points(i, 2) * original_points(i, 1),
                projected_points(i, 2) * original_points(i, 2), projected_points(i, 2) * original_points(i, 3),
                0.0, 0.0, 0.0, 0.0, -projected_points(i, 0) * original_points(i, 0), -projected_points(i, 0) * original_points(i, 1),
                -projected_points(i, 0) * original_points(i, 2), -projected_points(i, 0) * original_points(i, 3);
        j += 2;
    }

    return matrix;
}

void Camera::normalize_matrix()
{
    auto element = m_camera_matrix(0, 0);
    for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 4; ++j)
            m_camera_matrix(i, j) /= element;
}

void Camera::print_matrix()
{
    int max_length = 0;
    size_t n = m_camera_matrix.rows();
    size_t m = m_camera_matrix.cols();

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            std::string number = std::to_string(m_camera_matrix(i,j));
            max_length = (int)number.size() > max_length ? (int)number.size() : max_length;
        }
    }

    for (int i = 0; i < n; ++i) {
        std::cout << "[";
        for (int j = 0; j < m; ++j) {
            std::string number = std::to_string(m_camera_matrix(i,j));
            if (j > 0) {
                std::cout << std::string(max_length - number.size() + 1, ' ') << number;
            } else {
                std::cout << number << std::string(max_length - number.size(), ' ');
            }
        }
        std::cout << "]\n";
    }

    std::cout << std::endl;
}

void normalize_matrix_snd(Eigen::Matrix3d A)
{
    double element = A(2, 2);
    for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j)
            A(i, j) /= element;
}