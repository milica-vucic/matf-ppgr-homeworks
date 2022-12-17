#include <iostream>

#include <Eigen/Eigen>
#include <Eigen/QR>

void normalize_matrix(Eigen::Matrix3d &A);
void print_matrix(const Eigen::MatrixXd &matrix);
void normalize_matrix_snd(Eigen::Matrix<double, 3, 4> &A);

void camera_parameters(Eigen::Matrix<double, 3, 4> &T, Eigen::Matrix3d &K, Eigen::Matrix3d &A, Eigen::Vector3d &C);
Eigen::Matrix<double, 3, 4> camera_DLT(Eigen::Matrix<double, 6, 4> &original_points, Eigen::Matrix<double, 6, 3> &projected_points);
Eigen::Matrix<double, 12, 12> matrix_for_DLT(Eigen::Matrix<double, 6, 4> &original_points, Eigen::Matrix<double, 6, 3> &projected_points);

int main ()
{
    Eigen::Matrix3d K, A;
    Eigen::Vector3d C;

    double n = 3.0;
    std::cout << "n = " << std::to_string(n) << std::endl;
    Eigen::Matrix<double, 3, 4> T;
    T << 5, -1 - 2 * n, 3, 18 - 3 * n,
         0, -1, 5, 21,
         0, -1, 0, 1;

    camera_parameters(T, K, A, C);
    std::cout << "Camera center: ";
    for (unsigned i = 0; i < 3; ++i)
        std::cout << C(i) << " ";
    std::cout << std::endl;

    std::cout << "Matrix K: " << std::endl;
    print_matrix(K);

    std::cout << "Matrix A: " << std::endl;
    print_matrix(A);

    Eigen::Matrix<double, 6, 4> original_points;
    original_points << 460, 280, 250, 1,
                       50, 380, 350, 1,
                       470, 500, 100, 1,
                       380, 630, 50 * n, 1,
                       30 * n, 290, 0, 1,
                       580, 0, 130, 1;

    Eigen::Matrix<double, 6, 3> projected_points;
    projected_points << 288, 251, 1,
                        79, 510, 1,
                        470, 440, 1,
                        520, 590, 1,
                        365, 388, 1,
                        365, 20, 1;

    Eigen::Matrix<double, 3, 4> DLT_matrix = camera_DLT(original_points, projected_points);
    std::cout << "DLT matrix: " << std::endl;
    normalize_matrix_snd(DLT_matrix);
    print_matrix(DLT_matrix);

    return 0;
}

void camera_parameters(Eigen::Matrix<double, 3, 4> &T, Eigen::Matrix3d &K, Eigen::Matrix3d &A, Eigen::Vector3d &C)
{
    Eigen::Matrix3d T0 = T.block<3, 3>(0, 0);
    if (T0.determinant() < 0)
        T = -1 * T;

    Eigen::Matrix3d tmp;
    double c4;

    tmp << T.col(1), T.col(2), T.col(3);
    C(0) = (tmp.transpose()).determinant();

    tmp << T.col(0), T.col(2), T.col(3);
    C(1) = -(tmp.transpose()).determinant();

    tmp << T.col(0), T.col(1), T.col(3);
    C(2) = (tmp.transpose()).determinant();

    c4 = -T0.determinant();

    C(0) = C(0) / c4;
    C(1) = C(1) / c4;
    C(2) = C(2) / c4;

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

    A = Q.transpose();
    K = R.inverse();

    normalize_matrix(K);
}

Eigen::Matrix<double, 12, 12> matrix_for_DLT(Eigen::Matrix<double, 6, 4> &original_points, Eigen::Matrix<double, 6, 3> &projected_points)
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

Eigen::Matrix<double, 3, 4> camera_DLT(Eigen::Matrix<double, 6, 4> &original_points, Eigen::Matrix<double, 6, 3> &projected_points)
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

void normalize_matrix(Eigen::Matrix3d &A)
{
    double n = A(2, 2);
    for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j)
            A(i, j) /= n;
}

void normalize_matrix_snd(Eigen::Matrix<double, 3, 4> &A)
{
    double n = A(0, 0);
    for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 4; ++j)
            A(i, j) /= n;
}

void print_matrix(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &matrix)
{
    int max_length = 0;
    size_t n = matrix.rows();
    size_t m = matrix.cols();

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            std::string number = std::to_string(matrix(i,j));
            max_length = (int)number.size() > max_length ? (int)number.size() : max_length;
        }
    }

    for (int i = 0; i < n; ++i) {
        std::cout << "[";
        for (int j = 0; j < m; ++j) {
            std::string number = std::to_string(matrix(i,j));
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