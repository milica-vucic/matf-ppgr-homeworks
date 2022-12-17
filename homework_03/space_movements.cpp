#include <iostream>
#include <cmath>
#include <Eigen/Eigen>

void Euler2A(const double &phi, const double &theta, const double &psi, Eigen::Matrix3d &A);
void AxisAngle(const Eigen::Matrix3d &A, Eigen::Vector3d &p, double &phi);
void Rodriguez(Eigen::Vector3d &p, const double &phi, Eigen::Matrix3d &A);
void A2Euler(const Eigen::Matrix3d &A, double &phi, double &theta, double &psi);
void AxisAngle2Q(Eigen::Vector3d p, const double &phi, Eigen::Vector4d &q);
void Q2AxisAngle(Eigen::Vector4d &q, double &phi, Eigen::Vector3d &p);
void print_matrix(const Eigen::MatrixXd &matrix);

int main ()
{
    // double phi = 0.0; // -atan(1.0 / 4.0);
    // double theta = M_PI / 2.0; // -asin(8.0 / 9.0);
    // double psi = 0.0; //atan(4.0);

    double phi = ((2 * M_PI) / 6) * (363 % 5 + 1);
    double theta = (M_PI / 17) * (363 % 8 + 1);
    double psi = ((2 * M_PI) / 8) * (363 % 7 + 1);

    double phi_1, theta_1, psi_1;

    Eigen::Vector3d p;
    Eigen::Vector4d q;
    Eigen::Matrix3d A, R;
    double angle;

    std::cout << "--- Euler2A ---" << std::endl;
    Euler2A(phi, theta, psi, A);
    print_matrix(A);

    std::cout << "--- AxisAngle ---" << std::endl;
    AxisAngle(A, p, angle);
    std::cout << "Vector p: ";
    for (unsigned i = 0; i < 3; ++i)
        std::cout << p(i) << " ";
    std::cout << std::endl << "Φ: " << angle << std::endl << std::endl;

    std::cout << "--- Rodriguez ---" << std::endl;
    Rodriguez(p, angle, R);
    print_matrix(R);

    std::cout << "--- A2Euler ---" << std::endl;
    A2Euler(A, phi_1, theta_1, psi_1);
    std::cout << "Φ: " << phi_1 << std::endl
              << "Θ: " << theta_1 << std::endl
              << "Ψ: " << psi_1 << std::endl << std::endl;

    std::cout << "--- AxisAngle2Q ---" << std::endl;
    AxisAngle2Q(p, angle, q);
    std::cout << "q = " << q(0) << "i + " << q(1) << "j + " << q(2) << "k + " << q(3) << std::endl << std::endl;

    std::cout << "--- Q2AxisAngle ---" << std::endl;
    Q2AxisAngle(q, angle, p);
    std::cout << "Vector p: ";
    for (unsigned i = 0; i < 3; ++i)
        std::cout << p(i) << " ";
    std::cout << std::endl << "Φ: " << angle << std::endl << std::endl;

    return 0;
}

void Euler2A(const double &phi, const double &theta, const double &psi, Eigen::Matrix3d &A)
{
    Eigen::Matrix3d R_x, R_y, R_z;

    R_x << 1.0, 0.0,       0.0,
           0.0, cos(phi), -sin(phi),
           0.0, sin(phi),  cos(phi);

    R_y << cos(theta),  0.0, sin(theta),
           0.0,         1.0, 0.0,
           -sin(theta), 0.0, cos(theta);

    R_z << cos(psi), -sin(psi), 0.0,
           sin(psi),  cos(psi), 0.0,
           0.0     ,  0.0     , 1.0;

    A = R_z * R_y * R_x;
}

void AxisAngle(const Eigen::Matrix3d &A, Eigen::Vector3d &p, double &phi)
{
    bool isDeterminantOne = round(A.determinant()) == 1;
    bool isIdentityMatrix = (A.transpose() * A).isIdentity();

    Eigen::Matrix3d S = A.transpose() * A;
    for (unsigned i = 0; i < 3; ++i) {
        for (unsigned j = 0; j < 3; ++j) {
            std::cout << S(i, j) << " ";
        }
        std::cout << std::endl;
    } 

    if (!isDeterminantOne || !isIdentityMatrix) {
        std::cerr << "A is not movement matrix!" << std::endl;
        return;
    }

    Eigen::Matrix3d E;
    E << 1.0, 0.0, 0.0,
         0.0, 1.0, 0.0,
         0.0, 0.0, 1.0;

    Eigen::Matrix3d X = A - E;
    Eigen::Vector3d zeroVector = {0.0, 0.0, 0.0};
    Eigen::Vector3d u = {X(0, 0), X(0, 1), X(0, 2)};
    Eigen::Vector3d v = {X(1, 0), X(1, 1), X(1, 2)};

    u.normalize();
    v.normalize();
    p = u.cross(v);
    if (p == zeroVector) {                    
        v = {X(2, 0), X(2, 1), X(2, 2)};
        v.normalize();
        p = u.cross(v);
    }

    p.normalize();

    Eigen::Vector3d up = A * u;
    phi = acos(u.dot(up));

    Eigen::Matrix3d D;
    D << p, u, up;
    if (D.determinant() < 0)
        p = -p;
}

void Rodriguez(Eigen::Vector3d &p, const double &phi, Eigen::Matrix3d &A)
{
    if (p.norm() != 1)
        p.normalize();

    Eigen::Matrix3d E;
    E << 1.0, 0.0, 0.0,
         0.0, 1.0, 0.0,
         0.0, 0.0, 1.0;

    Eigen::Matrix3d px;
    px << 0.0, -p(2), p(1),
          p(2), 0.0, -p(0),
         -p(1), p(0), 0.0;
    
    Eigen::Matrix3d ppt = p * p.transpose();

    A = ppt + cos(phi) * (E - ppt) + sin(phi) * px;
}

void A2Euler(const Eigen::Matrix3d &A, double &phi, double &theta, double &psi)
{
    if (round(A.determinant()) != 1) {
        std::cerr << "A is not movement matrix!" << std::endl;
        return;
    }

    if (A(2, 0) < 1) {
        if (A(2, 0) > -1) {
            psi   = atan2(A(1, 0), A(0, 0));
            theta = asin(-A(2, 0));
            phi   = atan2(A(2, 1), A(2, 2));
        } else {
            psi   = atan2(-A(0, 1), A(1, 1));
            theta = M_PI / 2;
            phi   = 0;
        }
    } else {
        psi   = atan2(-A(0, 1), A(1, 1));
        theta = -M_PI / 2;
        phi   = 0;
    }
}

void AxisAngle2Q(Eigen::Vector3d p, const double &phi, Eigen::Vector4d &q)
{
    double w = cos(phi / 2);
    p.normalize();

    auto x = sin(phi / 2) * p(0);
    auto y = sin(phi / 2) * p(1);
    auto z = sin(phi / 2) * p(2);

    q << x, y, z, w;
}

void Q2AxisAngle(Eigen::Vector4d &q, double &phi, Eigen::Vector3d &p)
{
    if (q.norm() != 1)
        q.normalize();

    if (q(3) < 0)
        q = -q;

    phi = 2 * acos(q(3));
    if (abs(q(3)) == 1) {
        p << 1.0, 0.0, 0.0;
    } else {
        p << q(0), q(1), q(2);
        p.normalize();
    }
}

void print_matrix(const Eigen::MatrixXd &matrix)
{
    int max_length = 0;
    int n = matrix.rows();
    int m = matrix.cols();

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
