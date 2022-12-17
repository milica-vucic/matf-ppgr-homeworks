#include <iostream>

#include "camera.h"

#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <Eigen/Eigen>
#include <Eigen/QR>

struct Point {
    double x, y;
    double z = 1.0;
};

void on_click(int event, int x, int y, int flags, void* user_data);
Eigen::Matrix<double, 6, 3> vector_to_matrix(std::vector<Point> points);

std::vector<Point> projected_points_3D;
cv::Mat image, image_resized;

int main ()
{
    Eigen::Matrix<double, 6, 4> original_points;
    Eigen::Matrix<double, 6, 3> projected_points;

    original_points << 25, 87, 208, 1,
                       45, 229, 0, 1,
                       188, 163, 0, 1,
                       113, 227, 139, 1,
                       82, 0, 14, 1,
                       6, 0, 279, 1;

    image = cv::imread("/home/milica/Desktop/PPGR/matf-ppgr-homeworks/homework_04/images/scene_pixels_01.jpg", cv::IMREAD_COLOR);
    if (image.empty()) {
        std::cout << "Can't find/load image!" << std::endl;
        return -1;
    }

    cv::resize(image, image_resized, cv::Size(700, 900), cv::INTER_AREA);

    cv::imshow("3D scene", image_resized);
    cv::setMouseCallback("3D scene", on_click);

    while (true) {
        if (cv::waitKey(1) && projected_points_3D.size() == 6)
            break;
    }
    cv::waitKey(0);

    Camera camera = Camera();
    projected_points = vector_to_matrix(projected_points_3D);
    Eigen::Matrix<double, 3, 4> DLT_matrix = camera.camera_DLT(original_points, projected_points);
    camera.set_camera_matrix(DLT_matrix);

    std::cout << "Camera matrix for 3D scene: " << std::endl;
    camera.normalize_matrix();
    camera.print_matrix();

    return 0;
}

void on_click(int event, int x, int y, int flags, void* user_data)
{
    if (event == cv::EVENT_LBUTTONDOWN) {
        Point point;
        point.x = x * 1.0;
        point.y = y * 1.0;
        projected_points_3D.push_back(point);

        cv::putText(image_resized, "(" + std::to_string(x) + ", " + std::to_string(y) + ")",
                    cv::Point(x, y), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::COLOR_BGR2BGR555, 2);
        cv::imshow("3D scene", image_resized);
    }
}

Eigen::Matrix<double, 6, 3> vector_to_matrix(std::vector<Point> points)
{
    Eigen::MatrixXd matrix(6, 3);
    unsigned i = 0;
    for (Point p : points)
        matrix.row(i++) << p.x, p.y, p.z;

    return matrix;
}