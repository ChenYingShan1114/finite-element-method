#include <phi_linear_tetrahedron.h>
#include <iostream>
//Input:
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  element - the 1x4 vertex indices for this tetrahedron
//  X - the position in the undeformed space at which to compute the energy density
//Output:
//  phi - the 4x1 values the basis functions
void phi_linear_tetrahedron(Eigen::Vector4d &phi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {

    // for (int i = 0; i < element.size(); i++) {
    //     std::cout << V.row(element(i)) << std::endl;
    // }

    Eigen::Matrix3d T;
    // T << V(element(1), 0) - V(element(0), 0), V(element(2), 0) - V(element(0), 0), V(element(3), 0) - V(element(0), 0),
    //      V(element(1), 1) - V(element(0), 1), V(element(2), 1) - V(element(0), 1), V(element(3), 1) - V(element(0), 1),
    //      V(element(1), 2) - V(element(0), 2), V(element(2), 2) - V(element(0), 2), V(element(3), 2) - V(element(0), 2);
    for (int i = 0; i < T.rows(); i++) {
        for (int j = 0; j < T.cols(); j++) {
            T(i, j) = V(element(j+1), i) - V(element(0), i);
        }
    }
    // std::cout << T << std::endl;
    Eigen::Matrix3d T_inv = T.inverse();
    // std::cout << T_inv << std::endl;

    Eigen::Vector3d dX;
    // dX << X(0) - V(element(0), 0), X(1) - V(element(0), 1), X(2) - V(element(0), 2);
    for (int i = 0; i < dX.size(); i++) {
        dX(i) = X(i) - V(element(0), i);
    }
    // std::cout << X << std::endl;
    // std::cout << dX << std::endl << std::endl;


    Eigen::Vector3d tmp_phi;
    tmp_phi = T_inv * dX;
    phi << 1 - tmp_phi(0) - tmp_phi(1) - tmp_phi(2), tmp_phi(0), tmp_phi(1), tmp_phi(2);
    // std::cout << phi << std::endl << std::endl;
}