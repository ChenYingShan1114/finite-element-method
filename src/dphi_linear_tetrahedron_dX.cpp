#include <dphi_linear_tetrahedron_dX.h>
#include <phi_linear_tetrahedron.h>
#include <iostream>

//Input:
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  element - the 1x4 vertex indices for this tetrahedron
//  X - the position in the undeformed space at which to compute the energy density
//Output:
//  dphi - the 4x3 gradient of the the basis functions wrt to X. The i'th row stores d phi_i/dX
void dphi_linear_tetrahedron_dX(Eigen::Matrix43d &dphi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
    // std::cout << "dphi linear tetragedron dX X: " << X << std::endl;

    Eigen::Matrix3d T;
    // T << V(element(1), 0) - V(element(0), 0), V(element(2), 0) - V(element(0), 0), V(element(3), 0) - V(element(0), 0),
    //      V(element(1), 1) - V(element(0), 1), V(element(2), 1) - V(element(0), 1), V(element(3), 1) - V(element(0), 1),
    //      V(element(1), 2) - V(element(0), 2), V(element(2), 2) - V(element(0), 2), V(element(3), 2) - V(element(0), 2);
    for (int i = 0; i < T.rows(); i++) {
        for (int j = 0; j < T.cols(); j++) {
            T(i, j) = V(element(j+1), i) - V(element(0), i);
        }
    }
    Eigen::Matrix3d T_inv = T.inverse();    

    Eigen::Vector3d one;
    one << 1, 1, 1;

    Eigen::Vector3d tmp_phi0;
    tmp_phi0 = -one.transpose() * T_inv;

    dphi << tmp_phi0.transpose(), T_inv;
    // std::cout << tmp_phi0 << std::endl;
    // std::cout << T_inv << std::endl;
    // std::cout << dphi << std::endl;
    

    // method 2: use function phi_linear_tetrahedron and multiply (X-X0).inverse to get T
    // Eigen::Vector4d phi = Eigen::Vector4d::Zero();
    // phi_linear_tetrahedron(phi, V, element, X);


}