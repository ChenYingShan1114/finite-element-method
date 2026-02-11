#include <d2V_linear_tetrahedron_dq2.h>
#include <dphi_linear_tetrahedron_dX.h>
#include <d2psi_neo_hookean_dq2.h>
#include <quadrature_single_point.h>
#include <iostream>
//Input:
//  q - generalized coordinates for the FEM system
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  element - the 1x4 vertex indices for this tetrahedron
//  v0 - the undeformed tetrahedron volume
//  C,D - material parameters for the Neo-Hookean model
//Output:
//  dV - the 12x12 Hessian of the potential energy for a single tetrahedron
void d2V_linear_tetrahedron_dq2(Eigen::Matrix1212d &H, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D) {

    auto neohookean_linear_tet = [&](Eigen::Matrix1212d &dV, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
        
        //Code to compute non-integrated hessian matrix goes here
        // Selection matrix S
        // Eigen::SparseMatrixd S;
        // S.resize(V.rows(), element.size());
        // S.setZero();
        // for (int i = 0; i < element.size(); i++) {
        //     S.coeffRef(element(i), i) = 1;
        // }

        // calculate dphi for deformation gradient F
        Eigen::Matrix43d dphi = Eigen::Matrix43d::Zero();
        dphi_linear_tetrahedron_dX(dphi, V, element, X);
        // std::cout << dphi << std::endl << std::endl;

        // deformation gradient F
        // Eigen::Matrix3d F = Eigen::Matrix3d::Zero();
        // F = V.transpose() * S * dphi;
        // std::cout << F << std::endl;

        Eigen::Matrix3d F = Eigen::Matrix3d::Zero();
        Eigen::MatrixXd test;
        test.resize(3, 4);
        test.setZero();
        for (int i = 0; i < test.rows(); i++) {
            for (int j = 0; j < test.cols(); j++) {
                test(i, j) = q(3 * element(j) + i);
            }
        }
        // test << q(3 * element(0)    ), q(3 * element(1)    ), q(3 * element(2)    ), q(3 * element(3)    ),
        //         q(3 * element(0) + 1), q(3 * element(1) + 1), q(3 * element(2) + 1), q(3 * element(3) + 1),
        //         q(3 * element(0) + 2), q(3 * element(1) + 2), q(3 * element(2) + 2), q(3 * element(3) + 2);
        // std::cout << dphi << std::endl;
        // std::cout << V.row(element(0)) << " " << V.row(element(1)) << " " << V.row(element(2)) << " " << V.row(element(3)) << std::endl;
        // std::cout << test << std::endl;
        F = test * dphi;
        // std::cout << F << std::endl;

        // calculate second Piola-Kirchhoff
        // use F into the equation and put the vector dV
        Eigen::Matrix99d ddw_9 = Eigen::Matrix99d::Zero();
        d2psi_neo_hookean_dF2(ddw_9, F, C, D);
        // std::cout << ddw_9 << std::endl;

        // ddw_9 是針對deformation gradiant算出來的second piola-kirchhoff
        // 接下來要將dww 乘回對應的四面體上，利用dphi轉換

        // dV    = B.transpose() * ddw_9 * B ?
        // 12x12 =      12x9     * 9x9 * 9x12
        Eigen::MatrixXd B;
        B.resize(3 * dphi.cols(), element.size() * dphi.cols());
        B.setZero();
        // B.resize(ddw_9.rows(), dV.cols());
        // std::cout << B.rows() << " " << B.cols() << std::endl;

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 12; j++) {
                B(j % 3 + 3 * i, 3 * (j / 3) + i) = dphi(j / 3, j % 3);
            }
        }
        dV = B.transpose() * ddw_9 * B;
        // std::cout << B << std::endl;
        // std::cout << "neohookean_linear_tet dV: " << std::endl;
        // std::cout << dV << std::endl;

    };

    //integrate the non-integrated hessian across the tetrahedral element
    quadrature_single_point(H, q, element, volume, neohookean_linear_tet);  
    // std::cout << "d2V_linear_tetrahedron_dq2 dV: " << H << std::endl;

    //DO NOT REMOVE THIS CODE This code ensures that the hessian matrix is symmetric postive definite by projecting all
    //negative eigenvalues to small, postive values.
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix1212d> es(H);
    
    Eigen::MatrixXd DiagEval = es.eigenvalues().real().asDiagonal();
    Eigen::MatrixXd Evec = es.eigenvectors().real();
    
    for (int i = 0; i < 12; ++i) {
        if (es.eigenvalues()[i]<1e-6) {
            DiagEval(i,i) = 1e-3;
        }
    }

    H = Evec * DiagEval * Evec.transpose();
    // std::cout << "d2V_linear_tetrahedron_dq2 dV: " << std::endl;
    // std::cout << H << std::endl << std::endl;

}
