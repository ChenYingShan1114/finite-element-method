#include <dV_linear_tetrahedron_dq.h>

#include <dphi_linear_tetrahedron_dX.h>
#include <dpsi_neo_hookean_dF.h>
#include <quadrature_single_point.h>
#include <iostream>
//Input:
//  q - generalized coordinates for the FEM system
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  element - the 1x4 vertex indices for this tetrahedron
//  v0 - the undeformed tetrahedron volume
//  C,D - material parameters for the Neo-Hookean model
//Output:
//  dV - the 12x1 gradient of the potential energy for a single tetrahedron
void dV_linear_tetrahedron_dq(Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D) {

   auto neohookean_linear_tet = [&](Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
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
        // for (int i = 0; i < F.rows(); i++) {
        //     for (int j = 0; j < F.cols(); j++) {
        //         if (F(i, j) < 1e-14) F(i, j) = 0;
        //     }
        // }
        // std::cout << F << std::endl;

        // calculate gradient of Neohookean strain energy
        // use F into the equation and put the vector dV
        Eigen::Vector9d dw_9 = Eigen::Vector9d::Zero(); // 9x1
        dpsi_neo_hookean_dF(dw_9, F, C, D);

        Eigen::MatrixXd B;
        B.resize(dw_9.size(), dV.size());
        B.setZero();
        
        for (int i = 0; i < dV.rows() / dphi.rows(); i++) {
            for (int j = 0; j < dV.size(); j++) {
                B(j % 3 + 3 * i, 3 * (j / 3) + i) = dphi(j / 3, j % 3);
            }
        }
        dV = B.transpose() * dw_9;

        
        // std::cout << B << std::endl;
        // std::cout << dw_9 << std::endl;
        
        // std::cout << "neohookean linear dV: " << dV << std::endl;

    };

    quadrature_single_point(dV, q, element, volume, neohookean_linear_tet);  
    // std::cout << volume << " dV_linear_tetrahedron_dq dV: " << dV << std::endl;

    
}