#include <V_linear_tetrahedron.h>

#include <dphi_linear_tetrahedron_dX.h>
#include <psi_neo_hookean.h>
#include <quadrature_single_point.h>

#include <iostream>
//Input:
// q - generalized coordinates of FEM system
// V - vertex matrix for the mesh
// element - vertex indices of the element
// volume - volume of tetrahedron
// C,D - material parameters
//Output:
//  energy - potential energy of tetrahedron
void V_linear_tetrahedron(double &energy, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D) {

    auto neohookean_linear_tet = [&](double &e, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
        // std::cout << "neohookean linear tet X: " << X << std::endl;
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
        // V: rows 500, cols 3; S: rows 500, cols 4; V.transpose() * S: rows 3, cols 4
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

        // calculate Neohookean strain energy
        // use F into the equation and put the energy value into e
        psi_neo_hookean(e, F, C, D);
        // std::cout << "auto neohooken linear tet energy: " << e << std::endl;
    };

    quadrature_single_point(energy, q, element, volume, neohookean_linear_tet);  
    // std::cout << "V linear tetrahedron energy: " << energy << std::endl;
    
}