#include <d2psi_neo_hookean_dq2.h>
#include <iostream>
//Input:
//  F - the dense 3x3 deformation gradient
//  C,D - material parameters for the Neo-Hookean model
//Output:
//  psi - the 9x9 Hessian of the potential energy wrt to the deformation gradient
void d2psi_neo_hookean_dF2(Eigen::Matrix99d &ddw, Eigen::Ref<const Eigen::Matrix3d> F, double C, double D) {
    ddw.setZero();

    // principal stretches
    double lambda1 = F(0, 0), lambda2 = F(1, 1), lambda3 = F(2, 2);
  
    // left Cauchy-Green deformation tensor
    // Eigen::Matrix3d B;
    // B = F * F.transpose();

    // principal invariants
    // double I1 = pow(lambda1, 2) + pow(lambda2, 2) + pow(lambda3, 2);
    double I1 = F.squaredNorm();
    double I2 = pow(lambda1 * lambda2, 2) + pow(lambda1 * lambda3, 2) + pow(lambda2 * lambda3, 2);
    double I3 = pow(lambda1 * lambda2 * lambda3, 2);
    
    // Jacobian of F
    // double J = lambda1 * lambda2 * lambda3;
    double J = F.determinant();

    // F^{-T}
    Eigen::Matrix3d F_inv = F.inverse();
    Eigen::Matrix3d F_inv_T = F_inv.transpose();
    
    // std::cout << F << std::endl;
    // std::cout << lambda1 << " " << lambda2 << " " << lambda3 << " " << I1 << " " << I2 << " " << I3 << " " << J << std::endl;
    // std::cout << F_inv << std::endl;
    // std::cout << F_inv_T << std::endl;
    // std::cout << std::endl;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {

            // dF^{-T}_{ij}
            // std::cout << "ex. " << i << " " << j << std::endl;
            Eigen::Matrix3d dF_inv_T = Eigen::Matrix3d::Zero();
            dF_inv_T = -F_inv_T.col(j) * F_inv_T.row(i);
            // std::cout << dF_inv_T << std::endl;

            // delta
            Eigen::Matrix3d delta = Eigen::Matrix3d::Zero();
            delta(i, j) = 1;
            // std::cout << delta << std::endl;

            Eigen::Matrix3d dFF = Eigen::Matrix3d::Zero();
            dFF = 2 * C * pow(J, -2.0 / 3.0) * (-2.0 / 3.0 * F_inv_T * (F(i, j) - 1.0 / 3.0 * I1 * F_inv_T(i, j)) 
                            + (delta - 2.0 / 3.0 * F * F_inv_T(i, j) - 1.0 / 3.0 * I1 * dF_inv_T))
                + 2 * D * J * ((2 * J - 1) * F_inv_T * F_inv_T(i, j) + (J - 1) * dF_inv_T);
            // check when F = I, what is the result of dFF
            // std::cout << dFF << std::endl << std::endl;

            // put in ddw
            for (int k1 = 0; k1 < dFF.rows(); k1++) {
                for (int k2 = 0; k2 < dFF.cols(); k2++) {
                    // if (dFF(k1, k2) > 1e-9) {
                        ddw(3 * i + j, 3 * k1 + k2) = dFF(k1, k2);
                        // std::cout << k1 << " " << k2 << std::endl;
                        // std::cout << dFF << std::endl;
                    // }
                }
            }
            // std::cout << ddw << std::endl << std::endl;
            // std::cout << 2 * D + 8.0 / 3.0 * C << " " << 2 * D - 4.0 / 3.0 * C << " " << 2 * C << std::endl;
        }
    }
    // std::cout << ddw << std::endl << std::endl;
    // for (int i = 0; i < 9; i++) {
    //     for (int j = i + 1; j < 9; j++) {
    //         if (ddw(i, j) != ddw(j, i)) {
    //             std::cout << i << " " << j << " " << ddw(i, j) << " " << ddw(j, i) << " " << ddw(i, j) - ddw(j, i) << std::endl;
    //         }
    //     }
    // }

    // for (int i = 0; i < 3; i++) {
    //     for (int j = 0; j < 3; j++) {
    //         // ex. ij
    //         std::cout << "ex. " << i << " " << j << std::endl;
    //         int i1, i2, j1, j2;
    //         i1 = (i + 1)%3;
    //         i2 = (i + 2)%3;
    //         j1 = (j + 1)%3;
    //         j2 = (j + 2)%3;
    //         // std::cout << i1 << " " << i2 << " " << j1 << " " << j2 << std::endl;

    //         // Derivatives of cofactors that contain F(i, j)
    //         Eigen::Matrix3d tmp = Eigen::Matrix3d::Zero();
    //         tmp(i1, j1) = F(i2, j2);   // dC(i1,j1)/dF(i,j)
    //         tmp(i2, j2) = F(i1, j1);  // dC(i2,j2)/dF(i,j)
    //         tmp(i1, j2) = -F(i2, j1);  // dC(i1,j2)/dF(i,j)
    //         tmp(i2, j1) = -F(i1, j2);  // dC(i2,j1)/dF(i,j)
    //         // std::cout << tmp << std::endl << std::endl;

    //         Eigen::Matrix3d manual = Eigen::Matrix3d::Zero();
    //         manual = -F_inv_T * F_inv_T(i, j) + 1/J * tmp;
    //         // std::cout << F_inv_T << std::endl << std::endl;
    //         // std::cout << manual << std::endl << std::endl;

    //         // std::cout << F_inv_T.col(0) << " " << F_inv_T.row(2) << std::endl << std::endl;
    //         Eigen::Matrix3d identity = Eigen::Matrix3d::Zero();
    //         identity = -F_inv_T.col(j) * F_inv_T.row(i);
    //         // std::cout << -F_inv_T.col(j) * F_inv_T.row(i) << std::endl << std::endl;
    //         // std::cout << identity - manual << std::endl << std::endl;

    //         // std::cout << "check" << std::endl << std::endl;

    //     }
    // }

    // compressible neo-Hookean strain energy density -> second Piola-Kirchhoff
    // Eigen::Matrix3d dwdF00 = Eigen::Matrix3d::Zero();
    // Eigen::Matrix3d tmp00 = Eigen::Matrix3d::Zero();
    // tmp00(0, 0) = 1;
    // std::cout << F_inv_T << std::endl;
    // std::cout << F_inv_T.col(0) << " " << F_inv_T.row(0) << std::endl;
    // std::cout << F_inv_T.col(0) * F_inv_T.row(0) << std::endl;


    // std::cout << F_inv_T.row(0) * F_inv_T.col(0) << std::endl;
    // dwdF00 = 2 * C * (-2/3) * pow(J, -5/3) * (F(0, 0) - 1/3 * I1 * F_inv_T(0, 0)) * J * F_inv_T 
    //         + 2 * C * pow(J, -2/3) * (tmp00 - 1/3 * I1 * (-F_inv_T.col(0) * F_inv_T.row(0)) - 1/3 * 2 * F * F_inv_T(0, 0))
    //         + 2 * D * J * F_inv_T * (J - 1) * F_inv_T(0, 0)
    //         + 2 * D * J * J * F_inv_T * F_inv_T(0, 0)
    //         + 2 * D * J * (J - 1) * (-F_inv_T.col(0) * F_inv_T.row(0));
    
    // std::cout << dwdF00 << std::endl;





    // TEST Inverse Transpose Derivative
    // ex. 00

    // Eigen::Matrix3d F_test = Eigen::Matrix3d::Zero();
    // F_test << 1, 2, 4, 3, 5, 6, 2, 3, 6;
    // Eigen::Matrix3d M = F_test.inverse().transpose();
    // double J_test = F_test(0, 0) * F_test(1, 1) * F_test(2, 2) + F_test(0, 1) * F_test(1, 2) * F_test(2, 0) + F_test(0, 2) * F_test(1, 0) * F_test(2, 1)
    //                  - F_test(0, 0) * F_test(1, 2) * F_test(2, 1) - F_test(0, 1) * F_test(1, 0) * F_test(2, 2) - F_test(0, 2) * F_test(1, 1) * F_test(2, 0);
    // std::cout << J_test << std::endl << std::endl;

    // std::cout << F_test << std::endl << std::endl;
    // Eigen::Matrix3d tmp = Eigen::Matrix3d::Zero();
    // tmp << 0, 0, 0, 0, F_test(2, 2), -F_test(2, 1), 0, -F_test(1, 2), F_test(1, 1);
    // std::cout << tmp << std::endl << std::endl;
    
    // Eigen::Matrix3d test = Eigen::Matrix3d::Zero();
    // test = -M * M(0, 0) + 1/J_test * tmp;
    // std::cout << M << std::endl << std::endl;
    // std::cout << test << std::endl << std::endl;

    // std::cout << -M.col(0) * M.row(0) << std::endl << std::endl;

    // ex. 12
    // Eigen::Matrix3d F_test = Eigen::Matrix3d::Zero();
    // F_test << 1, 2, 4, 3, 5, 6, 2, 3, 6;
    // Eigen::Matrix3d M = F_test.inverse().transpose();
    // double J_test = F_test(0, 0) * F_test(1, 1) * F_test(2, 2) + F_test(0, 1) * F_test(1, 2) * F_test(2, 0) + F_test(0, 2) * F_test(1, 0) * F_test(2, 1)
    //                  - F_test(0, 0) * F_test(1, 2) * F_test(2, 1) - F_test(0, 1) * F_test(1, 0) * F_test(2, 2) - F_test(0, 2) * F_test(1, 1) * F_test(2, 0);
    // std::cout << J_test << std::endl << std::endl;

    // std::cout << F_test << std::endl << std::endl;
    // Eigen::Matrix3d tmp = Eigen::Matrix3d::Zero();
    // tmp << -F_test(2, 1), F_test(2, 0), 0, 0, 0, 0, F_test(0, 1), -F_test(0, 0), 0;
    // std::cout << tmp << std::endl << std::endl;
    
    // Eigen::Matrix3d test = Eigen::Matrix3d::Zero();
    // test = -M * M(1, 2) + 1/J_test * tmp;
    // std::cout << M << std::endl << std::endl;
    // std::cout << test << std::endl << std::endl;

    // std::cout << -M.col(2) * M.row(1) << std::endl << std::endl;

    // ex. 20
    // Eigen::Matrix3d F_test = Eigen::Matrix3d::Zero();
    // F_test << 1, 2, 4, 
    //           3, 5, 6, 
    //           9, 7, 8;
    // Eigen::Matrix3d M = F_test.inverse().transpose();
    // double J_test = F_test(0, 0) * F_test(1, 1) * F_test(2, 2) + F_test(0, 1) * F_test(1, 2) * F_test(2, 0) + F_test(0, 2) * F_test(1, 0) * F_test(2, 1)
    //                  - F_test(0, 0) * F_test(1, 2) * F_test(2, 1) - F_test(0, 1) * F_test(1, 0) * F_test(2, 2) - F_test(0, 2) * F_test(1, 1) * F_test(2, 0);
    // std::cout << J_test << std::endl << std::endl;

    // std::cout << F_test << std::endl << std::endl;
    // Eigen::Matrix3d tmp = Eigen::Matrix3d::Zero();
    // tmp << 0, F_test(1, 2), -F_test(1, 1), 0, -F_test(0, 2), F_test(0, 1), 0, 0, 0;
    // std::cout << tmp << std::endl << std::endl;
    
    // Eigen::Matrix3d test = Eigen::Matrix3d::Zero();
    // test = -M * M(2, 0) + 1/J_test * tmp;
    // std::cout << M << std::endl << std::endl;
    // std::cout << test << std::endl << std::endl;

    // std::cout << M.col(0) << " " << M.row(2) << std::endl << std::endl;

    // std::cout << -M.col(0) * M.row(2) << std::endl << std::endl;

    // Eigen::Matrix3d F_test = Eigen::Matrix3d::Zero();
    // F_test << 1, 2, 4, 3, 5, 6, 9, 7, 8;
    // Eigen::Matrix3d M = F_test.inverse().transpose();
    // double J_test = F_test(0, 0) * F_test(1, 1) * F_test(2, 2) + F_test(0, 1) * F_test(1, 2) * F_test(2, 0) + F_test(0, 2) * F_test(1, 0) * F_test(2, 1)
    //                  - F_test(0, 0) * F_test(1, 2) * F_test(2, 1) - F_test(0, 1) * F_test(1, 0) * F_test(2, 2) - F_test(0, 2) * F_test(1, 1) * F_test(2, 0);
    // std::cout << J_test << std::endl << std::endl;
    // std::cout << F_test << std::endl << std::endl;


    // for (int i = 0; i < 3; i++) {
    //     for (int j = 0; j < 3; j++) {
    //         // ex. ij
    //         std::cout << "ex. " << i << " " << j << std::endl;
    //         int i1, i2, j1, j2;
    //         i1 = (i + 1)%3;
    //         i2 = (i + 2)%3;
    //         j1 = (j + 1)%3;
    //         j2 = (j + 2)%3;
    //         // std::cout << i1 << " " << i2 << " " << j1 << " " << j2 << std::endl;

    //         // Derivatives of cofactors that contain F(i, j)
    //         // Eigen::Matrix3d tmp = Eigen::Matrix3d::Zero();
    //         // tmp(i1, j1) = F_test(i2, j2);   // dC(i1,j1)/dF(i,j)
    //         // tmp(i2, j2) = F_test(i1, j1);  // dC(i2,j2)/dF(i,j)
    //         // tmp(i1, j2) = -F_test(i2, j1);  // dC(i1,j2)/dF(i,j)
    //         // tmp(i2, j1) = -F_test(i1, j2);  // dC(i2,j1)/dF(i,j)
    //         // std::cout << tmp << std::endl << std::endl;

    //         Eigen::Matrix3d manual = Eigen::Matrix3d::Zero();
    //         manual = -M * M(i, j) + 1/J_test * tmp;
    //         // std::cout << M << std::endl << std::endl;
    //         // std::cout << manual << std::endl << std::endl;

    //         // std::cout << M.col(0) << " " << M.row(2) << std::endl << std::endl;
    //         Eigen::Matrix3d identity = Eigen::Matrix3d::Zero();
    //         identity = -M.col(j) * M.row(i);
    //         // std::cout << -M.col(j) * M.row(i) << std::endl << std::endl;
    //         std::cout << identity - manual << std::endl << std::endl;

    //         // std::cout << "check" << std::endl << std::endl;

    //     }
    // }


}