#include <T_linear_tetrahedron.h>
#include <mass_matrix_linear_tetrahedron.h>
#include <iostream>
//Input:
// qdot - generalied velocity of FEM system
// element - vertex indices of the element
// density - material density
// volume - volume of tetrahedron
//Output:
//  T - kinetic energy of tetrahedron
void T_linear_tetrahedron(double &T, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::RowVectorXi> element, double density, double volume) {
    
    Eigen::Matrix1212d M0 = Eigen::Matrix1212d::Zero();
    mass_matrix_linear_tetrahedron(M0, qdot, element, density, volume);
    
    // Method 1: build up qdot_select (haven't test yet)
    // Eigen::Vector12d qdot_select = Eigen::Vector12d::Zero();
    // std::cout << qdot_select << std::endl;
    // for (int i = 0; i < element.size(); i++) {
    //     qdot_select(3 * i    ) = qdot(3 * element(i)    );
    //     qdot_select(3 * i + 1) = qdot(3 * element(i) + 1);
    //     qdot_select(3 * i + 2) = qdot(3 * element(i) + 2); 
    // }
    // std::cout << qdot_select << std::endl;
    // std::cout <<element << std::endl;
    // std::cout << qdot.size() << " " << M0.size() << std::endl;
    // std::cout << qdot.rows() << " " << M0.rows() << std::endl;
    // T = 0.5 * qdot_select.transpose() * M0 * qdot_select;

    // Method 2: make selection matrix 
    Eigen::SparseMatrixd S;
    S.resize(3 * element.size(), qdot.size());
    S.setZero();
    // std::cout << S.rows() << " " << S.cols() << std::endl;
    for (int i = 0; i < element.size(); i++) {
        S.coeffRef(3 * i, 3 * element(i)) = 1;
        S.coeffRef(3 * i + 1, 3 * element(i) + 1) = 1;
        S.coeffRef(3 * i + 2, 3 * element(i) + 2) = 1;
    }
    // std::cout << (S*qdot).size() << std::endl;
    T = 0.5 * qdot.transpose() * S.transpose() * M0 * S * qdot;
    
}