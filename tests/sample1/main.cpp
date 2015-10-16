#include <lp-engine/lpengine.h>
#include <iostream>

int main() {
  lpproblem problem;
//  Eigen::MatrixXd A(3,3);
//  A << 2.0f, 3.0f, 1.0f,
//       4.0f, 1.0f, 2.0f,
//       3.0f, 4.0f, 2.0f;

//  Eigen::VectorXd b(3);
//  b << 5.0f, 11.0f, 8.0f;

//  Eigen::VectorXd c(3);
//  c << 5.0f, 4.0f, 3.0f;

  Eigen::MatrixXd A(3,2);
  A << 1.0f, -1.0f,
       2.0f, -1.0f,
       0.0f, 1.0f;

  Eigen::VectorXd b(3);
  b << 1.0f, 3.0f, 5.0f;

  Eigen::VectorXd c(2);
  c << 4.0f, 3.0f;

  problem.Set_A(A);
  problem.Set_b(b);
  problem.Set_c(c);

  problem.UpdateProblem();

//  Eigen::MatrixXd AA(5,5  );
//  AA.setRandom(AA.rows(),AA.cols());
//  std::cout << "inv:" << AA.inverse() << std::endl;

  lpengine my_engine(problem);
  my_engine.Solve();

  return 0;
}
