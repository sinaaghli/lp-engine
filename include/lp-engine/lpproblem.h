#ifndef LPPROBLEM_H_
#define LPPROBLEM_H_

#include <Eigen/Eigen>
#include <iostream>
#include <vector>
#include <lp-engine/lpdictionary.h>

enum LPSolverType { NONE, LP_SIMPLEX_PRIMAL, LP_SIMPLEX_DUAL };
enum LPPRoblemStatus {
  STATUS_SOLVING,
  STATUS_UNBOUNDED,
  STATUS_INFEASIBLE,
  STATUS_OPTIMAL
};

class lpproblem {
 public:
  lpproblem();
  ~lpproblem();
  void Set_A(Eigen::MatrixXd A_mat);
  void Set_b(Eigen::VectorXd b_vec);
  void Set_c(Eigen::VectorXd c_vec);
  bool UpdateProblem();
  bool Set_solver(LPSolverType solver_type);
  void PrintDictionary();

  // private:
  Eigen::MatrixXd A_start_;
  Eigen::VectorXd b_start_;
  Eigen::VectorXd c_start_;
  LPSolverType solver_type_;

  Eigen::MatrixXd a_b_;
  Eigen::MatrixXd a_n_;
  Eigen::VectorXd x_b_;
  Eigen::VectorXd x_n_;
  Eigen::VectorXd c_b_;
  Eigen::VectorXd c_n_;
  Eigen::VectorXd b_;
  Eigen::VectorXd z_n_;
  Eigen::VectorXd x_b_hat_;
  Eigen::VectorXd z_n_hat_;
  Eigen::VectorXd delta_x_b_;
  Eigen::VectorXd delta_z_n_;
  // Primal step size
  double t_;
  // Dual Step size
  double s_;
  // Latest chosen leaving variable
  int leaving_var_;
  // Latest chosen entering variable
  int entering_var_;
  std::vector<int> basic_indices_;
  std::vector<int> nonbasic_indices_;

  int dim_n_;
  int dim_m_;

  int dic_number_;
};

#endif  // LPPROBLEM_H__
