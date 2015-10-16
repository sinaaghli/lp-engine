#ifndef LP_ENGINE_H_
#define LP_ENGINE_H_

#include <iostream>
#include <Eigen/Eigen>
#include <lp-engine/lpproblem.h>
#include <cstddef>

class lpengine{
 public:
  lpengine(lpproblem& problem);
  ~lpengine();

  LPPRoblemStatus Solve();

 private:
  bool ReadDicSolution();
  bool CheckDicOptimal();
  bool CheckDicBounded();
  bool IsSemiPositive(Eigen::VectorXd& vec);
  bool CheckDicfeasible(Eigen::VectorXd& x_b_hat);
  int ChooseEnteringVar(Eigen::VectorXd& vec);
  int ChooseLeavingVar(Eigen::VectorXd& delta_x_b, Eigen::VectorXd& x_b_hat);
  LPPRoblemStatus SimplexStep();
  void SwapCols(Eigen::MatrixXd& mat1, int index_mat1, Eigen::MatrixXd& mat2, int index_mat2);
  void PrintStatusMsg(LPPRoblemStatus& status);
  lpproblem p_;
};

#endif  // LP_ENGINE_H_
