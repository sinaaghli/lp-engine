#include <lp-engine/lpengine.h>

lpengine::lpengine(lpproblem& problem) : p_(problem) {}

lpengine::~lpengine() {}

bool lpengine::CheckDicOptimal() { return false; }

bool lpengine::CheckDicBounded() { return true; }

bool lpengine::CheckDicfeasible(Eigen::VectorXd& x_b_hat) {
  for (int ii = 0; ii < x_b_hat.rows(); ii++) {
    if (!(x_b_hat[ii] >= 0)) {
      return false;
    }
  }
  return true;
}

int lpengine::ChooseEnteringVar(Eigen::VectorXd& vec) {
  for (int ii = 0; ii < vec.rows(); ii++) {
    if (vec[ii] < 0) {
      return ii;
    }
  }
  std::cout << "-- Error Didn't find entring variable. Halting ..."
            << std::endl;
  // TODO: this should never happen but figure out how to handle this situation.
  exit(EXIT_FAILURE);
}

int lpengine::ChooseLeavingVar(Eigen::VectorXd& delta_x_b,
                               Eigen::VectorXd& x_b_hat) {
  Eigen::VectorXd tmp(delta_x_b.rows());
  for (int ii = 0; ii < delta_x_b.rows(); ++ii) {
    if (x_b_hat[ii] == 0) {
      tmp[ii] = 0;
    } else {
      tmp[ii] = delta_x_b[ii] / x_b_hat[ii];
    }
  }
  std::ptrdiff_t max_index;
  tmp.maxCoeff(&max_index);

  // calculate primal step size
  p_.t_ = 1 / tmp[max_index];

  return max_index;
}

void lpengine::PrintStatusMsg(LPPRoblemStatus& status) {
  std::cout << "-- Problem Status : ";
  switch (status) {
    case STATUS_INFEASIBLE:
      std::cout << "INFEASIBLE" << std::endl;
      break;
    case STATUS_OPTIMAL:
      std::cout << "OPTIMAL" << std::endl;
      break;
    case STATUS_SOLVING:
      std::cout << "SOLVING" << std::endl;
      break;
    case STATUS_UNBOUNDED:
      std::cout << "UNBOUNDED" << std::endl;
      break;
  }
}

bool lpengine::IsSemiPositive(Eigen::VectorXd& vec) {
  for (int ii = 0; ii < vec.rows(); ii++) {
    if (vec[ii] < 0) {
      return false;
    }
  }
  return true;
}

void lpengine::SwapCols(Eigen::MatrixXd& mat1, int index_mat1,
                        Eigen::MatrixXd& mat2, int index_mat2) {
  Eigen::VectorXd tmp(mat1.rows());
  tmp = mat1.col(index_mat1);
  mat1.col(index_mat1) = mat2.col(index_mat2);
  mat2.col(index_mat2) = tmp;
}

LPPRoblemStatus lpengine::SimplexStep() {
  // Step1: Check if Dual is feasible then dictionary is optimal
  p_.dic_number_++;

  std::cout << "a_b_\n" << p_.a_b_ << std::endl;
  std::cout << "a_n_\n" << p_.a_n_ << std::endl;
  std::cout << "x_b_hat\n" << p_.x_b_hat_ << std::endl;
  std::cout << "z_n_hat\n" << p_.z_n_hat_ << std::endl;

  Eigen::MatrixXd tmp(p_.dim_m_, p_.dim_n_);
  tmp = p_.a_b_.inverse() * p_.a_n_;
  p_.z_n_hat_ = tmp.transpose() * p_.c_b_ - p_.c_n_;
  if (IsSemiPositive(p_.z_n_hat_)) {
    // current dictionary is optimal
    return LPPRoblemStatus::STATUS_OPTIMAL;
  }

  // Step2: Select Entring Variable
  p_.entering_var_ = ChooseEnteringVar(p_.z_n_hat_);

  std::cout << "entering\n" << p_.entering_var_ << std::endl;

  // Step3: compute Primal step direction delta(x_b_)
  Eigen::VectorXd ej(p_.dim_m_);
  ej.setZero(ej.rows());
  ej[p_.entering_var_] = 1;
  p_.delta_x_b_ = tmp * ej;

  std::cout << "delta x_b\n" << p_.delta_x_b_ << std::endl;

  // Step4&5: Compute primal step length and leaving var
  p_.leaving_var_ = ChooseLeavingVar(p_.delta_x_b_, p_.x_b_hat_);

  std::cout << "leaving var\n" << p_.leaving_var_ << std::endl;
  std::cout << "t is\n" << p_.t_ << std::endl;

  // Step6: choose dual Step Direction. dual_enteringvar=dual_leavingvar
  Eigen::VectorXd ei(p_.dim_n_);
  ei.setZero(ej.rows());
  ei[p_.leaving_var_] = 1;
  p_.delta_z_n_ = -tmp.transpose() * ei;

  std::cout << "temp trans" << tmp.transpose() << std::endl;
  std::cout << "ei\n" << ei << std::endl;
  std::cout << "delta z_n\n" << p_.delta_z_n_ << std::endl;

  // Step7: Compute dual step length dual_leavingvar=primal_enteringvar
  p_.s_ = p_.z_n_hat_[p_.entering_var_] / p_.delta_z_n_[p_.entering_var_];

  std::cout << "s\n" << p_.s_ << std::endl;

  // Step8: Update current primal and dual solutions
  p_.x_b_hat_ = p_.x_b_hat_ - p_.t_ * p_.delta_x_b_;
  p_.z_n_hat_ = p_.z_n_hat_ - p_.s_ * p_.delta_z_n_;

  std::cout << "x_b_hat\n" << p_.x_b_hat_ << std::endl;
  std::cout << "z_n_hat\n" << p_.z_n_hat_ << std::endl;

  // Step9: Update Basis
  SwapCols(p_.a_n_,p_.entering_var_,p_.a_b_,p_.leaving_var_);

  std::cout << "a_b_\n" << p_.a_b_ << std::endl;
  std::cout << "a_n_\n" << p_.a_n_ << std::endl;

  p_.x_b_hat_[p_.leaving_var_] = p_.t_;
  p_.z_n_hat_[p_.entering_var_] = p_.s_;

  std::cout << "x_b_hat\n" << p_.x_b_hat_ << std::endl;
  std::cout << "z_n_hat\n" << p_.z_n_hat_ << std::endl;

  exit(EXIT_FAILURE);
  return LPPRoblemStatus::STATUS_SOLVING;
}

LPPRoblemStatus lpengine::Solve() {
  // Do a initial feasibility check and boundedness check
  p_.x_b_hat_ = p_.a_b_.inverse() * p_.b_;
  if (!CheckDicfeasible(p_.x_b_hat_)) {
    std::cout << "-- Initial Dictionary is not feasible. Exiting Solver ..."
              << std::endl;
    return LPPRoblemStatus::STATUS_INFEASIBLE;
  }

  if (!CheckDicBounded(/*TODO:implement*/)) {
    std::cout << "-- Initial Dictionary is not bounded. Exiting Solver ..."
              << std::endl;
    return LPPRoblemStatus::STATUS_UNBOUNDED;
  }

  // Problem is feasible and bounded at this point then do Simplex steps
  LPPRoblemStatus current_status = LPPRoblemStatus::STATUS_SOLVING;
  while (current_status == LPPRoblemStatus::STATUS_SOLVING) {
    current_status = SimplexStep();
  }

  // Pivot
  ////////////////////////////

  return LPPRoblemStatus::STATUS_OPTIMAL;
}
