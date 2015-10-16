#include <lp-engine/lpproblem.h>

lpproblem::lpproblem()
    : solver_type_(LPSolverType::NONE), dic_number_(0) {}
lpproblem::~lpproblem() {}

void lpproblem::Set_A(Eigen::MatrixXd A_mat) {
  A_start_.resize(A_mat.rows(), A_mat.cols());
  A_start_ << A_mat;
}

void lpproblem::Set_b(Eigen::VectorXd b_vec) {
  b_start_.resize(b_vec.rows());
  b_start_ << b_vec;
}

void lpproblem::Set_c(Eigen::VectorXd c_vec) {
  c_start_.resize(c_vec.rows());
  c_start_ << c_vec;
}

bool lpproblem::UpdateProblem() {
  dim_n_ = A_start_.cols();
  dim_m_ = A_start_.rows();

  if ((b_start_.rows()) != dim_m_ || (c_start_.rows() != dim_n_)) {
    std::cout << "-- Problem Error: Matrix/vector dimentions doesn't match"
              << std::endl;
    return false;
  }

  a_n_.resize(dim_m_, dim_n_);
  a_n_ = A_start_;
  a_b_.resize(dim_m_, dim_m_);
  a_b_.setIdentity();
  c_b_.resize(dim_m_);
  c_b_.setZero(dim_m_);
  c_n_.resize(dim_n_);
  c_n_ = c_start_;
  b_.resize(dim_m_);
  b_ = b_start_;
  x_b_hat_.resize(dim_m_);
  x_b_hat_.setZero(x_b_hat_.rows());
  x_n_.resize(dim_m_);
  x_n_.setZero(x_n_.rows());
  x_b_.resize(dim_m_);
  x_b_ = b_;  // this is also x_hat knowing Xn is zero
  z_n_hat_.resize(dim_n_);
  z_n_hat_ = -c_n_;
  z_n_.resize(dim_n_);
  delta_x_b_.resize(dim_m_);
  for (int ii = 1; ii <= dim_n_; ii++) {
    nonbasic_indices_.push_back(ii);
  }

  for (int ii = dim_n_ + 1; ii <= dim_n_ + dim_m_; ii++) {
    basic_indices_.push_back(ii);
  }

  //  std::cout << "a-basic:\n" << a_b_ << std::endl;
  //  std::cout << "a-nonbasic:\n" << a_n_ << std::endl;
  //  std::cout << "c-basic:\n" << c_b_ << std::endl;
  //  std::cout << "c-nonbasic:\n" << c_n_ << std::endl;
  //  std::cout << "b:\n" << b_ << std::endl;

  return true;
}

void lpproblem::PrintDictionary() {
  std::cout << "\tDic# " << dic_number_ << "\t|"
            << "\tDic# " << dic_number_ << "\t|"
            << std::endl;
}

bool lpproblem::Set_solver(LPSolverType solver_type) {
  solver_type_ = solver_type;
  return true;
}
