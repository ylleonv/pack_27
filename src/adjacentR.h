#ifndef ADJACENTR_H_
#define ADJACENTR_H_
#include "distribution.h"

class AdjacentR : virtual public Logistic, Normal, Cauchit, Student, Gumbel, Gompertz{
public:

  AdjacentR();

  virtual Eigen::VectorXd inverse_logistic(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative_logistic(const Eigen::VectorXd& eta) const ;

  virtual Eigen::VectorXd inverse_normal(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative_normal(const Eigen::VectorXd& eta) const ;

  virtual Eigen::VectorXd inverse_cauchit(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative_cauchit(const Eigen::VectorXd& eta) const ;

  virtual Eigen::VectorXd inverse_gompertz(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative_gompertz(const Eigen::VectorXd& eta) const ;

  virtual Eigen::VectorXd inverse_gumbel(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative_gumbel(const Eigen::VectorXd& eta) const ;

  // List GLMadj(std::string response,
  //             StringVector explanatory_complete,
  //             StringVector explanatory_proportional,
  //             std::string distribution,
  //             SEXP categories_order,
  //             DataFrame dataframe);
};

#endif
