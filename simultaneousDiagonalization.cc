#include"simultaneousDiagonalization.h"
#include<iostream>

simultaneousDiagonalization::simultaneousDiagonalization(
    std::initializer_list<Eigen::MatrixXd> matrices, double eps ): 
  matrices_(matrices), 
  r_ ( matrices_[0].rows() ), c_ ( matrices_[0].cols() ) 
{
  if(r_ != c_){
    throw std::runtime_error("Sizes inconsist.");
  }
  for(auto && i : matrices_){
    if(i.rows() != r_ || i.cols() != c_){
      throw std::runtime_error("Sizes inconsist.");
    }
  }
  eigenvectors_.setIdentity(r_,c_);
  compute_(eps); // solve the problem on construction...
}

double simultaneousDiagonalization::off_(const Eigen::MatrixXd & A){
  double res=0.0;
  for(int i=0; i!=A.rows(); ++i){
    for(int j=i; j!=A.rows(); ++j){
      if(i != j) res += pow(A(i,j),2);
    }
  }
  return res;
}

double simultaneousDiagonalization::offsum_()const{
  double res=0.0;
  for(auto && i : this -> matrices_){
    res += off_(i);
  }
  return res;
}

Eigen::Matrix2d simultaneousDiagonalization::G_(int i, int j) const {
  Eigen::Matrix2d res=Eigen::Matrix2d::Zero();
  for (auto && A: this -> matrices_){
    Eigen::Vector2d h;
    h << A(i,i)-A(j,j), A(i,j)+A(j,i);
#ifndef NDEBUG
    std::cout << "h   " << std::endl << h << std::endl;
#endif
    res += h*h.adjoint();
  }
  return res;
}

Eigen::MatrixXd simultaneousDiagonalization::R_(int i, int j) const{
  Eigen::MatrixXd res = Eigen::MatrixXd::Identity(r_, c_);
  Eigen::Vector2d xy = 
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d>(G_(i,j))
    .eigenvectors().col(1);
  double x = xy[0];
  double y = xy[1];
  double c = sqrt((x+1.)/2.);
  double s = y/sqrt(2.*(x+1.));
  res(i,i) = res(j,j) = c;
  res(i,j) = -s; 
  res(j,i) =  s;
  return res;
}

void simultaneousDiagonalization::compute_(double eps){
  for(int iter = 0; iter < 1000; ++iter){
    if(offsum_()>eps){
      for(int i=0; i<r_; ++i){
        for(int j=i+1; j<r_; ++j){
          Eigen::MatrixXd Rmat = R_(i,j);
          for(auto && A : matrices_){
            A= Rmat.adjoint() * A * Rmat;
          }
          eigenvectors_ = eigenvectors_ * Rmat;
        }
      }
    }else{
      std:: cout << "CONVERGED." <<std::endl;
      return;
    }
  }
  for(int i = 0; i != 3; ++i) std::cout 
    << "-!- ### NOT CONVERGED WITHIN 1000 STEPS ### -!-" <<std::endl;
  return;
}
