#include"simultaneousDiagonalization.h"
#include<iostream>

#define N 5
int main(){
  // First, build a system that are Self-adjoint and get its eigenvectors
  Eigen::MatrixXd H(N,N);
  for(int i =0; i<N; ++i){
    for(int j =i; j <N; ++j){
      H(i,j)=H(j,i)=sin(i*0.1*M_PI)*exp(j*0.1);
    }
  }
  Eigen::MatrixXd eigenvectors = 
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(H)
    .eigenvectors();
  std::cout << "eigenvectors"<< std::endl;
  std::cout << eigenvectors << std::endl << std::endl;

  // Second, build some diagonal matrices that hold the eigenvectors
  Eigen::MatrixXd A(N,N);
  A.setIdentity(N,N);
  for(int i =0; i<N; ++i){
    A(i,i)=1.*i;
  }

  std::cout << "A eigen values " << A.diagonal().adjoint() << std::endl ;
  std::cout << std::endl;

  // Back to full matrices

  A=eigenvectors*A*eigenvectors.adjoint();

  simultaneousDiagonalization sd({A});
  std::cout << "sd.eigenvectors()"<< std::endl<< std::endl;
  std::cout << sd.eigenvectors()<< std::endl<< std::endl;
  std::cout << "sd.eigenvalues()"<< std::endl<< std::endl;
  std::cout << "A eigen values "<< sd.eigenvalues(0).adjoint()<< std::endl;

  return 0;
}

