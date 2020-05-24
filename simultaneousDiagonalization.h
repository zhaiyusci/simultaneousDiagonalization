#ifndef SIMULTANEOUSDIAGONALIZATION_H_
#define SIMULTANEOUSDIAGONALIZATION_H_
#include<initializer_list>
#include<vector>
#include<eigen3/Eigen/Dense>
class simultaneousDiagonalization{
  public:
    simultaneousDiagonalization(
        std::vector<Eigen::MatrixXd> matrices,
        double eps = 1.e-8);
    inline Eigen::VectorXd eigenvalues(int i) const {
      return matrices_[i].diagonal(); 
    }
    inline Eigen::MatrixXd eigenvectors() const {
      return eigenvectors_; // all matrices share the same eigen vectors
    }
  private:
    std::vector<Eigen::MatrixXd> matrices_;
    Eigen::MatrixXd eigenvectors_;
    int r_;
    int c_;
    static double off_(const Eigen::MatrixXd & A);
    double offsum_() const;
    Eigen::Matrix2d G_(int i, int j) const;
    Eigen::MatrixXd R_(int i, int j) const;
    void compute_(double eps);
};
#endif // SIMULTANEOUSDIAGONALIZATION_H_
