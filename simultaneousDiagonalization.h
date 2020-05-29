#ifndef SIMULTANEOUSDIAGONALIZATION_H_
#define SIMULTANEOUSDIAGONALIZATION_H_
#include<initializer_list>
#include<vector>
#include<eigen3/Eigen/Dense>
#include<eigen3/Eigen/Sparse>
class simultaneousDiagonalization{
  public:
    simultaneousDiagonalization(
        std::vector<Eigen::MatrixXd> matrices,
        double eps = 1.e-8);
    inline Eigen::VectorXd eigenvalues(int i) const {
      return matrices_[i].diagonal(); 
    }
    inline Eigen::MatrixXd remains(int i) const {
      return matrices_[i]; 
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
    // Eigen::SparseMatrix<double> R_(int i, int j) const;
    std::tuple<double, double> R_(int i, int j) const;
    void compute_(double eps);
    std::vector<double> scale_;
    void deal_fixed_point_();
    void remove_scaling_();
    void sweep_();
};
#endif // SIMULTANEOUSDIAGONALIZATION_H_
