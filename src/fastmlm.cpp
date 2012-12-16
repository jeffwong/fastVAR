// [[Rcpp::depends(RcppEigen)]]

#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::ArrayXd;
using Eigen::Lower;
    
typedef Map<MatrixXd> MapMatd;

inline ArrayXd Dplus(const ArrayXd& d) {
    ArrayXd di(d.size());
    double comp(0.000000001);
    for (int j = 0; j < d.size(); ++j) di[j] = (d[j] < comp) ? 0. : 1./d[j];
    return di;
}

inline MatrixXd AtA(const MapMatd& A) {
    int n(A.cols());
    return MatrixXd(n,n).setZero().selfadjointView<Lower>()
    .rankUpdate(A.adjoint());
}

// [[Rcpp::export]]
List fast_mlm(NumericMatrix Xr, NumericMatrix Yr) {
    const MapMatd  X(as<MapMatd >(Xr));
    const MapMatd  Y(as<MapMatd >(Yr));
    const Eigen::SelfAdjointEigenSolver<MatrixXd> VLV(AtA(X));
    const ArrayXd Dp(Dplus(VLV.eigenvalues()).sqrt());
    const int r((Dp > 0).count());
    const MatrixXd VDp(VLV.eigenvectors() * Dp.matrix().asDiagonal());
    const MatrixXd betahat(VDp * VDp.adjoint() * X.adjoint() * Y);
    return List::create(Named("coefficients") = betahat);
}
