# include <RcppEigen.h>
# include <Rcpp.h>

// [[Rcpp::depends(RcppEigen)]]

Eigen::MatrixXd makeB(Eigen::MatrixXd V, Eigen::MatrixXd W, double A) {
    Eigen::MatrixXd B = A * W + V;
    return B;
}

Eigen::VectorXd eigenY(Eigen::MatrixXd B) {
    Eigen::EigenSolver<Eigen::MatrixXd> es(B);
    Eigen::VectorXd vals = es.eigenvalues().real();
    int j = 0;
    for(int i = 0; i < vals.size(); i++) {
        if(vals(i) == 1) {
            j = i;
        }
    }
    Eigen::VectorXd Y = es.eigenvectors().col(j).real();
    Y = Y / Y.sum();
    return Y;
}

// [[Rcpp::export]]
SEXP simA1(const Eigen::Map<Eigen::MatrixXd> V,
           const Eigen::Map<Eigen::MatrixXd> W,
           const Eigen::Map<Eigen::VectorXd> D,
           const Eigen::Map<Eigen::VectorXd> X,
           const Eigen::Map<Eigen::VectorXd> inputA) {
    Eigen::VectorXd A(X.size());
    A(0) = inputA(0);
    Eigen::MatrixXd B = makeB(V, W, inputA(0));
    Eigen::VectorXd Y = eigenY(B);
    for(int i = 1; i < X.size(); i++) {
        double numer = X[i] - D.transpose() * V * Y;
        double denom = D.transpose() * W * Y;
        A[i] = (numer / denom);
        B = makeB(V, W, A[i]);
        Y = B * Y;
    }
    return Rcpp::wrap(A);
}

// [[Rcpp::export]]
SEXP simA2(const Eigen::Map<Eigen::MatrixXd> V,
           const Eigen::Map<Eigen::MatrixXd> W,
           const Eigen::Map<Eigen::VectorXd> D,
           const Eigen::Map<Eigen::VectorXd> X,
           const Eigen::Map<Eigen::VectorXd> inputA) {

    Eigen::MatrixXd B = makeB(V, W, inputA(0));

    Eigen::VectorXd Y = eigenY(B);

    B = makeB(V, W, inputA(1));

    Y = B * Y;

    Eigen::VectorXd A(X.size() - 1);
    A(0) = inputA(1);
    for(int i = 2; i < X.size(); i++) {
        double numer = X[i] - D.transpose() * V * V * Y;
        double denom = D.transpose() * V * W * Y;
        A[i - 1] = (numer / denom);
        B = makeB(V, W, A[i - 1]);
        Y = B * Y;
    }
    return Rcpp::wrap(A);
}
