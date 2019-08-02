# include <RcppEigen.h>
# include <Rcpp.h>

// [[Rcpp::depends(RcppEigen)]]

Eigen::MatrixXd makeBc (double A, double Q, double rho) {
    Eigen::MatrixXd B(14, 14);
    B.col(0) << 1- A, A, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    B.col(1) << 0, 0, (1 - rho) * A, 0, 0, 0, 0, (1 - rho) * (1 - A), 0, 0, 0, 0, rho, 0;
    B.col(2) << 0, 0, (1 - rho) * A, 0, 0, 0, 0, (1 - rho) * (1 - A), 0, 0, 0, 0, rho, 0;
    B.col(3) << 0, 0, (1 - rho) * A, 0, 0, 0, 0, (1 - rho) * (1 - A), 0, 0, 0, 0, rho, 0;
    B.col(4) << 0, 0, (1 - rho) * A, 0, 0, 0, 0, (1 - rho) * (1 - A), 0, 0, 0, 0, rho, 0;
    B.col(5) << 0, 0, (1 - rho) * A, 0, 0, 0, 0, (1 - rho) * (1 - A), 0, 0, 0, 0, rho, 0;
    B.col(6) << 0, 0, (1 - rho) * A, 0, 0, 0, 0, (1 - rho) * (1 - A), 0, 0, 0, 0, rho, 0;
    B.col(7) << (1- A) * (1 - Q), A * (1 - Q), 0, A * Q, 0, 0, 0, 0, (1 - A) * Q, 0, 0, 0, 0, 0;
    B.col(8) << (1- A) * (1 - Q), A * (1 - Q), 0, 0, A * Q, 0, 0, 0, 0, (1 - A) * Q, 0, 0, 0, 0;
    B.col(9) << (1- A) * (1 - Q), A * (1 - Q), 0, 0, 0, A * Q, 0, 0, 0, 0, (1 - A) * Q, 0, 0, 0;
    B.col(10) << (1- A) * (1 - Q), A * (1 - Q), 0, 0, 0, 0, A * Q, 0, 0, 0, 0, (1 - A) * Q, 0, 0;
    B.col(11) << (1- A) * (1 - Q), A * (1 - Q), 0, 0, 0, 0, A * Q, 0, 0, 0, 0, (1 - A) * Q, 0, 0;
    B.col(12) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;
    B.col(13) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    return B;
}

Eigen::MatrixXd makeVc(double Q, double rho) {
    Eigen::MatrixXd V = makeBc(0, Q, rho);
    return V;
}

Eigen::MatrixXd makeWc(double Q, double rho) {
    Eigen::MatrixXd V = makeVc(Q, rho);
    Eigen::MatrixXd B = makeBc(0.1, Q, rho);
    Eigen::MatrixXd W = (B - V) / 0.1;
    return W;
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
SEXP simA1(const double Q,
           const Eigen::Map<Eigen::VectorXd> D,
           const Eigen::Map<Eigen::VectorXd> X,
           const Eigen::Map<Eigen::VectorXd> Tx,
           const Eigen::Map<Eigen::VectorXd> inputA) {
    Eigen::VectorXd A(X.size());
    A(0) = inputA(0);
    Eigen::MatrixXd B = makeBc(inputA(0), Q, Tx(0));
    Eigen::VectorXd Y = eigenY(B);
    for(int i = 1; i < X.size(); i++) {
        Eigen::MatrixXd V = makeVc(Q, Tx(i - 1));
        Eigen::MatrixXd W = makeWc(Q, Tx(i - 1));
        double numer = X[i] - D.transpose() * V * Y;
        double denom = D.transpose() * W * Y;
        A[i] = (numer / denom);
        B = makeBc(A[i], Q, Tx(i));
        Y = B * Y;
    }
    return Rcpp::wrap(A);
}

// [[Rcpp::export]]
SEXP simA2(const double Q,
           const Eigen::Map<Eigen::VectorXd> D,
           const Eigen::Map<Eigen::VectorXd> X,
           const Eigen::Map<Eigen::VectorXd> Tx,
           const Eigen::Map<Eigen::VectorXd> inputA) {

    // Equilibrium associated with first AR
    Eigen::MatrixXd B = makeBc(inputA(0), Q, Tx(0));
    Eigen::VectorXd Y = eigenY(B);


    B = makeBc(inputA(1), Q, Tx(1));
    Y = B * Y;

    Eigen::VectorXd A(X.size() - 1);
    A(0) = inputA(1);
    for(int i = 2; i < X.size(); i++) {
        Eigen::MatrixXd V1 = makeVc(Q, Tx(i - 1));
        Eigen::MatrixXd W = makeWc(Q, Tx(i - 1));
        Eigen::MatrixXd V2 = makeVc(Q, Tx(i));
        double numer = X[i] - D.transpose() * V2 * V1 * Y;
        double denom = D.transpose() * V2 * W * Y;
        A[i - 1] = (numer / denom);
        B = makeBc(A[i - 1], Q, Tx(i - 1));
        Y = B * Y;
    }
    return Rcpp::wrap(A);
}
