// #include <Rcpp.h>
// // #include <RcppEigen.h>
// using namespace Rcpp;
//
// // using Eigen::VectorXd;
// // using Eigen::MatrixXd;
// // using Eigen::Map;
// // using Eigen::SelfAdjointEigenSolver;
//
// // [[Rcpp::export]]
// NumericMatrix matrixMult(NumericMatrix A, NumericMatrix B) {
//     int n = A.nrow();
//     NumericMatrix C(n, n);
//     for(int i = 0; i < n; i++) {
//         for(int j = 0; j < n; j++) {
//             double sum = 0;
//             for(int k = 0; k < n; k++) {
//                 sum += A(i, k) * B(k, j);
//             }
//             C(i, j) = sum;
//         }
//     }
//     return C;
// }
//
// // [[Rcpp::export]]
// NumericMatrix scalarMult(double A, NumericMatrix B) {
//     int n = B.nrow();
//     NumericMatrix C(n, n);
//     for(int i = 0; i < n; i++) {
//         for(int j = 0; j < n; j++) {
//             C(i, j) = A * B(i, j);
//         }
//     }
//     return C;
// }
//
// // [[Rcpp::export]]
// NumericMatrix matrixAdd(NumericMatrix A, NumericMatrix B) {
//     int n = A.nrow();
//     NumericMatrix C(n, n);
//     for(int i = 0; i < n; i++) {
//         for(int j = 0; j < n; j++) {
//             C(i, j) = A(i, j) + B(i, j);
//         }
//     }
//     return C;
// }
//
// //
// // // [[Rcpp::export]]
// // MatrixXd makeB(MatrixXd W, MatrixXd V, double A) {
// //     MatrixXd B = A * W + V;
// //     return B;
// // }
// //
// // // [[Rcpp::export]]
// // VectorXd eigenY(Map<MatrixXd> W, Map<MatrixXd> V, double A) {
// //     MatrixXd B = makeB(W, V, A);
// //     SelfAdjointEigenSolver<MatrixXd> es(B);
// //     MatrixXd ev = es.eigenvectors();
// //     VectorXd first = ev.col(1) / ev.col(1).sum();
// //     return ev.col(1);
// // }
//
// // [[Rcpp::export]]
// NumericVector sim_forward(NumericMatrix W, NumericMatrix V, double A0, double A1,
//                           NumericMatrix D, NumericVector X, NumericVector Y) {
//     int n = X.size();
//     NumericVector Avec(n);
//     Avec(0) = A0;
//     NumericMatrix B = matrixAdd(scalarMult(A1, W), V);
//     Avec(1) = A1;
//     Y = matrixMult(B, Y);
//     for(int i = 2; i < (n + 1); i++) {
//         NumericMatrix VY = matrixMult(V, Y);
//         NumericMatrix VVY = matrixMult(V, VY);
//         NumericMatrix WY = matrixMult(W, Y);
//         NumericMatrix VWY = matrixMult(VWY);
//         double A = (X(i) - matrixMult(D, VVY) / matrixMult(D, VWY);
//         Avec(i) = A;
//         B = scalarMult(A, W) + V;
//         Y = matrixMult(B, Y);
//     }
//     return Avec;
// }
//
// // // [[Rcpp::export]]
// // VectorXd forwardA(Map<MatrixXd> W, Map<MatrixXd> V, double startA, Map<VectorXd> X) {
// //     VectorXd Avec;
// //     for(int i = 0; i < X.size(); ++i) {
// //         Avec(i) = X(i);
// //     }
// //     return(Avec);
// // }
//
