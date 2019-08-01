library(Rcpp)
library(RcppEigen)
sourceCpp("src/forward_sim.cpp")

PAR = list()
r = 1 / 200
dt = 10
PAR$A = 0.1
PAR$Q = exp(-r *dt)
PAR$In = 5
PAR$Cn = 2
PAR$rho = 0.1


Bfn = makeBdrugs_age
Y = findYeq(PAR, Bfn)
D = makeD(PAR)
W = makeW(PAR, Bfn)
V = makeV(PAR, Bfn)
X = rnorm(18, 1, 0.05) * seq(0.5, 0.7, length.out = 18)
inputA = optim(fn = fn2, par = c(0, 0.1), PAR = PAR, Bfn = Bfn, X = X)$par

pct <- proc.time()
n = 4000
Amat1 = matrix(nrow = 17, ncol = n)
for(i in 1:n) {
    Amat1[, i] = simA2(V, W, D, X, inputA)
}
(proc.time() - pct)

pct <- proc.time()
Amat2 = matrix(nrow = 17, ncol = n)
for(i in 1:n) {
    Amat2[, i] = simAR(X, PAR, makeBdrugs_age, cpp = F)$A
}
(proc.time() - pct)

sim_forward(W, V, A, X)


etest <- cxxfunction(signature(tm="NumericMatrix",
                               tm2="NumericMatrix"),
                     plugin="RcppEigen",
                     body="
NumericMatrix tm22(tm2);
NumericMatrix tmm(tm);

const Eigen::Map<Eigen::MatrixXd> ttm(as<Eigen::Map<Eigen::MatrixXd> >(tmm));
const Eigen::Map<Eigen::MatrixXd> ttm2(as<Eigen::Map<Eigen::MatrixXd> >(tm22));

Eigen::MatrixXd prod = ttm*ttm2;
return(wrap(prod));
                 ")

set.seed(123)
M1 <- matrix(sample(1e3),ncol=50)
M2 <- matrix(sample(1e3),nrow=50)
etest(M1, M2)
