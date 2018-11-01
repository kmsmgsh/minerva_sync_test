library(Rcpp)
incltxt <- '
int fibonacci(const int x) {
  if (x == 0) return(0);
  if (x == 1) return(1);
  return fibonacci(x - 1) + fibonacci(x - 2);
}'
cppFunction(incltxt)
fibRcpp <- cxxfunction(signature(xs="int"),
                       plugin="Rcpp",
                       incl=incltxt,
                       body='
                       int x = Rcpp::as<int>(xs);
                       return Rcpp::wrap( fibonacci(x) );')
library(inline)
fibR <- function(n) {
  if (n == 0) return(0)
  if (n == 1) return(1)
  return (fibR(n - 1) + fibR(n - 2))}
