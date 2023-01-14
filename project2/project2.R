tridiagonal_matrix <- function(below, diagonal, over, func) {
  n <- length(diagonal)
  
  over[1] <- over[1] / diagonal[1]
  func[1] <- func[1] / diagonal[1]
  
  for(i in 2:(n - 1)) {
    over[i] <- over[i] / (diagonal[i] - below[i - 1] * over[i - 1])
    func[i] <- (func[i] - below[i - 1] * func[i - 1]) /
      (diagonal[i] - below[i - 1] * over[i - 1])
  }
  func[n] <- (func[n] - below[n - 1] * func[n - 1]) /
    (diagonal[n] - below[n - 1] * over[n - 1])

  x <- vector(length = n)
  x[n] <- func[n]
  for(i in (n - 1):1)
    x[i] <- func[i] - over[i] * x[i + 1]
  
  return(x)
}

f <- function(x){
  return (x*x + 3)
}

solve_equation <- function(n){
  h = 1 / n
  
  diag <- vector(length = n)
  below <- vector(length = n-1)
  over <- vector(length = n-1)
  func <- vector(length = n)
  
  diag_coof <- 1 - 2/h^2
  below_coof <- 1/h^2 + 1/h
  over_coof <- 1/h^2 - 1/h
  
  for (i in 1:(n-1)){
    func[i] <- f(i * h)
    diag[i] <- diag_coof
    below[i] <- below_coof
    over[i] <- over_coof
  }
  
  func[1] = func[1] - below_coof
  
  diag[n] <- diag_coof + over_coof*(-h)
  below[n-1] <- below_coof + over_coof
  func[n] <- f(1) - over_coof * 2 * h
  
  output <- tridiagonal_matrix(below, diag, over, func)
  output <- c(1, output)
  return (output)
}

plot_approximate_sollution <- function(eq, color="violet"){
  n = length(eq) - 1
  h = 1 / n
  a = 0
  b = 1
  x_vals = seq(a, b, h)
  
  plot(x_vals, eq, type="l", col=color, lwd="4",
       main=paste("Wykres funkcji u(x) na przedziale [0, 1] przy n =", n), xlab="x", ylab="y")
  legend(x="topright", legend=c("Rozwiązanie przybliżone"),
         col=c("violet"), lwd=c(4), text.font=4, bg='lightblue')
}

exact_u <- function(x){
  (24*(exp(1) - 1)/(5 * exp(1)))*x*exp(x) - 8*exp(x) + x^2 + 4*x + 9
}

compare_with_exact <- function(approximation, color="violet"){
  n = length(approximation) - 1
  h = 1 / n
  a = 0
  b = 1
  x_vals = seq(a, b, h)
  
  plot(x_vals, approximation, type="l", col=color, lwd="4",
       main="Wykres funkcji u(x)", xlab="x", ylab="y")
  curve(exact_u, from=a, to=b, lwd=2, add=TRUE)
  legend(x="topright", legend=c("Rozwiązanie dokładne",
                                paste("Rozwiazanie przyblizone dla n =",n)),
         col=c("black", "violet"), lwd=c(2, 4), text.font=4, bg='lightblue')
  
  exact_vals = exact_u(x_vals)
  error = abs(exact_vals - approximation)
  global_error <- max(error)
  return(global_error)
}

approx10 <- solve_equation(10)
approx50 <- solve_equation(50)

print(approx10)
print(approx50)

plot_approximate_sollution(approx10)
plot_approximate_sollution(approx50)

compare_with_exact(approx10)
compare_with_exact(approx50)

