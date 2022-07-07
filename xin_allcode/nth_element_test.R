library(rbenchmark)
library(Rcpp)
sourceCpp('C:\\Users\\Xin Wang\\Desktop\\2Rcpp.cpp')
set.seed(10)



z = rnorm(100000)
x = rnorm(100)

n <- 25000




# check that stl_partial_sort is equal to nth_partial_sort
stopifnot(all.equal(stl_partial_sort(x, 50)[1:50], 
                    nth_partial_sort(x, 50)[1:50]))
system.time(sort(x))
# benchmark stl_partial_sort, nth_element_sort, and sort
benchmark(stl_partial_sort(z, n),
          nth_partial_sort(z, n),
          stl_nth_element(z, n),
          sort(z, partial=1:n),
          order(z,partial=1:n),
          order="relative")[,1:4]
stl_nth_element(x,10)
