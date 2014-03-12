# NOTE: Just a basic sanity check for now, do something like in tabletablr later
library(testthat)

set.seed(3438925)
test.data <- data.frame(a=runif(10000), b=runif(10000), c=runif(10000))
test.query <- data.frame(a=runif(3), b=runif(3), c=runif(3))

# computed from a reference implementation
expected.index <- c(4315, 5907, 8497, 7656, 8060, 4361, 1294, 7269, 1500, 8445, 9070, 972, 5718, 7803, 2928, 1302, 9730, 4643, 8245, 8414, 7420, 7813, 8977, 2735, 5832, 2813, 8098, 1948, 7822, 2930)
#expected.distance <- c(0.023777863, 0.02939512, 0.03467671, 0.03632639, 0.04516371, 0.05262052, 0.05309608, 0.05621594, 0.05970151, 0.06230218, 0.012517320, 0.02543950, 0.03570553, 0.03713558, 0.04281672, 0.04307229, 0.04547842, 0.06392136, 0.06534280, 0.06538633, 0.008412378, 0.03306358, 0.03375904, 0.03614532, 0.04117221, 0.04257510, 0.04294894, 0.04683583, 0.04753219, 0.05025776)
expected.distance <- c(0.023777863, 0.029395117, 0.034676710, 0.036326386, 0.045163708, 0.052620522, 0.053096077, 0.056215940, 0.059701512, 0.062302181, 0.012517320, 0.025439502, 0.035705529, 0.037135576, 0.042816720, 0.043072293, 0.045478418, 0.063921361, 0.065342797, 0.065386328, 0.008412378, 0.033063581, 0.033759037, 0.036145323, 0.041172207, 0.042575103, 0.042948938, 0.046835826, 0.047532191, 0.050257756)

actual.observed <- knn.lookup.rows(flanner(test.data), test.query, 10)

expect_equivalent(actual.observed, expected.index)
expect_equivalent(sqrt(attr(actual.observed, "distance")), expected.distance)

actual.observed <- knn.lookup.rows(flanner(test.data), test.query, 10, square.distance=FALSE)
expect_equivalent(actual.observed, expected.index)
expect_equivalent(attr(actual.observed, "distance"), expected.distance)

test.query2 <- data.frame(test.query, d=1:3)
actually.observed <- knn.lookup(flanner(test.data), test.query2, 10, distance.name="ddd")
expect_equivalent(actually.observed$d, rep(1:3, each=10))
expect_equivalent(actually.observed[,1:3], test.data[expected.index[1:30],])
expect_equivalent(sqrt(actually.observed$ddd), expected.distance)
expect_true(all(colnames(actually.observed) %in% c(colnames(test.data), colnames(test.query2), "ddd")))
