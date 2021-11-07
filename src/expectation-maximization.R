data = runif(10000)

MAX_ITERATIONS = 5000
CONVERGENCE_THRESHOLD = 1e-08
dimensions = 4

# modified sum only considers finite values
sum.finite = function (data) {
	sum(data[is.finite(data)])
}

# Initialize using K-Means
km = kmeans(data, dimensions)$cluster
means = sapply(1 : dimensions, function (i) {
  mean(data[km == i])
})
sigmas = sapply(1 : dimensions, function (i) {
  sd(data[km == i])
})
probs = sapply(1 : dimensions, function (i) {
  sum(km == i) / length(km)
})

# Initial value of expected log likelihood
log_likelihood = 0
log_likelihood[2] = sum(sapply(1 : dimensions, function (i) {
  sum.finite(log(probs[i]) + log(dnorm(data, means[i], sigmas[i])))
}))

index = 2


# Run EM Algorithm till covergence threshold or max iterations is reached
while (abs(log_likelihood[index] - log_likelihood[index-1]) >= CONVERGENCE_THRESHOLD && index <= MAX_ITERATIONS) {
	
	# Expectation step
	comps = sapply(1 : dimensions, function (i) {
		probs[i] * dnorm(data, means[i], sigmas[i])
	})
	comp_sum = rowSums(comps)
	ps = sapply(1 : dimensions, function (i) {
		comps[ ,i] / comp_sum
	})
	
	# Maximization step
	probs = sapply(1 : dimensions, function (i) {
		sum.finite(ps[ ,i]) / length(data)
	})
	means = sapply(1 : dimensions, function (i) {
		sum.finite(ps[ ,i] * data) / sum.finite(ps[, i])
	})
	sigmas = sapply(1 : dimensions, function (i) {
		sqrt(sum.finite(ps[ ,i] * (data - means[i]) ^ 2) / sum.finite(ps[, i]))
	})
	
	index = index + 1
	log_likelihood[index] = sum(log(comp_sum))
}

cat("Dimensions: ", dimensions, "\n")
cat("Iterations: ", index, "\n")
cat("Computed Cluster Means: ", means, "\n")
cat("Computed Cluster Standard Deviations: ", sigmas, "\n")

# Check final results using the mixtools library
library(mixtools)
em = normalmixEM (  x = data, 
                    k = dimensions,
                    mu = means,
                    sigma = sigmas,
                    maxit = MAX_ITERATIONS,
                    epsilon = CONVERGENCE_THRESHOLD
                )
cat("Library Cluster Means: ", em$mu, "\n")
cat("Library Cluster Standard Deviations: ", em$sigma, "\n")
