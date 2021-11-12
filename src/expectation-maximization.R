# R program to simulate the d-dimensional Expectation Maximization Clustering Algorithm

# Team Members:

# Animesh Sharma                - as3592
# Yaniv Bronshtein              - yb262
# Wanying Mo                    - wm318
# Venkata Krishnam Raju Dalta   - vkd20
# Vipul Gharde                  - vig4
# Aditya Maheshwari             - am2971
# Toshitt Ahuja                 - ta498
# Fan Shen                      - fs470

# Generate random input data from a uniform distribution
data = runif(100)

MAX_ITERATIONS = 5000
CONVERGENCE_THRESHOLD = 1e-08
dimensions = 5

# Modified sum function only considers finite values
finite_sum = function (data) {
	return (sum(data[is.finite(data)]))
}

# Initialize parameters using K-Means
km = kmeans(data, dimensions)$cluster
means = sapply(1 : dimensions, function (i) {
	mean(data[km == i])
})
sigmas = sapply(1 : dimensions, function (i) {
	sd(data[km == i])
})
probabilities = sapply(1 : dimensions, function (i) {
	sum(km == i) / length(km)
})
print(probabilities)

# Initial value of expected log likelihood
log_likelihood = 0
log_likelihood[2] = sum(log(sapply(1 : dimensions, function (i) {
	finite_sum(probabilities[i] * dnorm(data, means[i], sigmas[i]))
})))

index = 2

# Run EM Algorithm till covergence threshold or max iterations is reached
while (abs(log_likelihood[index] - log_likelihood[index-1]) >= CONVERGENCE_THRESHOLD && index <= MAX_ITERATIONS) {
	
	# Expectation step
	comps = sapply(1 : dimensions, function (i) {
		probabilities[i] * dnorm(data, means[i], sigmas[i])
	})
	comp_sum = rowSums(comps)
	ps = sapply(1 : dimensions, function (i) {
		comps[ ,i] / comp_sum
	})
	
	# Maximization step
	# Update values of the parameters
	probabilities = sapply(1 : dimensions, function (i) {
		finite_sum(ps[ ,i]) / length(data)
	})
	means = sapply(1 : dimensions, function (i) {
		finite_sum(ps[ ,i] * data) / finite_sum(ps[, i])
	})
	sigmas = sapply(1 : dimensions, function (i) {
		sqrt(finite_sum(ps[ ,i] * (data - means[i]) ^ 2) / finite_sum(ps[, i]))
	})
	
	index = index + 1
	log_likelihood[index] = sum(log(comp_sum))
}

cat("Dimensions: ", dimensions, "\n")
cat("Iterations: ", index, "\n")
cat("Computed Cluster Means: ", means, "\n")
cat("Computed Cluster Standard Deviations: ", sigmas, "\n")

clusters = matrix(, nrow = length(data), ncol = dimensions)

for(i in 1:length(data)) {
    index <- arrayInd(which.max(ps[i,]), dimensions)
    clusters[[i, index]] = data[i]
}

print(ps)
print(rowSums(ps))

print("Clusters:\n")
print(clusters)

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



