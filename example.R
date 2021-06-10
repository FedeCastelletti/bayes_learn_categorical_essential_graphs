## Generate a DAG D and a DAG-dependent categorical dataset Y

source("gen_data.r")

q = 5

out_data = gen_dag_data(i = 1, n = 1000, q = q, p = 0.25)

Y = out_data$Y
D = out_data$D

head(Y)

plot(D)

## Load post_ess_graphs function and run it

source("mcmc_ess.r")

out_mcmc = post_ess_graphs(Y, m = 2*q, T = 200)


# Compute posterior probabilities of edge inclusion

post_probs = matrix(rowMeans(out_mcmc$chain.G), 5, 5)

post_probs

# Obtain an EG estimate as follows

median_model = as(round(post_probs), "graphNEL")

median_pr_model = dag2essgraph(median_model)

# Compare true and estimated EG

par(mfrow = (c(1,2)))

plot(dag2essgraph(D))
plot(median_pr_model)

