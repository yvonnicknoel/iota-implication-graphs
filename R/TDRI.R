tdri = read.csv2("TDRI.csv")
dim(tdri)
# [1] 1803   56

stages = c("Pre-Operational","Primary","Concrete","Abstract","Formal","Systematic","Meta-systematic")
Q = matrix(0, nrow=ncol(tdri), ncol=length(stages))

# Eight items for each stage
Q[1:8,1] = 1
Q[9:16,2] = 1
Q[17:24,3] = 1
Q[25:32,4] = 1
Q[33:40,5] = 1
Q[41:48,6] = 1
Q[49:56,7] = 1
colnames(Q) = stages

source("iota17.R")
Dt = DiG$new(tdri)

# Examine cross-tables
# D1$crosstab_report(output_file = "crosstab_report.tex", title = "Binary Variable Cross-Tabulation Report")

# Inference: compute logOR, Asymmetry, iota, p-values
Dt$compute(sym.dist="hybrid")
Dt$compute(sym.dist="minusmin")
Dt$compute(sym.dist="minshift")

# Artficially relocate two extreme items, only for plotting purposes
Dt$coordinates[5,2] = -2
Dt$coordinates[8,2] = 1.5

# SBM cluster IDs are arbitrary; remap them so that cluster k contains the items
# theoretically belonging to stage k (cluster 1 = items 1-8 = Pre-Operational, ...,
# cluster 7 = items 49-56 = Meta-systematic). This makes cluster_labels = stages
# directly usable and the legend read in theoretical order. The cached metagraph
# (adjacency + centroids) is also permuted so that meta-arrows track the new IDs.
remap_clusters_to_stages = function(Dt) {
  remap = integer(length(stages))
  for(k in seq_along(stages)) {
    q_items = ((k - 1) * 8 + 1):(k * 8)
    sbm_id  = as.integer(names(sort(table(Dt$communities[q_items]), decreasing = TRUE))[1])
    remap[sbm_id] = k
  }
  Dt$communities = remap[Dt$communities]
  inv = order(remap)
  if(!is.null(Dt$metagraph[["adjacency"]])) {
    A = Dt$metagraph[["adjacency"]][inv, inv]
    rownames(A) = colnames(A) = paste0("Cluster", seq_along(remap))
    Dt$metagraph[["adjacency"]] = A
  }
  if(!is.null(Dt$metagraph[["centroids"]])) {
    Dt$metagraph[["centroids"]] = Dt$metagraph[["centroids"]][inv, , drop=FALSE]
  }
}
remap_clusters_to_stages(Dt)

Dt$plot.tikz(file="TDRI-full-meta-graph-minshift.tex", plot_communities=TRUE, community_type="ellipses", plot_metagraph=TRUE, meta.bend_right=1:7, meta.bend_angle=60, only_directed=FALSE, only_symmetric=FALSE, bend_angle=30, dilate_cluster=0.8, scale_factor=5, cluster_labels=stages, ylim=c(-5,5))

#------------------------------- Delta - logOR decomposition ---------------------------
iota_star = Dt$IOTA_STAR

p <- nrow(iota_star)

# Compute min and max matrices
min_mat <- pmin(iota_star, t(iota_star))
max_mat <- pmax(iota_star, t(iota_star))

# Extract pure components
# |Delta| = (max - min) / 2
# log(OR) = (max + min) / 2
delta_mat <- (max_mat - min_mat) / 2  # This is |Delta|, always positive
logOR_mat <- (max_mat + min_mat) / 2  # This is log(OR)

# Set diagonal to 0 for distance matrix
diag(delta_mat) <- 0
diag(logOR_mat) <- 0

# 1. Use |Delta| matrix for difficulty reconstruction via 1D MDS
delta_dist <- as.dist(delta_mat)
mds_delta <- cmdscale(delta_dist, k = 1, eig = TRUE)
difficulty <- mds_delta$points[, 1]

# 2. Use log(OR) for scale detection
# Compute mean log(OR) profile for each item
mean_logOR <- rowMeans(logOR_mat)

max_logOR <- max(logOR_mat)
logOR_dist <- as.dist(max_logOR - logOR_mat)
mds_logOR <- cmdscale(logOR_dist, k = 1, eig = TRUE)

association <- mds_logOR$points[, 1]

# 3. Plot results
plot(difficulty, association, pch = 19, cex = 1.2, xlab = "Difficulty (from |Δ|)", ylab = "Association (from log OR)", main = "Decomposition")
text(difficulty, association, labels = 1:p, pos = 3, cex = 0.6)

Dt$coordinates = cbind(difficulty,association)
Dt$plot.tikz(file="TDRI-full-meta-graph-decomposition.tex", plot_communities=TRUE, community_type="ellipses", plot_metagraph=TRUE, meta.bend_right=5:7, meta.bend_angle=60, only_directed=FALSE, only_symmetric=FALSE, bend_angle=30, dilate_cluster=0.8, scale_factor=5, cluster_labels=stages)

# Les log-OR diminuent avec les stades
for(i in 1:7) {
  cat(stages[i],"\n")
  start = 1+(i-1)*8
  end = start+7
  print(round(Dt$logOR[start:end,start:end],3))
}

# A corriger : les logOR diagonaux ne sont pas mis à NA

#------------- Eigenvalues
x11()
plot(D1$eig[1:40],type="h",xlab="Dimension",ylab="Eigenvalue")
dev.copy2pdf(file="TDRI-eigenvalues.pdf")

# Approche simple : Strongly connected components
g1 = graph_from_adjacency_matrix(D1$G, mode="directed", weighted=NULL)
scc = components(g1, mode = "strong")
scc$membership

table(scc$membership)
#  1  2  3  4 
# 16 16  8 16                  # Les 4 et 5 sont fusionnés (sans doute à cause des deux liens entre eux)


#------------------------------- Degree-corrected Stochastic Block Model (DC-SBM) --------------------------------
library(Matrix)
library(greed)

# A = your 0/1 adjacency (N x N), directed, no self-loops
A = as.matrix(D1$G)
A = Matrix(A, sparse = TRUE)

fit = greed(A, model = DcSbm(type = "directed"))  # DC-SBM, directed
z   = clustering(fit)                              # block labels (length N)
Khat= K(fit)                                       # chosen # of blocks (ICL)
par = coef(fit)                                    # parameters (MAP)
Theta = par$thetakl        # K x K normalized block intensities
gamma_out = par$gammaout   # node out-degree factors
gamma_in  = par$gammain    # node in-degree factors

#------------------------------- Max-Dicut approach --------------------------------
# On cherche la partition en deux qui maximise la difference de flux sortant / entrant
max_dicut = function(A, idx) {
  B = A[idx, idx]; diag(B) = 0
  p = nrow(B); best = -Inf; bestS = NULL
  for (mask in 1:(2^p - 2)) {                     # exclude empty and full sets
    S = which(as.logical(intToBits(mask)[1:p]))
    T = setdiff(seq_len(p), S)
    val = sum(B[S, T]) - sum(B[T, S])
    if (val > best) { best = val; bestS = S }
  }
  # shrink T to a minimal set that preserves the score
  S = bestS; T = setdiff(seq_len(p), S)
  repeat {
    drop = T[sapply(T, function(t) {
      T2 = setdiff(T, t)
      (sum(B[S, T2]) - sum(B[T2, S])) == best
    })]
    if (length(drop) == 0) break
    T = setdiff(T, drop)
  }
  list(score = best, S = idx[S], T = idx[T])
}

# usage:
A = D1$G
res = max_dicut(A, 17:32)
res
# $score
# [1] 10

# $S
# [1] 25 26 27 29 31

# $T
# [1] 17 18 19 20

library(igraph)
max_dicut_score = function(B){
  p = nrow(B); best = -Inf
  for (mask in 1:(2^p - 2)) {
    S = which(as.logical(intToBits(mask)[1:p])); T = setdiff(seq_len(p), S)
    best = max(best, sum(B[S,T]) - sum(B[T,S]))
  }
  best
}
B = A[17:32,17:32]; diag(B) = 0
obs = max_dicut_score(B)

g0 = graph_from_adjacency_matrix(B, mode="directed")
od = degree(g0, mode="out"); id = degree(g0, mode="in")
null = replicate(1000, {
  g = sample_degseq(od, id, method="simple")
  max_dicut_score(as.matrix(as_adj(g)))
})
pval = mean(null >= obs)
pval

#------------------------------- Stochastic Block Model ------------------------------

library(blockmodels)                   # fast variational SBM

#-- Graphe complet
Tr = D1$G
mod = BM_bernoulli("SBM", Tr, plotting="", verbosity=0, explore_min=2, explore_max=10)
mod$estimate()                        # fits K = 2...10 blocks
Kbest = which.max(mod$ICL)            # choose by Integrated Completed-Likelihood
probs   = mod$memberships[[Kbest]]$Z     # hard assignment of items to blocks
cl = apply(probs, 1, which.max)
names(cl) = colnames(tdri)
cl # 5 clusters
#  i1  i2  i3  i4  i5  i6  i7  i8  i9 i10 i11 i12 i13 i14 i15 i16 i17 i18 i19 i20 
#   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   2   2   2   2 
# i21 i22 i23 i24 i25 i26 i27 i28 i29 i30 i31 i32 i33 i34 i35 i36 i37 i38 i39 i40 
#   2   2   2   2   2   2   2   2   2   2   2   2   4   4   4   4   4   4   4   4 
# i41 i42 i43 i44 i45 i46 i47 i48 i49 i50 i51 i52 i53 i54 i55 i56 
#   1   1   1   1   1   1   1   1   5   5   5   5   5   5   5   5 

#-- Graphe univoque seulement
Tr = D1$G.uni
mod = BM_bernoulli("SBM", Tr, plotting="", verbosity=0, explore_min=2, explore_max=10)
mod$estimate()                        # fits K = 4...7 blocks
Kbest = which.max(mod$ICL)            # choose by Integrated Completed-Likelihood
probs   = mod$memberships[[Kbest]]$Z     # hard assignment of items to blocks
cl = apply(probs, 1, which.max)
names(cl) = colnames(tdri)
cl # 5 clusters aussi
#  i1  i2  i3  i4  i5  i6  i7  i8  i9 i10 i11 i12 i13 i14 i15 i16 i17 i18 i19 i20 
#   2   2   2   2   2   2   2   2   2   2   2   2   2   2   2   2   5   5   5   5 
# i21 i22 i23 i24 i25 i26 i27 i28 i29 i30 i31 i32 i33 i34 i35 i36 i37 i38 i39 i40 
#   5   5   5   5   5   5   5   5   5   5   5   5   3   3   3   3   3   3   3   3 
# i41 i42 i43 i44 i45 i46 i47 i48 i49 i50 i51 i52 i53 i54 i55 i56 
#   1   1   1   1   1   1   1   1   4   4   4   4   4   4   4   4 

#-- Graphe biunivoque seulement
Tr = D1$G.bi
mod = BM_bernoulli("SBM", Tr, plotting="", verbosity=0, explore_min=2, explore_max=10)
mod$estimate()                        # fits K = 4...7 blocks
Kbest = which.max(mod$ICL)            # choose by Integrated Completed-Likelihood
probs   = mod$memberships[[Kbest]]$Z     # hard assignment of items to blocks
cl = apply(probs, 1, which.max)
names(cl) = colnames(tdri)
cl # même résultat
#  i1  i2  i3  i4  i5  i6  i7  i8  i9 i10 i11 i12 i13 i14 i15 i16 i17 i18 i19 i20 
#   4   4   4   4   4   4   4   4   4   4   4   4   4   4   4   4   3   3   3   3 
# i21 i22 i23 i24 i25 i26 i27 i28 i29 i30 i31 i32 i33 i34 i35 i36 i37 i38 i39 i40 
#   3   3   3   3   3   3   3   3   3   3   3   3   2   2   2   2   2   2   2   2 
# i41 i42 i43 i44 i45 i46 i47 i48 i49 i50 i51 i52 i53 i54 i55 i56 
#   1   1   1   1   1   1   1   1   5   5   5   5   5   5   5   5 

#------------------------------- Latent Block Model ------------------------------

library(blockmodels)

#-- Graphe complet
Tr = D1$G
mod = BM_bernoulli("LBM", Tr, plotting="", verbosity=0, explore_min=2, explore_max=10)
mod$estimate()                        # fits K = 4...7 blocks
Kbest = which.max(mod$ICL)            # choose by Integrated Completed-Likelihood
probs   = mod$memberships[[Kbest]]$Z1 # Ou Z2 pour homogénéité par cibles
cl = apply(probs, 1, which.max)
names(cl) = colnames(tdri)
cl # 5 clusters
# Toujours les mêmes clusters

#-- Graphe univoque seulement
Tr = D1$G.uni
mod = BM_bernoulli("LBM", Tr, plotting="", verbosity=0, explore_min=2, explore_max=10)
mod$estimate()                        # fits K = 4...7 blocks
Kbest = which.max(mod$ICL)            # choose by Integrated Completed-Likelihood
probs   = mod$memberships[[Kbest]]$Z1
cl = apply(probs, 1, which.max)
names(cl) = colnames(tdri)
cl # 5 clusters aussi

#-- Sous-graphe unidirectionnel classe 2
sel2 = which(cl==4)
subg = D1$G.uni[sel2,sel2]

mod2 = BM_bernoulli("SBM", subg, plotting="", verbosity=0, explore_min=2, explore_max=10)
mod2$estimate()                        # fits K = 4...7 blocks
Kbest = which.max(mod2$ICL)            # choose by Integrated Completed-Likelihood
probs2   = mod2$memberships[[Kbest]]$Z     # hard assignment of items to blocks
cl2 = apply(probs2, 1, which.max)
names(cl2) = colnames(tdri)[sel2]
cl2
# i17 i18 i19 i20 i21 i22 i23 i24 i25 i26 i27 i28 i29 i30 i31 i32 
#   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 

#-- Graphe biunivoque seulement
Tr = D1$G.bi
mod = BM_bernoulli("LBM", Tr, plotting="", verbosity=0, explore_min=2, explore_max=10)
mod$estimate()                        # fits K = 4...7 blocks
Kbest = which.max(mod$ICL)            # choose by Integrated Completed-Likelihood
probs   = mod$memberships[[Kbest]]$Z2
cl = apply(probs, 1, which.max)
names(cl) = colnames(tdri)
cl # une seule classe, que ce soit en envoi ou en réception ??
#  i1  i2  i3  i4  i5  i6  i7  i8  i9 i10 i11 i12 i13 i14 i15 i16 i17 i18 i19 i20 
#   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
# i21 i22 i23 i24 i25 i26 i27 i28 i29 i30 i31 i32 i33 i34 i35 i36 i37 i38 i39 i40 
#   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
# i41 i42 i43 i44 i45 i46 i47 i48 i49 i50 i51 i52 i53 i54 i55 i56 
#   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 


#-- Sous-graphe bidirectionnel classe 2
sel2 = which(cl==2)
subg = D1$G.bi[sel2,sel2]

mod2 = BM_bernoulli("SBM", subg, plotting="", verbosity=0, explore_min=2, explore_max=10)
mod2$estimate()                        # fits K = 4...7 blocks
Kbest = which.max(mod2$ICL)            # choose by Integrated Completed-Likelihood
probs2   = mod2$memberships[[Kbest]]$Z     # hard assignment of items to blocks
cl2 = apply(probs2, 1, which.max)
names(cl2) = colnames(tdri)[sel2]
cl2
# i17 i18 i19 i20 i21 i22 i23 i24 i25 i26 i27 i28 i29 i30 i31 i32 
#   2   1   1   2   2   2   2   2   1   1   2   2   1   2   2   2 

#-- Sous-graphe unidirectionnel classe 2
sel2 = which(cl==2)
subg = D1$G.uni[sel2,sel2]

mod2 = BM_bernoulli("SBM", subg, plotting="", verbosity=0, explore_min=2, explore_max=10)
mod2$estimate()                        # fits K = 4...7 blocks
Kbest = which.max(mod2$ICL)            # choose by Integrated Completed-Likelihood
probs2   = mod2$memberships[[Kbest]]$Z     # hard assignment of items to blocks
cl2 = apply(probs2, 1, which.max)
names(cl2) = colnames(tdri)[sel2]
cl2

# On essaie les LBM qui classient les items par similarité de profil envoyeur ou receveur, séparés (deux classifications, ligne et colonne)
mod3 = BM_bernoulli("LBM", subg, plotting="", verbosity=0, explore_min=2, explore_max=10)
mod3$estimate()                        
Kbest = which.max(mod3$ICL)            # choose by Integrated Completed-Likelihood
probs3   = mod3$memberships[[Kbest]]$Z1
cl3 = apply(probs3, 1, which.max)
names(cl3) = colnames(tdri)[sel2]

# Pas de distinction
cl3
# i17 i18 i19 i20 i21 i22 i23 i24 i25 i26 i27 i28 i29 i30 i31 i32 
#   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 

# Essai avec pondération 2 pour les liens bidirectionnels (un peu arbitraire)
Tr.w = Tr + D1$G.bi                      # This gives bidir weights of 2, and unidir weights of 1
mod.w = BM_poisson("SBM", Tr.w, plotting="", verbosity=0, explore_min=4, explore_max=10)
mod.w$estimate()                        # fits K = 4...7 blocks
Kbest = which.max(mod$ICL)            # choose by Integrated Completed-Likelihood
cl = mod$memberships[[Kbest]]$Z     # hard assignment of items to blocks

fit_nested = estimateSimpleSBM(Tr,
                                model = "bernoulli",
                                directed = TRUE,
                                nested = TRUE,
                                estimOptions = list(exploreMin = 2,
                                                    exploreMax = 3))  # sub-blocks per parent

# Essai de Walktrap
Tr = D1$G.bi
g = graph_from_adjacency_matrix(Tr, mode="undirected", weighted=NULL)
communities = cluster_walktrap(g, steps=10)
membership_vec = membership(communities)
max(membership_vec)
modularity(communities)

# Essai de Leiden
Tr = D1$G.bi
g = graph_from_adjacency_matrix(Tr, mode="undirected", weighted=NULL)
communities = cluster_leiden(g, objective_function="modularity")
membership_vec = membership(communities)
max(membership_vec)
modularity(communities)


# Exploration du Cluster 2
subg   = induced_subgraph(g, which(membership(communities)==2))  # cluster 2 only

subcl  = cluster_leiden(subg, objective_function = "modularity", resolution = 1.0)         # >1 ⇒ smaller groups
table(subcl$membership)
subcl$membership

subcl  = cluster_leiden(subg, objective_function = "modularity", resolution = 1.5)         # >1 ⇒ smaller groups
table(subcl$membership)

subcl  = cluster_leiden(subg, objective_function = "modularity", resolution = 2.0)         # >1 ⇒ smaller groups
table(subcl$membership)


# Essai de infomap
g = graph_from_adjacency_matrix(D1$G.bi, mode="undirected", weighted=NULL)
communities = cluster_infomap(g, nb.trials = 1000)
membership_vec = membership(communities)
n_communities = max(membership_vec)
modularity(communities)
[1] 0.7023028

membership_vec
#  i1  i2  i3  i4  i5  i6  i7  i8  i9 i10 i11 i12 i13 i14 i15 i16 i17 i18 i19 i20 
#   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   2   2   2   2 
# i21 i22 i23 i24 i25 i26 i27 i28 i29 i30 i31 i32 i33 i34 i35 i36 i37 i38 i39 i40 
#   2   2   2   2   2   2   2   2   2   2   2   2   3   3   3   3   3   3   3   3 
# i41 i42 i43 i44 i45 i46 i47 i48 i49 i50 i51 i52 i53 i54 i55 i56 
#   4   4   4   4   4   4   4   4   5   5   5   5   5   5   5   5 

# CLustering métrique par GMM 
D1$cluster("VVV")
D1$plot.tikz(file="TDRI.tex", color_univoque="black!80", color_biunivoque="gray!60", plot_ellipses=TRUE)
D1$plot.tikz(file="TDRI-biunivoque.tex", only_symmetric=TRUE)
D1$plot.tikz(file="TDRI-univoque.tex", only_directed=TRUE)

# Inférence iota IUT
D1$compute(inference="combined")
D1$plot.tikz(file="TDRI-comb.tex", plot_coefs=FALSE)
D1$plot.tikz(file="TDRI-biunivoque-comb.tex", only_symmetric=TRUE, plot_coefs=FALSE)
D1$plot.tikz(file="TDRI-univoque-comb.tex", only_directed=TRUE, plot_coefs=FALSE)

# Inférence iota bayésienne
D1$compute(inference="bayesian")
D1$plot.tikz(file="TDRI-bayes.tex", plot_coefs=FALSE)
D1$plot.tikz(file="TDRI-biunivoque-bayes.tex", only_symmetric=TRUE, plot_coefs=FALSE)
D1$plot.tikz(file="TDRI-univoque-bayes.tex", only_directed=TRUE, plot_coefs=FALSE)
