source("iota17.R")
google.colors = c(blue="#4285f4", red="#db4437", yellow="#f4b400", green="#0f9d58", black="#333333")

#---------------------------------- Figure 2: Plot the implication graph for 7 items --------------------------------------
Dr = DiG$new(NULL)
Dr$generate_Rasch_scales(list(seq(-2,2,len=7)))
Dr$compute()

# Illustrating a simplified Rasch structure for 7 items
Dr$coordinates = cbind(seq(-2,2,len=7),rep(0,7))

# Plot (TikZ)
Dr$plot.tikz(filename="rasch7-tikz.tex", bend_right=1:7, bend_angle=30)

# The plot has then been manually tweaked

#----------------------- Figure 3: Delta vs. log-OR with varying latent variances and covariances ------------------------

pred_logOR = function(s1, s2, rho) {
  C = pi**2/3
  C*rho*s1*s2 / sqrt((s1**2 + C)*(s2**2+C))
}

plot_Delta_logOR <- function(Delta, logOR, scales,
                              col_within = c("#4285F4", "#DB4437"),
                              col_between = "grey50",
                              main = "Delta vs log-OR",
                              pch_within = c(21, 22),
                              pch_between = 24,
                              alpha = 0.4,
                              cex = 1, plot_legend=FALSE, levels = c(0, 0.77), ...) {
  # Get item pairs (upper triangle only to avoid duplicates)
  p <- nrow(Delta)
  pairs <- which(upper.tri(Delta), arr.ind = TRUE)

  # Extract pairwise values and pair types
  n_pairs <- nrow(pairs)
  delta_vals <- numeric(n_pairs)
  logOR_vals <- numeric(n_pairs)
  pair_type <- character(n_pairs)

  for (k in 1:n_pairs) {
    i <- pairs[k, 1]
    j <- pairs[k, 2]
    delta_vals[k] <- Delta[i, j]
    logOR_vals[k] <- logOR[i, j]

    # Determine pair type
    if (scales[i] == scales[j]) {
      pair_type[k] <- paste0("within", scales[i])
    } else {
      pair_type[k] <- "between"
    }
  }

  # Assign border colors, fill colors (with transparency), and shapes
  border_cols <- ifelse(pair_type == "within1", col_within[1],
                 ifelse(pair_type == "within2", col_within[2], col_between))
  fill_cols <- ifelse(pair_type == "within1", adjustcolor(col_within[1], alpha.f = alpha),
               ifelse(pair_type == "within2", adjustcolor(col_within[2], alpha.f = alpha),
                      adjustcolor(col_between, alpha.f = alpha)))
  pchs <- ifelse(pair_type == "within1", pch_within[1],
          ifelse(pair_type == "within2", pch_within[2], pch_between))

  # Plot Delta vs log-OR
  plot(delta_vals, logOR_vals, type = "n",
       xlab = expression(Delta),
       ylab = "log-OR",
       main = main, cex.main=1.3, cex.lab=1.2, ...)

  for(i in 1:length(levels)) {
    abline(h=levels[i], lty=2, col="darkgrey")
  }
  points(delta_vals, logOR_vals, pch = pchs, cex = cex, col = border_cols, bg = fill_cols)

  # Add legend
  if(plot_legend) {
     legend("topleft", inset=0.02, bty="n",
         legend = c("Within scale 1", "Within scale 2", "Between scales"),
         col = c(col_within, col_between),
         pt.bg = c(adjustcolor(col_within[1], alpha.f = alpha),
                   adjustcolor(col_within[2], alpha.f = alpha),
                   adjustcolor(col_between, alpha.f = alpha)),
         pch = c(pch_within, pch_between), cex = 1.2)
  }

  # Return data invisibly for further analysis
  invisible(data.frame(i = pairs[,1], j = pairs[,2],
                       Delta = delta_vals, logOR = logOR_vals,
                       type = pair_type, col = border_cols))
}

source("iota17.R")

x11(height=10, width=10)
par(mfrow=c(2,2))

p1 = 8; p2 = 8
rho_ex = 0.5

# rho=0, homogeneous
D1 = DiG$new(NULL)
D1$generate_Rasch_scales(deltas = list(seq(-2,2,len=p1), seq(-2,2,len=p2)),
                         sigmas = c(1.0, 1.0),
                         R = matrix(c(1,0,0,1), 2, 2),
                         Nobs = 2000,
                         seed = 123456)
D1$compute()
plot_Delta_logOR(D1$Delta, D1$logOR, D1$true.params$scale_membership, main=expression(paste("(a) ",sigma[1]==1, ", ", sigma[2]==1, ", ", rho==0)), ylim=c(-.5,2.5), plot_legend=TRUE, levels=c(pred_logOR(1,1,0),pred_logOR(1,1,1)))


# rho=0.5, homogeneous
D2 = DiG$new(NULL)
D2$generate_Rasch_scales(deltas = list(seq(-2,2,len=p1), seq(-2,2,len=p2)),
                         sigmas = c(1.0, 1.0),
                         R = matrix(c(1,rho_ex,rho_ex,1), 2, 2),
                         Nobs = 2000,
                         seed = 123456)
D2$compute()
plot_Delta_logOR(D2$Delta, D2$logOR, D2$true.params$scale_membership, main=expression(paste("(b) ",sigma[1]==1, ", ", sigma[2]==1, ", ", rho==.5)), ylim=c(-.5,2.5), levels=c(pred_logOR(1,1,1),pred_logOR(1,1,rho_ex)))

# rho=0, heterogeneous
D3 = DiG$new(NULL)
D3$generate_Rasch_scales(deltas = list(seq(-2,2,len=p1), seq(-2,2,len=p2)),
                         sigmas = c(2.0, 1.0),
                         R = matrix(c(1,0,0,1), 2, 2),
                         Nobs = 2000,
                         seed = 123456)
D3$compute()
plot_Delta_logOR(D3$Delta, D3$logOR, D3$true.params$scale_membership, main=expression(paste("(c) ",sigma[1]==2, ", ", sigma[2]==1, ", ", rho==0)), ylim=c(-.5,2.5), levels=c(pred_logOR(1,1,1),pred_logOR(2,2,1),pred_logOR(1,1,0)))

# rho=0.5, heterogeneous
D4 = DiG$new(NULL)
D4$generate_Rasch_scales(deltas = list(seq(-2,2,len=p1), seq(-2,2,len=p2)),
                         sigmas = c(2.0, 1.0),
                         R = matrix(c(1,rho_ex,rho_ex,1), 2, 2),
                         Nobs = 2000,
                         seed = 123456)
D4$compute()
plot_Delta_logOR(D4$Delta, D4$logOR, D4$true.params$scale_membership, main=expression(paste("(d) ",sigma[1]==2, ", ", sigma[2]==1, ", ", rho==.5)), ylim=c(-.5,2.5), levels=c(pred_logOR(1,1,1), pred_logOR(2,2,1), pred_logOR(2,1,rho_ex)))

dev.copy2pdf(file="iota-logOR.pdf")

#---------------------------- Figure 4: item coordinates for the 4 simulated datasets ----------------------------------
source("iota17.R")

# Helper function: draw significant implication arrows on item coordinate plots
add_sig_arrows = function(D, scale_membership, curve=0.25, lwd=0.9, arr.length=0.12, alpha_val=0.7,
                          scale2_items=9:16) {

  require(diagram)

  coords = D$coordinates
  p = nrow(coords)
  # scale_colors = c("steelblue", "firebrick")
  scale_colors = google.colors[1:2]

  # Helper: choose curve direction — scale-2 within-links bend left
  link_curve = function(i, j) {
    if (scale_membership[i] == 2 && scale_membership[j] == 2) -curve else curve
  }

  # Directed arrows (G.uni): single-headed
  for (i in 1:p) {
    for (j in 1:p) {
      if (i != j && D$G.uni[i,j] == 1) {
        s_i = scale_membership[i]
        s_j = scale_membership[j]
        between = (s_i != s_j)
        col = if (!between) adjustcolor(scale_colors[s_i], alpha.f=alpha_val)
              else "darkgrey"
        if (!between) {
          curvedarrow(from=coords[i,], to=coords[j,], curve=link_curve(i,j), lwd=lwd,
                      lcol=col, arr.col=col, arr.pos=0.7, arr.length=arr.length,
                      arr.type="triangle")
        } else {
          straightarrow(from=coords[i,], to=coords[j,], lwd=lwd,
                        lcol=col, arr.col=col, arr.pos=0.7, arr.length=arr.length,
                        arr.type="triangle")
        }
      }
    }
  }

  # Symmetric links (G.bi): upper triangle only
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      if (D$G.bi[i,j] == 1) {
        s_i = scale_membership[i]
        s_j = scale_membership[j]
        between = (s_i != s_j)
        if (!between) {
          # Within-scale: double-headed arrow, coloured by scale
          col = adjustcolor(scale_colors[s_i], alpha.f=alpha_val)
          curvedarrow(from=coords[i,], to=coords[j,], curve=link_curve(i,j), lwd=lwd,
                      lcol=col, arr.col=col, arr.pos=0.7, arr.length=arr.length,
                      arr.type="triangle", endhead=TRUE)
        } else {
          # Between-scale: straight grey segment, no arrowheads
          col = adjustcolor("black", alpha.f=alpha_val)
          segments(x0=coords[i,1], y0=coords[i,2],
                   x1=coords[j,1], y1=coords[j,2],
                   lwd=lwd, col=col)
        }
      }
    }
  }
}

x11(width=10,height=10)
par(mfrow=c(2,2))

# (a) sigma1=1, sigma2=1, rho=0
D1$cluster_logOR_levels()
lev_a = D1$logOR_levels$centers
plot(D1$coordinates[,1], D1$coordinates[,2], type="n",main=expression(paste("(a) ",sigma[1]==1, ", ", sigma[2]==1, ", ", rho==0)), xlim=c(-2.2,2.2), ylim=c(-2.1,2.1),xlab="Difficulty", ylab="Associative locations", cex.main=1.3, cex.lab=1.2)
abline(h=lev_a, lwd=2, lty=2, col="lightgray")
add_sig_arrows(D1, D1$true.params$scale_membership)
points(D1$coordinates[,1], D1$coordinates[,2], pch=21, col=gray(.3), bg="white")
text(D1$coordinates[,1], D1$coordinates[,2]+0.2, paste("i",1:16,sep=""))

# (b) sigma1=1, sigma2=1, rho=0.5
D2$cluster_logOR_levels()
lev_b = D2$logOR_levels$centers
plot(D2$coordinates[,1], D2$coordinates[,2], type="n",main=expression(paste("(b) ",sigma[1]==1, ", ", sigma[2]==1, ", ", rho==.5)), xlim=c(-2.2,2.2),ylim=c(-2.1,2.1),xlab="Difficulty", ylab="Associative locations", cex.main=1.3, cex.lab=1.2)
abline(h=lev_b, lwd=2, lty=2, col="lightgray")
add_sig_arrows(D2, D2$true.params$scale_membership)
points(D2$coordinates[,1], D2$coordinates[,2], pch=21, col=gray(.3), bg="white")
text(D2$coordinates[,1], D2$coordinates[,2]+0.2, paste("i",1:16,sep=""))
legend("topright", inset=0.02,
       legend = c("Within scale A", "Within scale B", "Between (symmetric)", "Between (directed)"),
       col = c(adjustcolor(google.colors[1], alpha.f=0.7), adjustcolor(google.colors[2], alpha.f=0.7),
               "black", adjustcolor("grey50", alpha.f=0.7)),
       lty = 1, lwd = 1.8, cex = 0.85, bty = "n")

# (c) sigma1=2, sigma2=1, rho=0
D3$cluster_logOR_levels()
lev_c = D3$logOR_levels$centers
plot(D3$coordinates[,1], D3$coordinates[,2], type="n",main=expression(paste("(c) ",sigma[1]==2, ", ", sigma[2]==1, ", ", rho==0)), xlim=c(-2.2,2.2), ylim=c(-2.1,2.1),xlab="Difficulty", ylab="Associative locations", cex.main=1.3, cex.lab=1.2)
abline(h=lev_c, lwd=2, lty=2, col="lightgray")
add_sig_arrows(D3, D3$true.params$scale_membership)
points(D3$coordinates[,1], D3$coordinates[,2], pch=21, col=gray(.3), bg="white")
text(D3$coordinates[,1], D3$coordinates[,2]+0.2, paste("i",1:16,sep=""))

# (d) sigma1=2, sigma2=1, rho=0.5
D4$cluster_logOR_levels()
lev_d = D4$logOR_levels$centers
plot(D4$coordinates[,1], D4$coordinates[,2], type="n",main=expression(paste("(d) ",sigma[1]==2, ", ", sigma[2]==1, ", ", rho==.5)), xlim=c(-2.2,2.2), ylim=c(-2.1,2.1),xlab="Difficulty", ylab="Associative locations", cex.main=1.3, cex.lab=1.2)
abline(h=lev_d, lwd=2, lty=2, col="lightgray")
add_sig_arrows(D4, D4$true.params$scale_membership)
points(D4$coordinates[,1], D4$coordinates[,2], pch=21, col=gray(.3), bg="white")
text(D4$coordinates[,1], D4$coordinates[,2]+0.2, paste("i",1:16,sep=""))

dev.copy2pdf(file="item-coordinates.pdf")


#---------------------------- Figure 6: Decomposition plot for the real data example (TDRI test) ----------------------------------
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

x11(height=8,width=8)

source("iota17.R")

# Construct the main analysis object
Dt = DiG$new(tdri)

# 7 communities detected
Dt$compute(rotate=FALSE)

# 5 logOR levels detected
Dt$cluster_logOR_levels()
# Dt$logOR_levels

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

# Plot decomposition
Dt$plot.decp(color.by="community", show.labels="clusters", xlim=c(-2,8), ylim=c(-2,8), community.labels=stages, show.pair.subsets = TRUE)
dev.copy2pdf(file="TDRI-decomposition.pdf")

#---------------------------- Figure 7: Final graph for the real data example (TDRI test) ----------------------------------

# Plot implication graph and metagraph
Dt$plot.tikz(file="TDRI-full-meta-graph-iota16-unrotated.tex", plot_communities=TRUE, community_type="ellipses", plot_metagraph=TRUE, meta.bend_right=1:7, meta.bend_angle=70, only_directed=FALSE, only_symmetric=FALSE, bend_angle=30, dilate_cluster=0.8, scale_factor=5, cluster_labels=stages, cluster_labels_below=TRUE, plot_logOR_levels=TRUE, plot_difficulty_gradient=TRUE)

# Optional: Recompute for best alignment with difficulty
Dt$compute(rotate=TRUE)
remap_clusters_to_stages(Dt)
Dt$plot.tikz(file="TDRI-full-meta-graph-iota16-rotated.tex", plot_communities=TRUE, community_type="ellipses", plot_metagraph=TRUE, meta.bend_right=1:7, meta.bend_angle=60, only_directed=FALSE, only_symmetric=FALSE, bend_angle=30, dilate_cluster=0.8, scale_factor=5, cluster_labels=stages) 
