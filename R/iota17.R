library(R6)
library(igraph)
library(diagram)
library(igraph)
library(mvtnorm)
library(mclust,quietly=TRUE)
library(blockmodels)                   # fast variational SBM
library(MASS)                          # for isoMDS (nonmetric MDS)
library(mirt)

# Dimensional iota Graph class
DiG = R6Class("DiG",

  public = list(

    mat       = NULL,
    names     = NULL,
    zero.val  = NULL,
    alpha     = 0.05,
    FDR       = TRUE,
    Delta     = NULL,
    SE_asym   = NULL,
    Z_asym    = NULL,
    P_asym    = NULL,
    logOR     = NULL,
    SE_OR     = NULL,
    Z_OR      = NULL,
    P_OR      = NULL,
    IOTA_STAR = NULL,
    P_IOTA    = NULL,
    P_IOTA_adj = NULL,
    dist_Delta = NULL,
    dist_logOR = NULL,
    postProb  = NULL,
    G         = NULL,
    G.uni     = NULL,
    G.bi      = NULL,
    iGraph    = NULL,
    Intensity = NULL,
    success.rates = NULL,
    full.distances = NULL,
    coordinates = NULL,
    cor2success = NULL,
    rotate = NULL,
    eig = NULL,
    ndim = NULL,
    clustfit = NULL,
    scales = NULL,
    gaussian.ellipses = NULL,
    communities = NULL,
    modularity = NULL,
    metagraph = list(),
    true.params = NULL,
    logOR_levels = NULL,

    # Constructor
    initialize = function(mat, names=NULL) {
      self$mat = mat
      if(is.null(names)) self$names = colnames(mat)
      else self$names = names
    },

    # Generate data from a multidimensional Rasch model with arbitrary scales and arbitrary correlations
    #
    # Arguments:
    #   deltas  - list of numeric vectors, one per scale (item difficulties)
    #   sigmas  - numeric vector of latent standard deviations (length = number of scales)
    #   R       - correlation matrix between latent traits (default: identity = independent scales)
    #   Nobs    - number of observations
    #
    # Example:
    #   D = DiG$new(NULL)
    #   D$generate_Rasch_scales(deltas = list(c(-2,-1,0,1,2), c(-1,0,1)),
    #                           sigmas = c(1.5, 1.0),
    #                           R = matrix(c(1,0.3,0.3,1), 2, 2))
    #   D$compute()

    generate_Rasch_scales = function(deltas, sigmas=1, R=NULL, Nobs=2000, seed=12345678) {

      set.seed(seed)

      n_scales = length(deltas)
      n_items  = sapply(deltas, length)
      p_total  = sum(n_items)

      # Recycle sigmas if scalar
      if(length(sigmas) == 1) sigmas = rep(sigmas, n_scales)

      # Default: independent scales
      if (is.null(R)) R = diag(n_scales)

      # Build covariance matrix from sigmas and correlation matrix
      S = diag(sigmas) %*% R %*% diag(sigmas)

      # Build discrimination matrix (p_total x n_scales): each item loads on its own scale only
      a = matrix(0, nrow = p_total, ncol = n_scales)
      d = numeric(p_total)
      row = 0
      for (k in 1:n_scales) {
        for (i in 1:n_items[k]) {
          row = row + 1
          a[row, k] = 1
          d[row] = -deltas[[k]][i]
        }
      }

      items = rep('2PL', p_total)
      X = simdata(a, d, Nobs, itemtype = items, sigma = S)

      # Store in mat field
      self$mat = X
      self$names = colnames(X)

      # Store true parameters
      self$true.params = list(
        deltas = deltas,
        sigmas = sigmas,
        R = R,
        Nobs = Nobs,
        n_scales = n_scales,
        n_items = n_items,
        scale_membership = rep(1:n_scales, n_items)
      )
    },

    # Compute the implication matrix
    compute = function(max_iter=5000, zero.val=1, alpha=0.05, FDR=TRUE, rotate=TRUE, inference="Wald", credibility=0.6, cluster.method="SBM", nonmetric=FALSE) {

      self$zero.val = zero.val
      self$alpha    = alpha
      self$FDR      = FDR
      self$rotate   = rotate

      # Case 1: Data is a matrix or data.frame of more than 2 binary variables
      if( any(class(self$mat) %in% c("data.frame","matrix")) && ncol(self$mat) > 2) {

        X1 = as.matrix(self$mat)
        self$success.rates = colMeans(X1,na.rm=TRUE)

        # Discard empty rows
        sel.rows = rowSums(is.na(X1)) != ncol(X1)
        X1 = X1[sel.rows,]

        X0 = 1-X1
        
        # Deal with missing values
        notmissing = crossprod(!is.na(X1),!is.na(X1))
        X1[is.na(X1)] = 0
        X0[is.na(X0)] = 0

        # Co-occurrences
        n11 = crossprod(X1,X1)
        n10 = crossprod(X1,X0)
        n01 = crossprod(X0,X1)
        n00 = crossprod(X0,X0)

        # Bayesian approach
        if(inference == "bayesian") {

          # [TODO: Manage missing values]
          # Compute posterior probabilities
          self$Intensity = self$postProb.matrix(n11,n10,n00)
          diag(self$Intensity) = 0
        
          # Construct graph
          self$G = (self$Intensity > credibility) + 0
        }

        else if(inference == "Wald") {

          # Pseudo-bayesian correction in case of any zero cell
          is.zero = (n11 == 0) | (n10 == 0) | (n01 == 0) | (n00 == 0)
          correction = is.zero * zero.val
          if(any(is.zero)){
            n11 = n11 + correction
            n10 = n10 + correction
            n01 = n01 + correction
            n00 = n00 + correction
          }

          # Compute joint frequencies (potentially adjusted for different totals)
          # f11 = n11 / (notmissing + 4*correction)
          # f10 = n10 / (notmissing + 4*correction)
          # f01 = n01 / (notmissing + 4*correction)
          # f00 = n00 / (notmissing + 4*correction)

          ## ODD RATIO TEST

          # Compute one-sided test of positive log-odds ratio
          # [TODO]: Detect negative log-odds ratios
          self$logOR = log(n11) + log(n00) - log(n10) - log(n01)
          diag(self$logOR) = NA
          self$SE_OR = sqrt((1/n11) + (1/n00) + (1/n01) + (1/n10))
          self$Z_OR = self$logOR / self$SE_OR
          self$P_OR = 1-pnorm(self$Z_OR)

          ## ASYMMETRY TEST

          # Compute two-sided test of asymmetry
          self$Delta = log(n01) - log(n10)
          self$SE_asym = sqrt((1/n01) + (1/n10))
          self$Z_asym = self$Delta / self$SE_asym
          self$P_asym = 2 * (1-pnorm(abs(self$Z_asym)))  # Two-sided p-value

          ## Union-intersection test on IOTA_STAR
          self$IOTA_STAR = self$logOR + self$Delta
          diag(self$IOTA_STAR) = 0

          # Exclude diagonal from p-value matrices (self-comparisons are not valid tests)
          diag(self$P_OR) = NA
          diag(self$P_asym) = NA

          self$P_IOTA = pmax(self$P_OR, self$P_asym)

          # Make P_IOTA directional: only report p-value where Delta > 0
          # (i.e., where row item implies column item)
          # Positions where Delta <= 0 are set to NA (not the inferred direction)
          self$P_IOTA[self$Delta <= 0] = NA

          # Implication intensity (complement of combined p-value)
          self$Intensity = 1 - self$P_IOTA

          # Apply Benjamini-Yekutieli correction for FDR control under arbitrary dependence
          # Hierarchical/gatekeeping approach:
          # 1. First gate: Apply BY to P_OR (lower triangle only, since logOR is symmetric)
          # 2. Second gate: Apply BY to P_asym ONLY for pairs with significant OR
          # This preserves power by not penalizing for tests that won't be used
          if(FDR) {
            # Step 1: Apply BY correction to P_OR on lower triangle only (logOR is symmetric)
            lower_mask = lower.tri(self$P_OR)
            p_OR_lower = self$P_OR[lower_mask]
            p_OR_lower_adj = p.adjust(p_OR_lower, method = "BY")

            # Create symmetric P_OR_adj matrix (mirror lower to upper)
            P_OR_adj_matrix = matrix(NA, nrow = nrow(self$P_OR), ncol = ncol(self$P_OR))
            dimnames(P_OR_adj_matrix) = dimnames(self$P_OR)
            P_OR_adj_matrix[lower_mask] = p_OR_lower_adj
            P_OR_adj_matrix[upper.tri(P_OR_adj_matrix)] = t(P_OR_adj_matrix)[upper.tri(P_OR_adj_matrix)]

            # Step 2: Identify pairs with significant OR in lower triangle (first gate)
            sig_OR_lower = !is.na(p_OR_lower_adj) & (p_OR_lower_adj < alpha)

            # Step 3: Apply BY correction to P_asym ONLY for significant OR pairs in lower triangle
            # Initialize P_asym_adj with NA (will be filled directionally)
            P_asym_adj_matrix = matrix(NA, nrow = nrow(self$P_asym), ncol = ncol(self$P_asym))
            dimnames(P_asym_adj_matrix) = dimnames(self$P_asym)

            if(sum(sig_OR_lower) > 0) {
              # Get p-values for significant OR pairs
              p_asym_subset = self$P_asym[lower_mask][sig_OR_lower]
              p_asym_subset_adj = p.adjust(p_asym_subset, method = "BY")

              # Get row/column indices for significant pairs in lower triangle
              lower_indices = which(lower_mask, arr.ind = TRUE)  # Returns matrix with row, col
              sig_indices = lower_indices[sig_OR_lower, , drop = FALSE]

              # Store directionally based on sign of Delta
              # Convention: row i implies column j when stored at [i,j]
              for(idx in seq_len(nrow(sig_indices))) {
                i = sig_indices[idx, 1]
                j = sig_indices[idx, 2]
                if(self$Delta[i, j] > 0) {
                  # i → j: store at [i,j]
                  P_asym_adj_matrix[i, j] = p_asym_subset_adj[idx]
                } else if(self$Delta[i, j] < 0) {
                  # j → i: store at [j,i]
                  P_asym_adj_matrix[j, i] = p_asym_subset_adj[idx]
                }
                # If Delta == 0, no direction, leave as NA
              }
            }

            # For graph construction, need symmetric version of P_asym_adj
            # Initialize with 1 (non-significant), fill with adjusted p-values, mirror to upper
            P_asym_adj_sym = matrix(1, nrow = nrow(self$P_asym), ncol = ncol(self$P_asym))
            dimnames(P_asym_adj_sym) = dimnames(self$P_asym)
            if(sum(sig_OR_lower) > 0) {
              P_asym_adj_sym[lower_mask][sig_OR_lower] = p_asym_subset_adj
            }
            P_asym_adj_sym[upper.tri(P_asym_adj_sym)] = t(P_asym_adj_sym)[upper.tri(P_asym_adj_sym)]
            diag(P_asym_adj_sym) = NA

            P_IOTA_adj_sym = pmax(P_OR_adj_matrix, P_asym_adj_sym, na.rm = TRUE)

            # P_IOTA_adj is directional: only where Delta > 0
            self$P_IOTA_adj = pmax(P_OR_adj_matrix, P_asym_adj_matrix, na.rm = TRUE)
            self$P_IOTA_adj[is.na(P_asym_adj_matrix)] = NA  # Keep directional structure

          } else {
            # No FDR correction (not recommended)
            P_OR_adj_matrix = self$P_OR
            P_IOTA_adj_sym = pmax(self$P_OR, self$P_asym)  # Keep symmetric for graph construction
            self$P_IOTA_adj = self$P_IOTA
          }

          # Two-step graph construction:
          # 1. G.bi (symmetric): significant association (OR) but NOT significant directed dependence
          #    → bidirectional links (association without clear directional preference)
          # 2. G.uni (directed): significant directed dependence (both OR and asymmetry)
          #    → unidirectional links (direction given by sign of Delta)
          # Note: use P_IOTA_adj_sym (symmetric) for graph construction

          sig_OR       = (P_OR_adj_matrix < alpha)
          sig_directed = (P_IOTA_adj_sym < alpha)

          # G.bi: symmetric association (significant OR but NOT significant directed dependence)
          # These are pairs with association but no clear directional preference
          self$G.bi = (sig_OR & !sig_directed) + 0
          self$G.bi[is.na(self$P_OR)] = 0
          diag(self$G.bi) = 0

          # G.uni: directed implication (significant directed dependence, direction by Delta sign)
          # For position [i,j]: Delta > 0 means A -> B (row i implies column j)
          # Note: sig_directed implies sig_OR (since P_IOTA = max(P_OR, P_asym))
          self$G.uni = (sig_directed & (self$Delta > 0)) + 0
          self$G.uni[is.na(self$P_IOTA)] = 0
          diag(self$G.uni) = 0

          # G: overall graph (union of symmetric and directed links)
          self$G = ((self$G.bi > 0) | (self$G.uni > 0)) + 0
        }

        else {
          stop("Please choose a valid inference.type argument: 'Wald' or 'bayesian'")
        }

        #-- Compute coordinates via two separate 1D MDS

        # 1. Dissimilarity from |Delta| (symmetric: |Delta[i,j]| = |Delta[j,i]|)
        self$dist_Delta = abs(self$Delta)
        diag(self$dist_Delta) = 0

        # 2. Dissimilarity from logOR: max(logOR) - logOR
        #    Higher logOR (stronger association) -> smaller distance
        max_logOR = max(self$logOR, na.rm = TRUE)
        self$dist_logOR = max_logOR - self$logOR
        diag(self$dist_logOR) = 0

        # Detect and store the list of variables with no missing inter-distances
        # (iteratively remove variables with large number of missing values)
        self$full.distances = 1:ncol(X1)
        D_check = self$dist_Delta
        while(any( (S = colSums(!is.finite(D_check))) > 0 )) {
          sel = which.max(S)
          D_check = D_check[-sel,-sel]
          self$full.distances = self$full.distances[-sel]
        }

        # Suppress warnings about negative eigenvalues (expected for non-Euclidean distances)
        options(warn = -1)

        # 1D MDS on |Delta| dissimilarities (difficulty axis)
        if(sum(!is.finite(self$dist_Delta)) == 0) {
          mds_delta = cmdscale(self$dist_Delta, k = 1, eig = TRUE)
          if(nonmetric) {
            mds_nm = isoMDS(self$dist_Delta, y = mds_delta$points, k = 1, trace = FALSE)
            coord_delta = mds_nm$points
          } else {
            coord_delta = mds_delta$points
          }
        } else {
          if(nonmetric) warning("Nonmetric MDS with missing distances not supported. Using metric MDS with EM imputation.")
          mds_delta = self$iterative_mds_em(self$dist_Delta, max_iter = max_iter, k = 1, eig = TRUE)
          coord_delta = mds_delta$points
        }

        # 1D MDS on max(logOR) - logOR dissimilarities (association axis)
        if(sum(!is.finite(self$dist_logOR)) == 0) {
          mds_logOR = cmdscale(self$dist_logOR, k = 1, eig = TRUE)
          if(nonmetric) {
            mds_nm = isoMDS(self$dist_logOR, y = mds_logOR$points, k = 1, trace = FALSE)
            coord_logOR = mds_nm$points
          } else {
            coord_logOR = mds_logOR$points
          }
        } else {
          if(nonmetric) warning("Nonmetric MDS with missing distances not supported. Using metric MDS with EM imputation.")
          mds_logOR = self$iterative_mds_em(self$dist_logOR, max_iter = max_iter, k = 1, eig = TRUE)
          coord_logOR = mds_logOR$points
        }

        # Store eigenvalues from both MDS (for checking 1D adequacy)
        self$eig = list(Delta = mds_delta$eig, logOR = mds_logOR$eig)

        # Report 1D MDS fit
        if(!is.null(self$eig$Delta)) {
          pos_eig = sum(self$eig$Delta[self$eig$Delta > 0])
          neg_eig = sum(abs(self$eig$Delta[self$eig$Delta < 0]))
          cat("MDS on |Delta|: 1st eigenvalue explains", round(100 * self$eig$Delta[1] / pos_eig, 1), "%",
              " | non-Euclidean ratio:", round(neg_eig / pos_eig, 3), "\n")
        }
        if(!is.null(self$eig$logOR)) {
          pos_eig = sum(self$eig$logOR[self$eig$logOR > 0])
          neg_eig = sum(abs(self$eig$logOR[self$eig$logOR < 0]))
          cat("MDS on logOR:   1st eigenvalue explains", round(100 * self$eig$logOR[1] / pos_eig, 1), "%",
              " | non-Euclidean ratio:", round(neg_eig / pos_eig, 3), "\n")
        }

        options(warn = 0)

        # Orient axes before rotation
        # i) Difficult items on the right (negative correlation with success rates)
        if(cor(as.vector(coord_delta), self$success.rates) > 0) {
          coord_delta = -coord_delta
        }
        # ii) Higher logOR at the top (positive correlation with mean logOR per item)
        mean_logOR = rowMeans(self$logOR, na.rm = TRUE)
        if(cor(as.vector(coord_logOR), mean_logOR) < 0) {
          coord_logOR = -coord_logOR
        }

        # Combine into 2D coordinates
        self$ndim = 2
        self$coordinates = cbind(as.vector(coord_delta), as.vector(coord_logOR))
        rownames(self$coordinates) = colnames(X1)
        colnames(self$coordinates) = c("Delta", "logOR")

        # Clustering (warning: Coordinates must have been computed first)
        self$detect_graphCommunities(method = cluster.method)

        # Rotate the coordinates to maximize the correlation between the first dimension and the success rate
        if(rotate) {

          C = self$coordinates
          y = self$success.rates

          # Step 1: Find the direction maximally correlated with y
          w1 = t(C) %*% y
          w1 = w1 / norm(w1, "2")

          # Step 2: Find the orthogonal direction (first eigenvector of the residual)
          C_residual = C - C %*% w1 %*% t(w1)
          w2 = eigen(t(C_residual) %*% C_residual)$vectors[, 1, drop = FALSE]

          # Step 3: Create rotation matrix and transform
          R = cbind(w1, w2)            # 2x2 rotation matrix
          self$coordinates = C %*% R   # rotated coordinates

          # Reverse (IRT convention: Easier items on the left)
          self$coordinates[,1] = -self$coordinates[,1]

          # Correlation between the first component and the success rate
          # (an inverse indicator of a nominal structure)
          self$cor2success = cor(self$coordinates[,1], y)
          cat("Correlation between first dimension and success rate / prevalence:\nR =", self$cor2success, "\n")
        }

      }

      # Case 2: mat is a Simple 2x2 count table
      else {

        # Data is a 2x2 count table
        if(any(class(self$mat) %in% c("table","xtabs"))) {

          if((ncol(self$mat) != 2) || (nrow(self$mat) != 2)) {
            cat("Table should have only two rows and two columns.\n")
            return()
          }

          if(!all(colnames(self$mat) %in% c("0","1"))) { 
            cat("Table should have '0' and '1' row names\n")
            return()
          }
          if(!all(rownames(self$mat) %in% c("0","1"))) { 
            cat("Table should have '0' and '1' row names\n")
            return()
          }

          # Table objects have no dimnames names set
          if(all(names(dimnames(self$mat)) == "")) {
            attr(self$mat,"dimnames") = list("Var1" = c(0,1), "Var2" = c(0,1))
          }
          tab = self$mat
          varlabels = names(dimnames(tab))
        }

        # Data is a 2x2 matrix
        else if(any(class(self$mat) == "matrix") && ncol(self$mat) == 2 && nrow(self$mat) == 2) {

          if(!all(colnames(self$mat) %in% c("0","1"))) { 
            cat("Table should have '0' and '1' row names\n")
            return()
          }
          if(!all(rownames(self$mat) %in% c("0","1"))) { 
            cat("Table should have '0' and '1' row names\n")
            return()
          }
          tab = self$mat
          varlabels = c("Var1","Var2")
        }

        # Case 3: mat is a 2-column matrix or data.frame
        if( (nrow(self$mat) > 2) && ((any(class(self$mat) == "matrix") && ncol(self$mat) == 2) || (any(class(self$mat) == "data.frame") && ncol(self$mat) == 2)) ) {

          tab = table(self$mat[,1],self$mat[,2])
          varlabels = colnames(self$mat)

          if(all(colnames(tab) != c("0","1"))) { 
            cat("Data should contains only '0' and '1' values.\n")
            return()
          }
          if(all(rownames(tab) != c("0","1"))) { 
            cat("Data should contains only '0' and '1' values.\n")
            return()
          }
        }

        cat("Implication table\n")
        cat(paste(varlabels,collapse=" -> "),"\n")
        print(tab[2:1,2:1])
        cat("\n")

        n11 = tab[2,2]
        n10 = tab[2,1]
        n00 = tab[1,1]

        # Deal with zero cells
        not.defined = (n11==0) || (n10==0) || (n00==0)
        
        if(not.defined) {

          if(is.null(zero.val)) {
            cat("Zero counts encountered. Log-counts are not computable. Please consider using the zero.val argument.\n")
            return()
          }

          # Pseudo-Bayesian correction
          else {
            n11 = n11 + zero.val
            n10 = n10 + zero.val
            n00 = n00 + zero.val
          }      
        }

        # Compute iota
        iota = log(n11)+log(n00)-2*log(n10)
        SE = sqrt((1/n00)+(4/n10)+(1/n11))

        z = iota/SE
        p = pnorm(z)

        data.frame(iota_star=iota,SE=SE,z=z,p.value=1-p,Intensity=p)
      }
    },

    # Compute posterior probabilities
    logspace_add = function(a, b) {

      ## numerically stable log(exp(a)+exp(b))
      if (is.infinite(a)) return(b)
      if (is.infinite(b)) return(a)
      if (a < b) { tmp = a; a = b; b = tmp }
      a + log1p(exp(b - a))
    },

    postProb0 = function(n11, n10, n00) {

      A = n11 + 1L          # shape parameters of the posterior
      B = n10 + 1L
      C = n00 + 1L

      logP = -Inf          # log-sum-exp accumulator
      log3 = log(3)

      for (i in 0:(A - 1L)) {
        log_i_fact = lgamma(i + 1)         # log(i!)
        for (j in 0:(C - 1L)) {
          log_term = lgamma(B + i + j) - lgamma(B) - log_i_fact - lgamma(j + 1) - (B + i + j) * log3
          logP = self$logspace_add(logP, log_term)
        }
      }
      exp(logP)
    },

    ## ------------------------------------------------------------------
    ##                 Compute posterior probabilities
    ## ------------------------------------------------------------------
    logSumExp = function(v) {
      m = max(v)
      m + log(sum(exp(v - m)))
    },

    postProb.scalar = function(n11, n10, n00) {

      A = n11 + 1L               # Dirichlet shapes
      B = n10 + 1L
      C = n00 + 1L
      log3 = log(3)

      ## vector of i = 0:(A-1) and j = 0:(C-1)
      i = 0:(A - 1L)
      j = 0:(C - 1L)

      ## build the (A × C) matrix of log-terms via outer()
      log_terms = outer(
          i, j,
          FUN = function(i, j)
            lgamma(B + i + j) - lgamma(B) -
            lgamma(i + 1)    - lgamma(j + 1) -
            (B + i + j) * log3
      )

      ## log-sum-exp over all entries then exponentiate
      exp( self$logSumExp(as.vector(log_terms)) )
    },

    postProb.matrix = function(n11, n10, n00) {

      if (!all(dim(n11) == dim(n10)) || !all(dim(n11) == dim(n00)))
        stop("Compute posterior probabilities: Input matrices must have identical dimensions.\n")

      ## convert the matrices to vectors, apply the scalar engine,
      ## then reshape back to a matrix
      cat("Computing exact posterior probabilities. Please be patient...\n")
      probs = mapply(self$postProb.scalar, n11 = as.vector(n11), n10 = as.vector(n10), n00 = as.vector(n00))

      matrix(probs, nrow = nrow(n11), ncol = ncol(n11), dimnames = dimnames(n11))
    },

    # Search for scales by model-based clustering
    # Options:
    # "EII" spherical, equal volume
    # "VII" spherical, unequal volume
    # "EEI" diagonal, equal volume and shape
    # "VEI" diagonal, varying volume, equal shape
    # "EVI" diagonal, equal volume, varying shape
    # "VVI" diagonal, varying volume and shape
    # "EEE" ellipsoidal, equal volume, shape, and orientation
    # "VEE" ellipsoidal, equal shape and orientation (*)
    # "EVE" ellipsoidal, equal volume and orientation (*)
    # "VVE" ellipsoidal, equal orientation (*)
    # "EEV" ellipsoidal, equal volume and equal shape
    # "VEV" ellipsoidal, equal shape
    # "EVV" ellipsoidal, equal volume (*)
    # "VVV" ellipsoidal, varying volume, shape, and orientatio
    cluster = function(model="VVE") {
    
      if(self$ndim == 1) {
        self$scales = rep(1,self$p)
        return()
      }

      # Subscale detection by model based clustering
      self$clustfit = Mclust(self$coordinates[,1:self$ndim],modelNames=model)
      self$scales = self$clustfit$classification
      names(self$scales) = rownames(self$coordinates)

      self$gaussian.ellipses = self$extract_gaussianEllipses(confidence = 0.95)
      cat(length(unique(self$scales)),"scales detected.\n")
    },

    # Get gaussian component parameters
    extract_gaussianEllipses = function(confidence = 0.95) {
      
      # Centers
      means = self$clustfit$parameters$mean
      
      # Variance and covariance matrices
      if (self$clustfit$modelName %in% c("EII", "VII", "EEI", "VEI", "EVI", "VVI")) {
        # Modèles avec covariances diagonales/sphériques
        sigmas = self$clustfit$parameters$variance$sigmasq
        if (length(dim(sigmas)) == 0) sigmas = array(sigmas, c(1,1,ncol(means)))
      } else {
        # Modèles avec covariances complètes
        sigmas = self$clustfit$parameters$variance$sigma  # (2 x 2 x G)
      }
      
      n_components = ncol(means)
      alpha = 1 - confidence
      chi2_quantile = qchisq(1 - alpha, df = 2)  # Pour 2D
      
      ellipses_data = data.frame()
      
      for (g in 1:n_components) {
        center_x = means[1, g]
        center_y = means[2, g]
        
        # Matrice de covariance pour cette composante
        if (length(dim(sigmas)) == 3) {
          cov_matrix = sigmas[, , g]
        } else {
          cov_matrix = diag(sigmas[g], 2)
        }
        
        # Décomposition en valeurs/vecteurs propres
        eigen_decomp = eigen(cov_matrix)
        eigenvals = eigen_decomp$values
        eigenvecs = eigen_decomp$vectors
        
        # Demi-axes de l'ellipse
        semi_major = sqrt(eigenvals[1] * chi2_quantile)
        semi_minor = sqrt(eigenvals[2] * chi2_quantile)
        
        # Angle de rotation (en degrés)
        angle_rad = atan2(eigenvecs[2, 1], eigenvecs[1, 1])
        angle_deg = angle_rad * 180 / pi
        
        ellipses_data = rbind(ellipses_data, data.frame(
          component = g,
          center_x = center_x,
          center_y = center_y,
          semi_major = semi_major,
          semi_minor = semi_minor,
          angle = angle_deg,
          stringsAsFactors = FALSE
        ))
      }
      
      return(ellipses_data)
    },

    # Another, non-metric approach: Detect communities or cliques
    detect_graphCommunities = function(method = "SBM") {
      
      if(method == "SBM") {

        mod = BM_bernoulli("SBM", self$G, plotting="", verbosity=0, explore_min=2, explore_max=10)
        mod$estimate()                                # fits K = 2...10 blocks
        Kbest = which.max(mod$ICL)                    # choose by Integrated Completed-Likelihood

        probs   = mod$memberships[[Kbest]]$Z          # assignment of items to blocks
        k = order(t(self$coordinates[,1]) %*% probs)  # cosmetic: order items by their location on the first dimension
        probs = probs[,k]
        if(is.vector(probs)) {
          membership_vec = probs
        } else {
          membership_vec = apply(probs, 1, which.max)
        }

        n_groups = length(unique(membership_vec))
        self$communities = membership_vec
        cat("Detected", n_groups, "communities.\n")
      }

      else if(method == "GMM") {
      }

      else {

        # Create the igraph
        g = graph_from_adjacency_matrix(self$G.bi, mode="undirected", weighted=NULL)
        
        # Apply selected method
        communities = switch(method,
          "louvain" = cluster_louvain(g),
          "walktrap" = cluster_walktrap(g),
          "infomap" = cluster_infomap(g),
          "leiden" = cluster_leiden(g),
          "betweenness" = cluster_edge_betweenness(g),
          stop("Unrecognized cluster method.\n")
        )
        
        # Cluster discovery
        membership_vec = membership(communities)      
        self$communities = membership_vec
        self$modularity = modularity(communities)
        n_groups = length(unique(membership_vec))

        cat("Detected", n_groups, "communities.\n")
        cat("Modularity:", self$modularity, "\n")
      }

      # Compute centroids (efficient for large matrices)
      sums = rowsum(self$coordinates, membership_vec)
      counts = table(membership_vec)
      self$metagraph[["centroids"]] = sums / as.vector(counts)

      # Adjacency matrix between communities
      self$metagraph[["adjacency"]] = matrix(0, n_groups, n_groups)
      rownames(self$metagraph[["adjacency"]]) = paste("Cluster", 1:n_groups, sep="")
      colnames(self$metagraph[["adjacency"]]) = paste("Cluster", 1:n_groups, sep="")
      
      # Count links between each pair of communities
      for (i in 1:nrow(self$G.uni)) {
        for (j in 1:ncol(self$G.uni)) {
          if (self$G.uni[i,j] == 1 && (i != j)) {
            comm_i = membership_vec[i]
            comm_j = membership_vec[j]
            self$metagraph[["adjacency"]][comm_i, comm_j] = self$metagraph[["adjacency"]][comm_i, comm_j] + 1
          }
        }
      }
    },

    # Plot the (log(OR), Delta) decomposition for all pairs
    # This visualizes where each pair falls in terms of symmetric association vs directional asymmetry
    # signed=TRUE keeps the sign of Delta (positive = i->j, negative = j->i)
    # color.by options: "significance", "community", "scale" (requires scales argument)
    #   Features:                                                                                                       
    #   - Plots each unique pair as a point with coordinates (log(OR), |Δ|)                                             
    #   - Color by significance (default): red = unidirectional, blue = bidirectional, grey = not significant           
    #   - Color by community: distinguishes within-community vs cross-community pairs                                   
    #   - Regime boundaries: draws the diagonal line |Δ| = log(OR) where the min-symmetrization switches behavior       
    #   (Regime 1 vs Regime 2 from our theoretical analysis)                                                            
    #   - Labels: optional pair labels with show.labels=TRUE                                                            
    #   - Returns a data frame with all pair data for further analysis                                                  
    plot.decp = function(color.by = "significance", scales = NULL, signed = TRUE,
                         sort.by.difficulty = TRUE,
                         show.labels = FALSE, show.regimes = TRUE,
                         pch = 19, cex = 1, cex.lab = 0.7, alpha = 0.7,
                         col.uni = "#E41A1C", col.bi = "#377EB8", col.ns = "grey60",
                         col.within = "#2CA02C", col.cross = "#9467BD",
                         community.labels = NULL,
                         communities = NULL,
                         show.pair.subsets = FALSE,
                         plot_parabolas = FALSE,
                         main = NULL,
                         xlab = NULL,
                         ylab = "log(OR) - Symmetric association",
                         xlim = NULL, ylim = NULL, ...) {

      # Check that compute() has been run
      if (is.null(self$logOR) || is.null(self$Delta)) {
        stop("Please run compute() first.")
      }

      # Resolve community assignment: user-provided overrides SBM-detected
      comm_vec = if (!is.null(communities)) communities else self$communities

      # Set default title and x-label based on signed option
      if (is.null(main)) {
        main = ifelse(signed,
                      expression(iota^"*"~"Decomposition:"~Delta~"vs log(OR)"),
                      expression(iota^"*"~"Decomposition: |"*Delta*"| vs log(OR)"))
      }
      if (is.null(xlab)) {
        xlab = ifelse(signed,
                      expression(Delta~"- Directional asymmetry"),
                      expression("|"*Delta*"|"~"- Directional asymmetry"))
      }

      n = nrow(self$logOR)

      # Sort items by difficulty if requested (easiest first)
      # This makes i = easier, j = harder in each pair (since upper triangle has i < j)
      if (sort.by.difficulty) {
        item_order = order(self$success.rates)  # increasing success rate = hardest first
        logOR_sorted = self$logOR[item_order, item_order]
        Delta_sorted = self$Delta[item_order, item_order]
        G_uni_sorted = if(!is.null(self$G.uni)) self$G.uni[item_order, item_order] else NULL
        G_bi_sorted  = if(!is.null(self$G.bi))  self$G.bi[item_order, item_order]  else NULL
        names_sorted = self$names[item_order]
        if (!is.null(scales)) scales = scales[item_order]
      } else {
        logOR_sorted = self$logOR
        Delta_sorted = self$Delta
        G_uni_sorted = self$G.uni
        G_bi_sorted  = self$G.bi
        names_sorted = self$names
      }

      # Extract upper triangle (each unique pair once)
      pairs_idx = which(upper.tri(logOR_sorted), arr.ind = TRUE)
      n_pairs = nrow(pairs_idx)

      # Get log(OR) and Delta for each pair
      logOR_vals = numeric(n_pairs)
      Delta_vals = numeric(n_pairs)
      pair_labels = character(n_pairs)
      link_type = character(n_pairs)  # "uni", "bi", or "ns" (not significant)
      direction = character(n_pairs)  # "i->j" or "j->i" for unidirectional links

      for (k in 1:n_pairs) {
        i = pairs_idx[k, 1]
        j = pairs_idx[k, 2]

        logOR_vals[k] = logOR_sorted[i, j]
        Delta_vals[k] = Delta_sorted[i, j]  # With hardest-first sorting: positive Delta = hard->easy (expected Rasch direction); negative Delta = easy->hard (anomalous, see plot annotations).
        pair_labels[k] = paste(names_sorted[i], "-", names_sorted[j])

        # Determine link type based on significance
        if (!is.null(G_uni_sorted) && !is.null(G_bi_sorted)) {
          is_uni_ij = (G_uni_sorted[i, j] == 1)  # i -> j
          is_uni_ji = (G_uni_sorted[j, i] == 1)  # j -> i
          is_bi = (G_bi_sorted[i, j] == 1)

          if (is_uni_ij || is_uni_ji) {
            link_type[k] = "uni"
            direction[k] = ifelse(is_uni_ij, "i->j", "j->i")
          } else if (is_bi) {
            link_type[k] = "bi"
            direction[k] = "both"
          } else {
            link_type[k] = "ns"
            direction[k] = "none"
          }
        } else {
          link_type[k] = "ns"
          direction[k] = "none"
        }
      }

      # Use signed or absolute Delta for plotting (Delta on x-axis)
      x_vals = if (signed) Delta_vals else abs(Delta_vals)
      y_vals = logOR_vals

      # Assign colors based on color.by argument
      if (color.by == "significance") {
        cols = ifelse(link_type == "uni", col.uni,
                     ifelse(link_type == "bi", col.bi, col.ns))
      } else if (color.by == "community" && !is.null(comm_vec)) {
        # Color and symbol by community membership of each pair
        comm_sorted = if(sort.by.difficulty) comm_vec[item_order] else comm_vec
        n_comm = length(unique(comm_sorted))

        # Palette and symbols for communities (10 entries to support up to 10 communities)
        comm_palette = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                         "#A65628", "#F781BF", "#999999", "#1B9E77", "#E6AB02")[1:n_comm]
        comm_pchs = c(19, 17, 15, 18, 8, 3, 4, 1, 6, 5)[1:n_comm]

        # For each pair: assign color/symbol based on community pair
        pair_comm_i = comm_sorted[pairs_idx[, 1]]
        pair_comm_j = comm_sorted[pairs_idx[, 2]]
        within_comm = (pair_comm_i == pair_comm_j)

        # Within-community pairs: color/symbol of that community
        # Cross-community pairs: grey, open circle
        cols = ifelse(within_comm, comm_palette[pair_comm_i], "grey60")
        pch_vals = ifelse(within_comm, comm_pchs[pair_comm_i], 1)
      } else if (color.by == "scale") {
        # Color by whether pair is within-scale or cross-scale (requires scales argument)
        if (is.null(scales)) {
          stop("color.by='scale' requires the 'scales' argument (vector of scale membership for each item)")
        }
        if (length(scales) != n) {
          stop(paste("scales vector length (", length(scales), ") must match number of items (", n, ")", sep=""))
        }
        within_scale = sapply(1:n_pairs, function(k) {
          i = pairs_idx[k, 1]
          j = pairs_idx[k, 2]
          scales[i] == scales[j]
        })
        cols = ifelse(within_scale, col.within, col.cross)
      } else {
        cols = rep(col.ns, n_pairs)
      }

      # Add transparency
      cols = adjustcolor(cols, alpha.f = alpha)

      # Set axis limits from data range with a small margin
      if (is.null(xlim)) {
        xr = range(x_vals, na.rm = TRUE)
        xmargin = diff(xr) * 0.1
        xlim = c(xr[1] - xmargin, xr[2] + xmargin)
      }
      if (is.null(ylim)) {
        yr = range(y_vals, na.rm = TRUE)
        ymargin = diff(yr) * 0.1
        ylim = c(yr[1] - ymargin, yr[2] + ymargin)
      }

      # Create empty plot first, then add grid, then points on top
      plot(x_vals, y_vals, type = "n",
           xlim = xlim, ylim = ylim,
           xlab = xlab, ylab = ylab, main = main, ...)

      # Add grid first (background)
      grid(col = "grey90", lty = 1)

      # Add reference lines
      abline(h = 0, col = "grey80", lty = 2)
      abline(v = 0, col = "grey80", lty = 2)

      # Show non-significance bands
      if (show.regimes) {
        # Compute band half-widths from median SE and critical values
        se_or_vals = self$SE_OR[lower.tri(self$SE_OR)]
        se_asym_vals = self$SE_asym[lower.tri(self$SE_asym)]
        median_se_or = median(se_or_vals, na.rm = TRUE)
        median_se_asym = median(se_asym_vals, na.rm = TRUE)

        # One-sided z-critical for association, two-sided for directionality
        z_crit_or   = qnorm(1 - self$alpha)       # e.g. 1.645
        z_crit_asym = qnorm(1 - self$alpha / 2)   # e.g. 1.96

        band_or   = z_crit_or * median_se_or
        band_asym = z_crit_asym * median_se_asym

        # Draw horizontal n.s. association band
        rect(xlim[1], -band_or, xlim[2], band_or,
             col = adjustcolor("grey", alpha.f = 0.15), border = NA)
        text(xlim[1] * 0.7, band_or * 0.4,
             "n.s. association", cex = 0.7, col = "grey50", font = 3)

        # Draw vertical n.s. directionality band
        rect(-band_asym, ylim[1], band_asym, ylim[2],
             col = adjustcolor("grey", alpha.f = 0.15), border = NA)
        text(band_asym * 0.6, ylim[1] * 0.7,
             "n.s. directionality", cex = 0.7, col = "grey50", font = 3, srt = 90)

        if (signed) {
          # Direction labels (at the bottom of the graph). Anchor each label to its
          # own side so asymmetric xlim ranges (e.g., c(-2, 8)) don't push the
          # negative-side label outside the visible plot area.
          y_pos = ylim[1] + 0.05 * diff(ylim)
          x_pos_right = xlim[2] - 0.05 * diff(xlim)
          x_pos_left  = xlim[1] + 0.05 * diff(xlim)
          if (sort.by.difficulty) {
            text(x_pos_right, y_pos,
                 expression(hard %->% easy),
                 cex = 0.9, col = "grey40", font = 2, adj = c(1, 0))
            text(x_pos_left, y_pos,
                 expression(easy %->% hard),
                 cex = 0.9, col = "grey40", font = 2, adj = c(0, 0))
          } else {
            text(x_pos_right, y_pos,
                 expression(i %->% j),
                 cex = 0.9, col = "grey40", font = 2, adj = c(1, 0))
            text(x_pos_left, y_pos,
                 expression(j %->% i),
                 cex = 0.9, col = "grey40", font = 2, adj = c(0, 0))
          }
        }
      }

      # Draw PCA first-component lines for each oriented community-pair subset
      if (show.pair.subsets && !is.null(comm_vec)) {
        comm_sorted_ps = if(sort.by.difficulty) comm_vec[item_order] else comm_vec

        # Build oriented pair type label for each unidirectional link
        pair_type = character(n_pairs)
        for (k in 1:n_pairs) {
          if (link_type[k] != "uni") next
          ci = comm_sorted_ps[pairs_idx[k, 1]]
          cj = comm_sorted_ps[pairs_idx[k, 2]]
          if (direction[k] == "i->j") {
            pair_type[k] = paste(ci, "->", cj)
          } else {
            pair_type[k] = paste(cj, "->", ci)
          }
        }

        # For each unique oriented pair type, compute PCA and draw PC1
        uni_mask = (link_type == "uni")
        types = unique(pair_type[uni_mask])

        for (tp in types) {
          idx = which(uni_mask & pair_type == tp)
          if (length(idx) < 2) next  # Need at least 2 points for PCA

          xy = cbind(x_vals[idx], y_vals[idx])
          pc = prcomp(xy, center = TRUE, scale. = FALSE)
          center = pc$center
          v1 = pc$rotation[, 1]  # First principal component direction

          # Project points onto PC1 to get range
          proj = (xy[, 1] - center[1]) * v1[1] + (xy[, 2] - center[2]) * v1[2]
          t_range = range(proj)

          # Extend slightly beyond the data
          t_range = t_range + c(-0.15, 0.15) * diff(t_range)

          # Line endpoints
          x0 = center[1] + t_range[1] * v1[1]
          y0 = center[2] + t_range[1] * v1[2]
          x1 = center[1] + t_range[2] * v1[1]
          y1 = center[2] + t_range[2] * v1[2]

          segments(x0, y0, x1, y1, col = adjustcolor("grey50", alpha.f = 0.5),
                   lwd = 1.5, lty = 1)
        }
      }

      # Add points on top of grid and reference lines (skip if glyphs/numbers will be drawn)
      if (!identical(show.labels, "clusters") && !identical(show.labels, "numbers")) {
        if (exists("pch_vals")) {
          points(x_vals, y_vals, pch = pch_vals, cex = cex, col = cols)
        } else {
          points(x_vals, y_vals, pch = pch, cex = cex, col = cols)
        }
      }

      # Add pair labels if requested
      if (identical(show.labels, TRUE)) {
        text(x_vals, y_vals, labels = pair_labels,
             pos = 3, cex = cex.lab, col = adjustcolor("black", alpha.f = 0.7))
      } else if (identical(show.labels, "clusters") && !is.null(comm_vec)) {
        # Cluster glyph labels: two community symbols linked by segment or arrow
        comm_sorted = if(sort.by.difficulty) comm_vec[item_order] else comm_vec
        n_comm = length(unique(comm_sorted))
        cl_palette = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                        "#A65628", "#F781BF", "#999999", "#1B9E77", "#E6AB02")[1:n_comm]
        cl_pchs = c(19, 17, 15, 18, 8, 3, 4, 1, 6, 5)[1:n_comm]

        # Offset in user coordinates (fraction of plot range)
        dx = diff(xlim) * 0.015
        glyph_cex = cex * 0.7

        for (k in 1:n_pairs) {
          if (link_type[k] == "ns") next  # Skip non-significant pairs

          xi = x_vals[k]
          yi = y_vals[k]
          ci = comm_sorted[pairs_idx[k, 1]]
          cj = comm_sorted[pairs_idx[k, 2]]

          # Draw left symbol (item i) and right symbol (item j)
          points(xi - dx, yi, pch = cl_pchs[ci], cex = glyph_cex,
                 col = cl_palette[ci])
          points(xi + dx, yi, pch = cl_pchs[cj], cex = glyph_cex,
                 col = cl_palette[cj])

          # Connect with arrow (unidirectional) or segment (bidirectional)
          if (link_type[k] == "uni") {
            if (direction[k] == "i->j") {
              arrows(xi - dx * 0.5, yi, xi + dx * 0.5, yi,
                     length = 0.04, col = "grey30", lwd = 0.8)
            } else {
              arrows(xi + dx * 0.5, yi, xi - dx * 0.5, yi,
                     length = 0.04, col = "grey30", lwd = 0.8)
            }
          } else if (link_type[k] == "bi") {
            segments(xi - dx * 0.5, yi, xi + dx * 0.5, yi,
                     col = "grey30", lwd = 0.8)
          }
        }
      } else if (identical(show.labels, "numbers") && !is.null(comm_vec)) {
        # Item number pairs, colored by cluster, linked by arrow or segment
        comm_sorted = if(sort.by.difficulty) comm_vec[item_order] else comm_vec
        n_comm = length(unique(comm_sorted))
        cl_palette = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                        "#A65628", "#F781BF", "#999999", "#1B9E77", "#E6AB02")[1:n_comm]

        # Extract numeric item indices (original, before sorting)
        if (sort.by.difficulty) {
          orig_idx = item_order  # maps sorted position -> original item number
        } else {
          orig_idx = 1:n
        }

        dx = diff(xlim) * 0.02
        num_cex = cex.lab * 0.9

        for (k in 1:n_pairs) {
          if (link_type[k] == "ns") next

          xi = x_vals[k]
          yi = y_vals[k]
          si = pairs_idx[k, 1]  # sorted index
          sj = pairs_idx[k, 2]
          ci = comm_sorted[si]
          cj = comm_sorted[sj]

          # Item numbers in original ordering
          ni = orig_idx[si]
          nj = orig_idx[sj]

          # Draw item numbers
          text(xi - dx, yi, labels = ni, cex = num_cex, col = cl_palette[ci], font = 2)
          text(xi + dx, yi, labels = nj, cex = num_cex, col = cl_palette[cj], font = 2)

          # Connect with arrow (unidirectional) or segment (bidirectional)
          if (link_type[k] == "uni") {
            if (direction[k] == "i->j") {
              arrows(xi - dx * 0.4, yi, xi + dx * 0.4, yi,
                     length = 0.03, col = "grey50", lwd = 0.6)
            } else {
              arrows(xi + dx * 0.4, yi, xi - dx * 0.4, yi,
                     length = 0.03, col = "grey50", lwd = 0.6)
            }
          } else if (link_type[k] == "bi") {
            segments(xi - dx * 0.4, yi, xi + dx * 0.4, yi,
                     col = "grey50", lwd = 0.6)
          }
        }
      }

      # Per-community quadratic fits log(OR) = a + b * Delta^2 on significant
      # within-community pairs, plus a cross-community reference fit. Drawn on
      # top of the cloud, after points and glyphs.
      if (plot_parabolas && color.by == "community" && !is.null(comm_vec)) {
        comm_for_fit = if (sort.by.difficulty) comm_vec[item_order] else comm_vec
        n_comm_fit = length(unique(comm_for_fit))
        fit_palette = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                        "#A65628", "#F781BF", "#999999", "#1B9E77", "#E6AB02")[1:n_comm_fit]

        ci = comm_for_fit[pairs_idx[, 1]]
        cj = comm_for_fit[pairs_idx[, 2]]
        sig = (link_type != "ns")
        x_grid = seq(min(x_vals, na.rm = TRUE), max(x_vals, na.rm = TRUE), length.out = 200)

        for (c in seq_len(n_comm_fit)) {
          idx = which(sig & ci == c & cj == c)
          if (length(idx) < 5) next
          fit = try(lm(y_vals[idx] ~ I(x_vals[idx]^2)), silent = TRUE)
          if (inherits(fit, "try-error")) next
          y_pred = predict(fit, newdata = data.frame(`I(x_vals[idx]^2)` = x_grid^2,
                                                     check.names = FALSE))
          # Use direct coefficient evaluation since predict() may complain about names
          a = coef(fit)[1]; b = coef(fit)[2]
          y_pred = a + b * x_grid^2
          lines(x_grid, y_pred, col = fit_palette[c], lwd = 2.4)
        }

        # Cross-community reference fit on significant pairs
        idx_x = which(sig & ci != cj)
        if (length(idx_x) >= 5) {
          fit_x = lm(y_vals[idx_x] ~ I(x_vals[idx_x]^2))
          a = coef(fit_x)[1]; b = coef(fit_x)[2]
          lines(x_grid, a + b * x_grid^2, col = "grey30", lwd = 2, lty = 2)
        }
      }

      # Helper: draw a rounded-corner box behind a legend
      legend_rounded = function(...) {
        # Draw legend invisibly to get bounding box
        lg = legend(..., plot = FALSE)
        r  = lg$rect
        pad = strwidth("M") * 0.3
        x0 = r$left - pad;  x1 = r$left + r$w + pad
        y0 = r$top - r$h - pad;  y1 = r$top + pad
        # Rounded rectangle via polygon
        rd = min(r$w, r$h) * 0.15  # corner radius
        angles = seq(0, pi/2, length.out = 20)
        cx = c(x1 - rd + rd*cos(angles),        # top-right
               x0 + rd - rd*cos(rev(angles)),    # top-left
               x0 + rd - rd*cos(angles),         # bottom-left
               x1 - rd + rd*cos(rev(angles)))    # bottom-right
        cy = c(y1 - rd + rd*sin(angles),
               y1 - rd + rd*sin(rev(angles)),
               y0 + rd - rd*sin(angles),
               y0 + rd - rd*sin(rev(angles)))
        polygon(cx, cy, col = adjustcolor("grey95", alpha.f = 0.85),
                border = "grey95", lwd = 0.6)
        # Re-draw the legend on top
        legend(..., plot = TRUE)
      }

      # Add legend
      if (color.by == "significance") {
        n_uni = sum(link_type == "uni")
        n_bi = sum(link_type == "bi")
        n_ns = sum(link_type == "ns")
        legend_rounded("topright", inset=.02, bty="n",
               legend = c(paste("Unidirectional (", n_uni, ")", sep = ""),
                         paste("Bidirectional (", n_bi, ")", sep = ""),
                         paste("Not significant (", n_ns, ")", sep = "")),
               col = c(col.uni, col.bi, col.ns),
               pch = pch, pt.cex = cex, cex = 0.8)
      } else if (color.by == "community") {
        comm_labels = if (!is.null(community.labels)) community.labels else paste("Community", 1:n_comm)
        if (identical(show.labels, "clusters") || identical(show.labels, "numbers")) {
          legend_rounded("topright", inset=.02, bty="n",
                 legend = comm_labels,
                 col = comm_palette,
                 pch = comm_pchs,
                 pt.cex = cex, cex = 0.8)
        } else {
          legend_rounded("topright", inset=.02, bty="n",
                 legend = c(comm_labels, "Cross-community"),
                 col = c(comm_palette, "grey60"),
                 pch = c(comm_pchs, 1),
                 pt.cex = cex, cex = 0.8)
        }
      } else if (color.by == "scale") {
        n_within = sum(within_scale)
        n_cross = sum(!within_scale)
        legend_rounded("topright", inset=.02, bty="n",
               legend = c(paste("Within-scale (", n_within, ")", sep = ""),
                         paste("Cross-scale (", n_cross, ")", sep = "")),
               col = c(col.within, col.cross),
               pch = pch, pt.cex = cex, cex = 0.8)
      }

      # Build return data frame
      result_df = data.frame(
        i = pairs_idx[, 1],
        j = pairs_idx[, 2],
        pair = pair_labels,
        logOR = logOR_vals,
        Delta = Delta_vals,
        link_type = link_type,
        direction = direction,
        stringsAsFactors = FALSE
      )

      # Add scale information if provided
      if (!is.null(scales)) {
        result_df$scale_i = scales[pairs_idx[, 1]]
        result_df$scale_j = scales[pairs_idx[, 2]]
        result_df$within_scale = within_scale
      }

      # Return data invisibly for further analysis
      invisible(result_df)
    },

    # Affichage du graphe dans une fenêtre R
    plot.diagram = function(adjacency=FALSE, plot.coefs = FALSE, plot.rates = TRUE, plot.arrows = TRUE, plot.axes = FALSE, pos = NULL, curve = 0.2, name = NULL, absent = 0, relsize = 0.8, lwd = 1, lcol = "black", box.size = 0.1, box.type = "circle", box.prop = 0.5, box.col = "white", box.lcol = lcol, box.lwd = lwd, shadow.size = 0.01, shadow.col = "grey", dr = 0.01, dtext = 0.3, self.lwd = 1, self.cex = 1, self.shiftx = box.size, self.shifty = NULL, self.arrpos = NULL, arr.lwd = lwd, arr.lcol = lcol, arr.tcol = lcol, arr.col = "black", arr.type = "curved", arr.pos = 0.5, arr.length = 0.4, arr.width = arr.length/2, endhead = FALSE, mx = 0.0, my = 0.0, box.cex = 1, txt.col = "black", txt.xadj = 0.5, txt.yadj = 0.5, txt.font = 1, prefix = "", cex = 1, cex.txt = cex, add = FALSE, main = "", cex.main = cex, segment.from = 0, segment.to = 1, latex = FALSE, ...) {

      g = matrix(sprintf("%.2f",self$Intensity), nrow=nrow(self$Intensity), ncol=ncol(self$Intensity))

      if(is.null(name)) name = self$names
      if(plot.rates) name = paste(name,"\n(",round(self$success.rates,2),")", sep="")

      # An adjacency representation is more appropriate for the non-cumulative case (DINA, unfolding)
      if(adjacency) {

        # Store the graph in igraph format
        ig = graph_from_adjacency_matrix(self$G)
        self$iGraph = ig

        # Compute the coordinates by force-directed layout
        coord = layout_with_fr(ig,dim=2)
        coord = self$normalize_coordinates(scale(coord))
        self$coordinates = coord
      }

      # Don't plot the coefficients by default
      if(!plot.coefs) {
        g[,] = ""
      }

      # Hide all links when arrows are disabled
      if(!plot.arrows) {
        g[,] = NA
      } else {
        # Don't plot the zero-valued links by default
        g[self$G == 0] = NA
      }

      # diagram wants coordinates in [0;1]
      normCoord = self$normalize_coordinates(self$coordinates[,1:2])

      # Draw axes before nodes so they appear behind
      if(plot.axes) {
        coords = self$coordinates[,1:2]
        x_min = min(coords[,1]); x_max = max(coords[,1])
        y_min = min(coords[,2]); y_max = max(coords[,2])
        x_range = x_max - x_min; y_range = y_max - y_min
        max_range = max(x_range, y_range)
        x_offset = (1 - x_range / max_range) / 2
        y_offset = (1 - y_range / max_range) / 2

        # Helper: map original coordinate to normalized [0,1] space
        to_norm_x = function(v) (v - x_min) / max_range + x_offset
        to_norm_y = function(v) (v - y_min) / max_range + y_offset

        # Origin in normalized space
        x0 = to_norm_x(0)
        y0 = to_norm_y(0)

        # Set up empty plot with plotmat's coordinate system
        plot(0, 0, type = "n", xlim = c(-0.05, 1.05), ylim = c(-0.05, 1.05),
             xlab = "", ylab = "", axes = FALSE, asp = 1, main = main, cex.main = cex.main)

        # Draw axis lines with arrows
        usr = par("usr")
        arrows(usr[1], y0, usr[2], y0, col = "grey60", lwd = 1, length = 0.1)
        arrows(x0, usr[3], x0, usr[4], col = "grey60", lwd = 1, length = 0.1)

        # Graduation: choose nice tick positions in original space
        x_ticks = pretty(c(x_min, x_max), n = 6)
        x_ticks = x_ticks[x_ticks > x_min & x_ticks < x_max]
        y_ticks = pretty(c(y_min, y_max), n = 6)
        y_ticks = y_ticks[y_ticks > y_min & y_ticks < y_max]

        tick_len = 0.01
        for(xt in x_ticks) {
          nx = to_norm_x(xt)
          segments(nx, y0 - tick_len, nx, y0 + tick_len, col = "grey60")
          text(nx, y0 - 2.5 * tick_len, labels = xt, cex = 0.6, col = "grey60")
        }
        for(yt in y_ticks) {
          ny = to_norm_y(yt)
          segments(x0 - tick_len, ny, x0 + tick_len, ny, col = "grey60")
          text(x0 - 2.5 * tick_len, ny, labels = yt, cex = 0.6, col = "grey60")
        }

        # Force plotmat to draw on top of axes
        add = TRUE
      }

      plotmat(t(g), pos=normCoord, curve = curve, name = name, absent = absent, relsize = relsize, lwd = lwd, lcol = lcol, box.size = box.size, box.type = box.type, box.prop = box.prop, box.col = box.col, box.lcol = lcol, box.lwd = lwd, shadow.size = shadow.size, shadow.col = shadow.col, dr = dr, dtext = dtext, self.lwd = self, self.cex = self.cex, self.shiftx = self.shiftx, self.shifty = self.shifty, self.arrpos = self.arrpos, arr.lwd = arr.lwd, arr.lcol = arr.lcol, arr.tcol = arr.tcol, arr.col = arr.col, arr.type = arr.type, arr.pos = arr.pos, arr.length = arr.length, arr.width = arr.width, endhead = endhead, mx = mx, my = my, box.cex = box.cex, txt.col = txt.col, txt.xadj = txt.xadj, txt.yadj = txt.yadj, txt.font = txt.font, prefix = prefix, cex = cex, cex.txt = cex.txt, add = add, main = main, cex.main = cex.main, segment.from = segment.from, segment.to = segment.to, latex = latex, ...)

    },

    # Normalisation des coordonnées du graphe entre 0 et 1 (pour le package diagram)
    normalize_coordinates = function(coords) {

      # Vérifier que l'entrée est une matrice à 2 colonnes
      if (!is.matrix(coords) || ncol(coords) != 2) {
        stop("Input argument should be a 2-column matrix.")
      }
      
      if (nrow(coords) == 0) {
        return(coords)
      }
      
      # Extraire les coordonnées x et y
      x = coords[, 1]
      y = coords[, 2]
      
      # Calculer les valeurs min et max pour chaque dimension
      x_min = min(x, na.rm = TRUE)
      x_max = max(x, na.rm = TRUE)
      y_min = min(y, na.rm = TRUE)
      y_max = max(y, na.rm = TRUE)
      
      # Calculer les étendues
      x_range = x_max - x_min
      y_range = y_max - y_min
      
      # Utiliser l'étendue maximale pour préserver l'aspect ratio
      max_range = max(x_range, y_range)
      
      # Éviter la division par zéro
      if (max_range == 0) {
        # Si tous les points sont identiques, les placer au centre
        normalized_coords = matrix(0.5, nrow = nrow(coords), ncol = 2)
        colnames(normalized_coords) = colnames(coords)
        return(normalized_coords)
      }
      
      # Normaliser les coordonnées
      x_normalized = (x - x_min) / max_range
      y_normalized = (y - y_min) / max_range
      
      # Centrer les coordonnées dans [0,1] pour que l'ensemble soit bien positionné
      x_center_offset = (1 - x_range / max_range) / 2
      y_center_offset = (1 - y_range / max_range) / 2
      
      x_final = x_normalized + x_center_offset
      y_final = y_normalized + y_center_offset
      
      # Créer la matrice de sortie
      normalized_coords = cbind(x_final, y_final)
      colnames(normalized_coords) = colnames(coords)
      
      return(normalized_coords)
    },

    # Find a subset with only one-way implications (cumulative structure)
    findCumulativeSubset = function(tol = 0) {

      ## ── helpers ──────────────────────────────────────────────────────────────
      pattern_ok = function(m) {
        upper = m[row(m) < col(m)]        # cells above the main diagonal
        lower = m[row(m) > col(m)]        # cells below the main diagonal
        all(upper >  tol) && all(lower < -tol)
      }

      violation_counts = function(m) {
        n  = nrow(m)
        v  = numeric(n)                   # violations attributable to each variable
        up = row(m) <  col(m)
        lo = row(m) >  col(m)

        viol = matrix(FALSE, n, n)
        viol[up] = m[up] <=  tol          # should be > +tol  above diag
        viol[lo] = m[lo] >= -tol          # should be < -tol  below diag

        # count violations appearing in a row *or* column of each variable
        for (k in seq_len(n)) v[k] = sum(viol[k, ]) + sum(viol[ ,k])
        v
      }

      ## ── set-up ───────────────────────────────────────────────────────────────
      k = order(self$success.rates)
      current   = self$IOTA_STAR[k,k]
      removed   = character(0)

      ## ── iterative pruning ────────────────────────────────────────────────────
      while (nrow(current) > 1 && !pattern_ok(current)) {
        vcount = violation_counts(current)
        idx    = which.max(vcount)                       # worst offender
        removed = c(removed, colnames(current)[idx])
        current = current[-idx, -idx, drop = FALSE]      # peel it off
      }

      list(filtered_matrix=current, removed_vars=removed)
    },
    #------------------------------ Estimate Y-axis band levels via clustering -------------------------
    # method: "mclust" (Gaussian mixture, BIC-based) or "kmeans" (silhouette-based)
    # max_k:  maximum number of levels to consider
    # k:      if provided, forces this number of levels (no automatic search)
    cluster_logOR_levels = function(method = "kmeans", max_k = 9, k = NULL) {

      y = self$coordinates[, 2]

      if (method == "mclust") {

        fit = Mclust(y, G = if (!is.null(k)) k else 1:max_k, modelNames = "V")
        best_k = fit$G
        membership = fit$classification
        centers = sort(fit$parameters$mean)

        self$logOR_levels = list(
          centers = centers,
          k = best_k,
          membership = membership,
          method = "mclust",
          fit = fit
        )

      } else if (method == "kmeans") {

        if (!is.null(k)) {
          # Fixed k
          best_k = k
        } else {
          # Automatic search via average silhouette width
          sil_scores = numeric(max_k - 1)
          for (kk in 2:max_k) {
            km = kmeans(y, centers = kk, nstart = 25)
            sil = cluster::silhouette(km$cluster, dist(y))
            sil_scores[kk - 1] = mean(sil[, 3])
          }
          best_k = which.max(sil_scores) + 1
        }

        km = kmeans(y, centers = best_k, nstart = 25)
        centers = sort(km$centers[, 1])
        membership = km$cluster

        self$logOR_levels = list(
          centers = centers,
          k = best_k,
          membership = membership,
          method = "kmeans",
          silhouettes = if (is.null(k)) sil_scores else NULL,
          fit = km
        )
      }

      invisible(self)
    },

    #----------------------------- Multidimensional scaling with missing distances imputation ----------
    iterative_mds_em = function(dist_matrix, max_iter = 5000, tol = 1e-6, init_method = "mean", k, eig=FALSE) {

      # Initialize the missing values in the distance matrix
      dist_init = dist_matrix
      
      missing.values = is.na(dist_matrix)

      # Step 1: Initialization
      if (init_method == "mean") {
        dist_init[missing.values] = mean(dist_matrix, na.rm = TRUE)
      } else if (init_method == "random") {
        dist_init[missing.values] = runif(sum(missing.values), min(dist_matrix, na.rm = TRUE), max(dist_matrix, na.rm = TRUE))
      }

      # Step 2: Initialize points via MDS
      mds_result = cmdscale(dist_init, k = k, eig=eig)
    
      if(eig) points = mds_result$points
      else points = mds_result
      
      prev_stress = Inf
      
      # Step 3: Iterate until convergence
      for (i in 1:max_iter) {

        # Calculate the distances from current points (M-step)
        d_hat = as.matrix(dist(points))

        # Replace missing distances with the estimated ones (E-step)
        dist_init[missing.values] = d_hat[missing.values]
        
        # Re-run MDS with the updated distances
        mds_result = cmdscale(dist_init, k = k, eig=eig)

        if(eig) points = mds_result$points
        else points = mds_result
        
        # Calculate current stress (can use other stress functions as well)
        stress = sum((d_hat - dist_init)^2, na.rm = TRUE)
        
        # Check for convergence
        if (abs(prev_stress - stress) < tol) {
          break
        }
        
        prev_stress = stress
      }
      
      cat("MDS converged at iteration", i, "\n")

      # Return the final configuration of points
      return(mds_result)
    },

    #----------------------------------- Compute a bounding curve around a set of points
    bounding_curve = function(points, dilate_factor = 1.5) {

      # Find the convex hull (note: chull() manage the ordering)
      hull_indices = chull(points[,1], points[,2])
      hull_points = points[hull_indices, ]

      # Compute the centroid of the original points
      centroid_x = mean(points[,1])
      centroid_y = mean(points[,2])

      # Dilate the points
      dx = hull_points[,1] - centroid_x
      dy = hull_points[,2] - centroid_y
      new_x = centroid_x + dx * dilate_factor
      new_y = centroid_y + dy * dilate_factor

      return(data.frame(x = new_x, y = new_y))

      # Trick to fake a graphical context (needed for xspline)
      fn = tempfile(fileext = ".png")
      grDevices::png(filename = fn, width = 200, height = 200, units = "px", bg = "transparent")
      on.exit({ try(grDevices::dev.off(), silent = TRUE); unlink(fn) }, add = TRUE)

      # >>> Élimine les marges et utilise tout le device
      op = graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(op), add = TRUE)
      graphics::par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i", plt = c(0, 1, 0, 1))

      # Crée la région de tracé
      graphics::plot.new()
      xr = range(new_x, finite = TRUE); yr = range(new_y, finite = TRUE)

      # Padding si étendue nulle (évite erreur fenêtre de tracé)
      if (!is.finite(diff(xr)) || diff(xr) == 0) xr = xr + c(-1, 1) * 1e-6
      if (!is.finite(diff(yr)) || diff(yr) == 0) yr = yr + c(-1, 1) * 1e-6
      graphics::plot.window(xlim = xr, ylim = yr)

      # draw=FALSE to get the points
      pts = xspline(new_x, new_y, shape = 1, open = FALSE, draw=FALSE)
      pts
    },

    #----------------------------------- Export graph to Tikz/LaTeX file

    plot.tikz = function(filename="DiG_plot.tex", adjacency=FALSE, skills=NULL, x.skills=NULL, y.skills=NULL, Q=NULL, count_processes=FALSE, node_groups=NULL, threshold=0.95, scale_factor=5, bend_angle=20, bend_right=NULL, meta.bend_right=NULL, meta.bend_angle=20, xlim=NULL, ylim=NULL, mainNodeWidth="2cm", mainNodeHeight="1cm", mainNodeTextWidth="1.5cm", 
    plot_coefs=FALSE, plot_band=FALSE, plot_success=TRUE, plot_communities=FALSE,  plot_ellipses=FALSE, plot_metagraph=FALSE, k=NULL, community_type="convexhull", dim1=1, dim2=2, legend=NULL, legend_title="Items", legend_x=0.75, legend_y=2.6, color_univoque="black!70", color_biunivoque="gray!60", width_univoque="1.5pt", width_biunivoque="2.5pt", opacity_univoque=0.9, opacity_biunivoque=1, only_directed=FALSE, only_symmetric=FALSE, confidence=0.95, from.list=NULL, to.list=NULL, only_internals=FALSE, dilate_cluster=1, cluster_labels=NULL, cluster_labels_below=FALSE, plot_logOR_levels=FALSE, plot_difficulty_gradient=FALSE) {
  
      n_items = nrow(self$IOTA_STAR)
 
      # Names and labels
      item_names = self$names
      if(is.null(item_names)) {
        item_names = rownames(self$IOTA_STAR)
      }

      if(is.null(item_names)) {
        item_names = paste0("Item", 1:n_items)
      }
      
      # Escape underscores for LaTeX display text
      tex_labels = gsub("_", "\\_", item_names, fixed = TRUE)

      if(!is.null(self$success.rates) && plot_success) {
        item_labels = paste(tex_labels," (",round(self$success.rates,2),")", sep="")
      } else {
        item_labels = tex_labels
      }

      if(!is.null(Q) && count_processes) {
        n_processes = rowSums(Q)
        item_labels = paste(tex_labels,"\\\\(",n_processes,")", sep="")
      }

      # Switch to adjacency representation if requested
      if(adjacency) {

        # Store the graph in igraph format
        ig = graph_from_adjacency_matrix(self$G)
        self$iGraph = ig
      
        # Compute the coordinates by force-directed layout
        coord = layout_with_fr(ig,dim=2)
        coord = self$normalize_coordinates(scale(coord))
        self$coordinates = coord
      }

      # TikZ code
      tikz_code = sprintf("\\begin{tikzpicture}[
        scale=%.2f,
        thick,
        main node/.style={ellipse, draw, font=\\small\\bfseries, minimum height=%s, minimum width=%s, inner sep=2pt, fill=white, fill opacity=0.7, text width=%s, align=center},
        skill node/.style={circle, fill=gray!75, minimum size=2mm, inner sep=2pt, label={[text width=2.5cm, align=center, text=black!80, font=\\large\\bfseries, fill=white, fill opacity=0.8,label distance=3mm]above:#1}, fill opacity=1},
        group/.style={draw=gray!50, fill=gray!20, rounded corners=5mm, inner sep=3mm},
        gaussian ellipse/.style={draw=blue!60, fill=blue!10, line width=1pt, opacity=0.7},
        community ellipse/.style={draw=gray!70, fill=lightgray!50, line width=1pt, opacity=0.6},
        auto-ellipse/.style={ellipse, draw=gray!70, fill=lightgray!50, line width=1pt, opacity=0.5, inner sep=0pt},
        arrow/.style={->, -{Stealth[length=3mm, width=2mm]},line width=%s, %s, opacity=%.2f, behind path},
        double arrow/.style={=>, {Stealth[length=3mm, width=2mm]}-{Stealth[length=3mm, width=2mm]}, line width=%s, %s, opacity=%.2f, behind path},
        meta arrow/.style={->, -{Stealth[length=7mm, width=5mm]}, line width=6pt, gray!50, opacity=0.7, behind path},
        weight label/.style={midway, font=\\scriptsize, fill=white, fill opacity=1, text opacity=1, inner sep=2pt, outer sep=0pt, rectangle, draw=none},
        meta label/.style={text=black, midway, font=\\large, fill=white, fill opacity=1, text opacity=1, inner sep=2pt, outer sep=0pt, rectangle, draw=none},
        axis style/.style={-, thick, gray!80},
        xtick node/.style={below=0.5mm, font=\\small, text=black!60},
        ytick node/.style={left=0.5mm, font=\\small, text=black!60}
      ]\n\n",scale_factor, mainNodeHeight, mainNodeWidth, mainNodeTextWidth, width_univoque, color_univoque, opacity_univoque, width_biunivoque, color_biunivoque, opacity_biunivoque)

      # Filter visible items based on xlim/ylim
      visible_items = 1:n_items
      if (!is.null(xlim)) {
        visible_items = visible_items[self$coordinates[visible_items, dim1] >= xlim[1] &
                                      self$coordinates[visible_items, dim1] <= xlim[2]]
      }
      if (!is.null(ylim)) {
        visible_items = visible_items[self$coordinates[visible_items, dim2] >= ylim[1] &
                                      self$coordinates[visible_items, dim2] <= ylim[2]]
      }

      # Defining nodes (only visible items)
      tikz_code = paste0(tikz_code, "% Node definitions\n")
      for (i in visible_items) {
        tikz_code = paste0(tikz_code,
                          sprintf("\\node[main node] (%s) at (%.2f,%.2f) {%s};\n",
                                  item_names[i],
                                  self$coordinates[i,dim1],
                                  self$coordinates[i,dim2],
                                  item_labels[i]))
      }
      
      # Define the skills
      if (!is.null(skills)) {

        n_skills = length(skills)
        skill_code = names(skills)
        if(is.null(skill_code)) skill_code = skills

        if(length(x.skills) == 1) x.skills = rep(x.skills,n_skills)
        if(length(y.skills) == 1) y.skills = rep(y.skills,n_skills)
        
        tikz_code = paste0(tikz_code, "\n% Skills\n")
        for (i in 1:n_skills) {
          # Filter skills based on xlim/ylim
          in_xlim = is.null(xlim) || (x.skills[i] >= xlim[1] && x.skills[i] <= xlim[2])
          in_ylim = is.null(ylim) || (y.skills[i] >= ylim[1] && y.skills[i] <= ylim[2])
          if (in_xlim && in_ylim) {
            tikz_code = paste0(tikz_code,
                                sprintf("\\node[skill node=%s] (%s) at (%.2f,%.2f) {};\n",
                                        skills[i],
                                        skill_code[i],
                                        x.skills[i],
                                        y.skills[i]))
          }
        }
      }

      # Collector for cluster labels (printed last, on top of everything)
      cluster_labels_code = ""

      # Background layer
      tikz_code = paste0(tikz_code,"\n\\begin{pgfonlayer}{background}\n")

      # Define groups of nodes, if relevant
      if(!is.null(node_groups)) {
        tikz_code = paste0(tikz_code, "\n% Groups\n")
        n_groups = length(node_groups)
        group_names = names(node_groups)
        group_code = paste("Grp", 1:n_groups, sep="")
        for (i in 1:n_groups) {
          group_nodes = paste("(",node_groups[[i]],")",sep="",collapse=" ")
          tikz_code = paste0(tikz_code,
                            # \node[group, fit=(Item_A) (Item_B) (Item_C), label={[font=\footnotesize]above:Groupe Alpha}] (alpha_group) {};
                            sprintf("\\node[group, fit=%s, label={[font=\\footnotesize]above:%s}] (%s) {};;\n",
                                    group_nodes,
                                    group_names[i],
                                    group_code[i])
                                    )
        }
      }

      # Display Gaussian mixture ellipses
      if(plot_communities && community_type == "Gaussian") {

        tikz_code = paste0(tikz_code, "% Gaussian ellipses for communities\n")
    
        n_groups = length(unique(self$communities))
        membership_vec = self$communities

        # Calculer centroïdes et ellipses par communauté
        for (comm in 1:n_groups) {

          # Items de cette communauté
          comm_items = which(membership_vec == comm)
          comm_coords = self$coordinates[comm_items, c(dim1,dim2)]
          
          if (nrow(comm_coords) == 1) {

            # Une seule variable dans le groupe
            centroid = comm_coords[1,]
            cov_matrix = diag(rep(0.1, ncol(comm_coords)))  # Ellipse minimale

          } else {

            # Centroïde (moyenne)
            centroid = colMeans(comm_coords)
            
            # Matrice de covariance
            cov_matrix = cov(comm_coords)

            # Assurer que la matrice est définie positive
            if (det(cov_matrix) <= 0) cov_matrix = cov_matrix + diag(rep(0.01, ncol(cov_matrix)))
          }
          
          # Paramètres de l'ellipse de confiance
          chi2_quantile = qchisq(confidence, df = 2)
          eigen_decomp = eigen(cov_matrix)
          eigenvals = pmax(eigen_decomp$values, 0.01)  # Éviter valeurs négatives
          eigenvecs = eigen_decomp$vectors
          
          # Demi-axes et angle de rotation
          semi_major = sqrt(eigenvals[1] * chi2_quantile)
          semi_minor = sqrt(eigenvals[2] * chi2_quantile)
          angle_rad = atan2(eigenvecs[2, 1], eigenvecs[1, 1])
          angle_deg = angle_rad * 180 / pi
          
          tikz_code = paste0(tikz_code,
                       sprintf("\\draw[community ellipse] (%.2f,%.2f) ellipse [x radius=%.2f, y radius=%.2f, rotate=%.1f];\n",
                              centroid[1],
                              centroid[2],
                              semi_major,
                              semi_minor,
                              angle_deg))
        }
      }

      # Display communities as (TikZ default) ellipses and draw meta-arrows
      if(plot_communities && community_type == "ellipses") {

        # Detect communities if not already done
        if(is.null(self$communities)) self$detect_graphCommunities()

        tikz_code = paste0(tikz_code, "% Gaussian ellipses for communities\n")
    
        n_groups = length(unique(self$communities))
        membership_vec = self$communities

        # Display community ellipses and meta-arrows
        for (comm in 1:n_groups) {

          # Community nodes
          comm_items = which(membership_vec == comm)
          comm_item_names = item_names[comm_items]
          comm_nodes = paste("(",comm_item_names,")",sep="",collapse=" ")

          # Use custom label if provided, otherwise default to "Cluster X"
          display_label = if(!is.null(cluster_labels) && comm <= length(cluster_labels)) cluster_labels[comm] else paste("Cluster", comm)

          tikz_code = paste0(tikz_code,
                            sprintf("\\node[auto-ellipse, fit=%s, scale=%.2f] (%s) {};\n",
                                    comm_nodes,
                                    dilate_cluster,
                                    paste("Cluster", comm, sep="")
                                    ))
          label_anchor = if(cluster_labels_below) "south" else "north"
          label_position = if(cluster_labels_below) "below" else "above"
          cluster_labels_code = paste0(cluster_labels_code,
                            sprintf("\\node[font=\\large\\bfseries, %s=3mm] at (%s.%s) {%s};\n",
                                    label_position,
                                    paste("Cluster", comm, sep=""),
                                    label_anchor,
                                    display_label))
        }

        # Display metagraph
        if(plot_metagraph) {
          tikz_code = paste0(tikz_code, "\n% Metagraph\n")
          for (i in 1:n_groups) {
            for (j in 1:n_groups) {
              if (self$metagraph[["adjacency"]][i,j] != 0 && i != j) {

                association_value = self$metagraph[["adjacency"]][i,j]
                from.node = paste("Cluster", i, sep="")
                to.node = paste("Cluster", j, sep="")

                # Draw arrow
                bend_direction = ifelse(i %in% meta.bend_right, "bend right", "bend left")
                      
                tikz_code = paste0(tikz_code,
                                  sprintf("\\draw[meta arrow, %s=%d] (%s) to[%s] (%s);\n",
                                          bend_direction,
                                          meta.bend_angle,
                                          from.node,
                                          bend_direction,
                                          to.node))

                # Draw arrow label (number of links)
                tikz_code = paste0(tikz_code,
                                  sprintf("\\path (%s) to[%s=%d] node[meta label] {%d} (%s);\n",
                                          from.node,
                                          bend_direction,
                                          meta.bend_angle,
                                          association_value,
                                          to.node))

              }
            }
          }
        }
      }

      # Display communities as convex hulls and optionally draw meta-arrows
      if(plot_communities && community_type == "convexhull") {

        # Detect communities if not already done
        if(is.null(self$communities)) self$detect_graphCommunities()

        tikz_code = paste0(tikz_code, "% Convex hulls for communities\n")

        n_groups = length(unique(self$communities))
        membership_vec = self$communities

        # Draw the community contours (always)
        for (comm in 1:n_groups) {

          # Community nodes
          comm_items = which(membership_vec == comm)
          comm_item_names = item_names[comm_items]
          coords = self$coordinates[comm_items, c(dim1,dim2)]

          # Compute convex hull
          hull_points = self$bounding_curve(coords, dilate_factor = dilate_cluster)

          # Find the topmost or bottommost point of the hull for label placement
          if (cluster_labels_below) {
            ref_idx = which.min(hull_points$y)
          } else {
            ref_idx = which.max(hull_points$y)
          }
          label_x = hull_points$x[ref_idx]
          label_y = hull_points$y[ref_idx]

          # Define and draw contours
          id = paste("Cluster", comm, sep="")
          coords_str = paste0("(", formatC(hull_points$x, digits=6, format="fg"),
                              ",", formatC(hull_points$y, digits=6, format="fg"), ")", collapse = " ")
          tikz_code = paste0(tikz_code,
                             sprintf("\n  %% Group %d\n", comm))
          # Create named path (for potential intersection calculations) AND draw it
          tikz_code = paste0(tikz_code,
                             sprintf("  \\draw[community ellipse, name path=%s] plot [smooth cycle, tension=0.6] coordinates {%s};\n",
                              id, coords_str))
          # Cluster label (deferred to main layer, printed last)
          display_label = if(!is.null(cluster_labels) && comm <= length(cluster_labels)) cluster_labels[comm] else paste("Cluster", comm)
          label_position = if(cluster_labels_below) "below" else "above"
          cluster_labels_code = paste0(cluster_labels_code,
                             sprintf("\\node[font=\\large\\bfseries, %s=3mm] at (%.6f,%.6f) {%s};\n",
                              label_position, label_x, label_y, display_label))

        }

        # Display metagraph arrows (optional)
        if(plot_metagraph) {

          if(is.null(self$metagraph[["adjacency"]])) self$create_metaGraph()

          tikz_code = paste0(tikz_code, "\n% Metagraph\n")

          # Group centers
          tikz_code = paste0(tikz_code, "\n% Group centers\n")
          for (i in 1:n_groups) {
            center_label = paste0("Center", i, sep="")
            tikz_code = paste0(tikz_code,
                                sprintf("\\coordinate (%s) at (%.6f,%.6f);\n",
                                center_label,
                                self$metagraph[["centroids"]][i,1],self$metagraph[["centroids"]][i,2]))
          }

          # Edges
          edge_id = 0L
          for (i in 1:n_groups) {
            for (j in 1:n_groups) {

              if (self$metagraph[["adjacency"]][i,j] != 0 && i != j) {

                edge_id = edge_id + 1L
                id1 = paste0("Cluster", i)
                id2 = paste0("Cluster", j)
                c1  = paste0("Center", i, sep="")
                c2  = paste0("Center", j, sep="")
                ln  = paste0("line", edge_id)
                pA  = paste0("pA", edge_id)
                pB  = paste0("pB", edge_id)

                association_value = self$metagraph[["adjacency"]][i,j]
                bend_direction = ifelse(i %in% meta.bend_right, "bend right", "bend left")

                # \path[name path=line12] (center1) -- (center2);
                # \path[name intersections={of=contour1 and line12, by=start12}];
                # \path[name intersections={of=contour2 and line12, by=end12}];
                # \draw[curved-arrow, blue!60, bend left=15] (start12) to (end12);

                # Define the line and compute the intersections with the group borders
                tikz_code = paste0(tikz_code,
                  sprintf("\n  %% Arrow %d: %s -> %s\n", edge_id, i, j),
                  sprintf("  \\path[name path=%s] (%s) -- (%s);\n", ln, c1, c2),
                  sprintf("  \\path[name intersections={of=%s and %s, by=%s}];\n", id1, ln, pA),
                  sprintf("  \\path[name intersections={of=%s and %s, by=%s}];\n", id2, ln, pB))

              }
            }
          }

        }
      }

      # Draw axes
      tikz_code = paste0(tikz_code, "\n% Axes\n")

      if(is.null(xlim)) xlim = c(min(c(self$coordinates[,dim1]),x.skills),max(c(self$coordinates[,dim1]),x.skills))
      if(is.null(ylim)) ylim = c(min(c(self$coordinates[,dim2]),y.skills),max(c(self$coordinates[,dim2]),y.skills))

      tikz_code = paste0(tikz_code,
                          sprintf("\\draw[axis style] (%.2f,0) -- (%.2f,0);\n",
                                  xlim[1],
                                  xlim[2]),
                          sprintf("\\draw[axis style] (0,%.2f) -- (0,%.2f);\n",
                                  ylim[1],
                                  ylim[2]))

      # Draw ticks
      tikz_code = paste0(tikz_code, "\n% Axis ticks\n")
      xticks = seq(-6,6,by=0.5)
      yticks = seq(-6,6,by=0.5)
      xticks = xticks[xticks > xlim[1] & xticks < xlim[2]]
      yticks = yticks[yticks > ylim[1] & yticks < ylim[2]]
      xtick.strings = paste(xticks,collapse=",")
      ytick.strings = paste(yticks,collapse=",")

      tikz_code = paste0(tikz_code,
                          sprintf("\\foreach \\x in {%s}\n  {\\draw (\\x,-%.2f) -- (\\x,%.2f); \\node[below] at (\\x,%.2f) {\\x};}\n",
                          xtick.strings,0.01,0.01,-0.01),
                          sprintf("\\foreach \\y in {%s}\n  {\\draw (-%.2f,\\y) -- (%.2f,\\y); \\node[left] at (%.2f,\\y) {\\y};}\n\n",
                          ytick.strings,0.01,0.01,-0.01))
      
      # Draw difficulty gradient line (behind everything)
      if(plot_difficulty_gradient) {
        C = self$coordinates[, c(dim1, dim2)]
        y = self$success.rates
        w = t(C) %*% y
        w = w / norm(w, "2")                     # unit direction vector
        center = colMeans(C)                      # centroid
        # Extended bounds (10% beyond ylim)
        ylim_ext = ylim + c(-0.1, 0.1) * diff(ylim)
        xlim_ext = xlim + c(-0.1, 0.1) * diff(xlim)
        # Find t values that reach extended boundaries
        t_candidates = c()
        if(w[1] != 0) t_candidates = c(t_candidates, (xlim_ext - center[1]) / w[1])
        if(w[2] != 0) t_candidates = c(t_candidates, (ylim_ext - center[2]) / w[2])
        # Keep only t values whose points fall within extended bounds
        pts = outer(t_candidates, as.vector(w)) + matrix(center, nrow=length(t_candidates), ncol=2, byrow=TRUE)
        inside = pts[,1] >= xlim_ext[1] - 0.01 & pts[,1] <= xlim_ext[2] + 0.01 &
                 pts[,2] >= ylim_ext[1] - 0.01 & pts[,2] <= ylim_ext[2] + 0.01
        t_inside = t_candidates[inside]
        t0 = min(t_inside); t1 = max(t_inside)
        # t0 = difficult end (away from w direction), t1 = easy end
        p_easy = center + t1 * as.vector(w)
        p_diff = center + t0 * as.vector(w)
        # Rotation angle for the label (in degrees)
        angle_deg = 180 + atan2(w[2], w[1]) * 180 / pi
        # Label position: 50% toward the difficult end
        p_label = center + 0.5 * t0 * as.vector(w)
        tikz_code = paste0(tikz_code, "\n% Difficulty gradient line\n")
        tikz_code = paste0(tikz_code,
                            sprintf("\\draw[gray!40, line width=2.5pt, dotted, -{Stealth[length=4mm, width=3mm]}] (%.4f,%.4f) -- (%.4f,%.4f);\n",
                                    p_easy[1], p_easy[2], p_diff[1], p_diff[2]))
        tikz_code = paste0(tikz_code,
                            sprintf("\\node[rotate=%.1f, font=\\large\\itshape, text=gray!40, above=1pt] at (%.4f,%.4f) {Difficulty gradient};\n",
                                    angle_deg, p_label[1], p_label[2]))
      }

      # Draw [-0.5;+0.5] interval band
      if(plot_band) {
        tikz_code = paste0(tikz_code, "\n% Interval band\n")
        tikz_code = paste0(tikz_code,
                            sprintf("\\draw[gray!80, dashed] (%.2f,-0.5) -- (%.2f,-0.5);\n",xlim[1],xlim[2]))
        tikz_code = paste0(tikz_code,
                            sprintf("\\draw[gray!80, dashed] (%.2f,0.5) -- (%.2f,0.5);\n",xlim[1],xlim[2]))
      }

      # Draw horizontal lines at detected logOR levels
      if(plot_logOR_levels && !is.null(self$logOR_levels)) {
        tikz_code = paste0(tikz_code, "\n% logOR level lines\n")
        for (lev in self$logOR_levels$centers) {
          tikz_code = paste0(tikz_code,
                              sprintf("\\draw[gray, dashed] (%.2f,%.4f) -- (%.2f,%.4f);\n",
                                      xlim[1], lev, xlim[2], lev))
        }
      }

      # Create bi-univoque links
      if(!only_directed) {

        if(is.null(from.list)) from.list = 2:n_items
        if(is.null(to.list)) to.list = 1:(n_items-1)

        tikz_code = paste0(tikz_code, "\n% Symmetric association links\n")
        for (i in from.list) {
          for (j in to.list) {
            # G.bi already encodes significance (P_OR_adj < alpha), no need for threshold filter
            # Only draw if both items are visible
            if ((self$G.bi[i, j] != 0) && (i %in% visible_items) && (j %in% visible_items)) {
              
              association_value = self$Intensity[i,j]
                          
              # Draw arrow
              tikz_code = paste0(tikz_code,
                                sprintf("\\draw[double arrow] (%s) to (%s);\n",
                                        item_names[i],
                                        item_names[j]))

              # draw label of the arc (implicative force)
              if(plot_coefs) {
                tikz_code = paste0(tikz_code,
                                  sprintf("\\path (%s) to node[weight label] {%.2f} (%s);\n",
                                          item_names[i],
                                          association_value,
                                          item_names[j]))
              }

            }

          }
        }
      }


      # Create univoque links (skip if metagraph is shown, as metagraph arrows summarize these)
      if(!only_symmetric && !plot_metagraph) {

        tikz_code = paste0(tikz_code, "\n% Univoque implication links\n")

        if(is.null(from.list)) from.list = 1:n_items
        if(is.null(to.list)) to.list = 1:n_items

        for (i in from.list) {
          for (j in to.list) {
            # G.uni already encodes significance (P_IOTA_adj < alpha), no need for threshold filter
            # Only draw if both items are visible
            if ((i != j) && (self$G.uni[i, j] != 0) && !(only_internals && self$communities[i] != self$communities[j]) && (i %in% visible_items) && (j %in% visible_items)) {
              
              opacity = opacity_univoque
              line_width = width_univoque

              association_value = self$Intensity[i,j]
              
              # Calcul de l'opacité et de l'épaisseur basées sur la force d'association
              # opacity = association_value
              # line_width = 0.2 + (association_value * 1.5)  # Entre 0.5 et 2
              
              # Déterminer la direction de courbure pour chaque source
              # Gauche par défaut, droit pour les items dans la liste bend_right
              bend_direction = ifelse(i %in% bend_right, "bend right", "bend left")
              
              # Couleur basée sur l'intensité
              # color_intensity = round(association_value * 100)
              # color_def = sprintf("black!%d", color_intensity)
              
              # Dessiner les arcs avec potentiellement des caractéristiques individuelles
              # tikz_code = paste0(tikz_code,
              #                   sprintf("\\draw[arrow, %s=%d, line width=%s, color=%s, opacity=%.2f] (%s) to[%s] (%s);\n",
              #                           bend_direction,
              #                           bend_angle,
              #                           width_univoque,
              #                           color_univoque,
              #                           opacity,
              #                           item_names[i],
              #                           bend_direction,
              #                           item_names[j]))

              tikz_code = paste0(tikz_code,
                                sprintf("\\draw[arrow, %s=%d] (%s) to[%s] (%s);\n",
                                        bend_direction,
                                        bend_angle,
                                        item_names[i],
                                        bend_direction,
                                        item_names[j]))

              # Dessiner le label de l'arc (force implicative)
              if(plot_coefs) {
                tikz_code = paste0(tikz_code,
                                  sprintf("\\path (%s) to[%s=%d] node[weight label] {%.2f} (%s);\n",
                                          item_names[i],
                                          bend_direction,
                                          bend_angle,
                                          association_value,
                                          item_names[j]))
              }

            }
          }
        }
      }

      # Generate the skill attributions
      if(!is.null(skills)) {
        for (i in 1:n_skills) {
          for (j in 1:n_items) {
            if(Q[j,i] != 0) {
              tikz_code = paste0(tikz_code,
                                sprintf("\\draw[-, dashed, line width=0.4, color=black!50, opacity=0.6] (%s) to (%s);\n",
                                        skill_code[i], item_names[j]))
            }
          }
        }
      }

      # Fin d'affichage sur le layer background
      tikz_code = paste0(tikz_code,"\\end{pgfonlayer}\n")

      # Cluster labels (printed last, on top of everything)
      if(nchar(cluster_labels_code) > 0) {
        tikz_code = paste0(tikz_code, "\n% Cluster labels (foreground)\n", cluster_labels_code)
      }

      # Fin du code TikZ
      tikz_code = paste0(tikz_code, "\n\\end{tikzpicture}")

      if(!is.null(legend)) {
        tikz_code = paste0(tikz_code, "\n% Legend\n")
        tikz_code = paste0(tikz_code,
                            sprintf("\\node[rectangle, draw=none, fill=black!5, rounded corners=8pt,anchor=north west, inner sep=12pt] at (%.2f,%.2f) {\n",x.legend,y.legend),
                            sprintf("\\begin{tabular}{l@{\\hspace{8mm}}l}
                              \\multicolumn{2}{c}{\\textbf{\\color{black!80}%s}} \\\\[3pt]
                              \\hline \\\\[-4pt]
                              %s
                            \\end{tabular}",legend.title, paste(legend,collapse="\\\\[3pt]\n")))
      }

      full_tex = paste0(
        "\\documentclass[border=1cm, tikz]{standalone}\n",
        "\\usepackage{tikz}\n",
        "\\usetikzlibrary{arrows,arrows.meta,fit,shapes,decorations.pathmorphing,positioning, calc, intersections}\n",
        "\\usepackage[utf8]{inputenc}\n",
        "\\usepackage[T1]{fontenc}\n\n",
        "\\pgfdeclarelayer{background}\n",
        "\\pgfsetlayers{background,main}\n",
        "\\begin{document}\n",
        tikz_code,
        "\n\\end{document}"
      )
      
      writeLines(full_tex, filename)
      cat("TikZ code saved as:", filename, "\n")
      cat("To compile: pdflatex", filename, "\n")
    },

    # Generate a LaTeX report of all 2x2 cross-tables
    crosstab_report = function(output_file = "crosstab_report.tex", title = "Binary Variable Cross-Tabulation Report") {
  
      data = self$mat

      if(!(class(data) %in% c("matrix", "data.frame"))) {
        stop("Data should be a matrix or a data.frame.")
      }

      # Input validation
      if (!is.data.frame(data)) {
        data = as.data.frame(data)
      }
      
      if (ncol(data) < 2) {
        stop("Data frame must have at least 2 variables for cross-tabulation")
      }
      
      # Get variable names and check binary nature
      vars = names(data)
      p = length(vars)
      
      # Validate binary variables
      for (var in vars) {
        unique_vals = sort(unique(data[[var]]))
        if (length(unique_vals) > 2) {
          stop(paste("Variable", var, "has more than 2 unique values:", 
                    paste(unique_vals, collapse = ", ")))
        }
        if (!all(unique_vals %in% c(0, 1, NA))) {
          warning(paste("Variable", var, "contains values other than 0 and 1:",
                      paste(unique_vals, collapse = ", ")))
        }
      }
      
      # Generate all ordered pairs (excluding self-pairs)
      pairs = expand.grid(var1 = vars, var2 = vars, stringsAsFactors = FALSE)
      pairs = pairs[pairs$var1 != pairs$var2, ]
      total_pairs = nrow(pairs)
      
      cat("Generating report for", p, "variables with", total_pairs, "ordered pairs...\n")
      
      # Initialize LaTeX document
      latex_content = c(
        "\\documentclass[10pt]{article}",
        "\\usepackage[margin=0.4in, top=0.6in, bottom=0.6in]{geometry}",
        "\\usepackage{array, booktabs, longtable}",
        "\\usepackage{amsmath, amsfonts}",
        "\\usepackage[table]{xcolor}",
        "\\renewcommand{\\arraystretch}{1.2}",
        "",
        "\\begin{document}",
        paste0("\\title{", title, "}"),
        "\\author{Generated by R Cross-Tabulation Function}",
        paste0("\\date{", Sys.Date(), "}"),
        "\\maketitle",
        "",
        paste0("\\section*{Summary}"),
        paste0("This report contains 2$\\times$2 cross-tabulations for all ordered pairs of the ", 
              p, " binary variables in the dataset. "),
        paste0("Total number of cross-tabulations: ", total_pairs, "."),
        "",
        "\\vspace{0.5cm}"
      )
      
      # Function to create a single formatted 2x2 table
      create_crosstab = function(var1_name, var2_name, data) {
        # Create contingency table
        ct = table(data[[var1_name]], data[[var2_name]], useNA = "no")
        
        # Ensure we have a proper 2x2 table with 1s and 0s as labels (1 first)
        full_ct = matrix(0, nrow = 2, ncol = 2, 
                          dimnames = list(c("1", "0"), c("1", "0")))
        
        # Fill in observed values
        row_names = rownames(ct)
        col_names = colnames(ct)
        
        for (i in seq_along(row_names)) {
          for (j in seq_along(col_names)) {
            r_idx = ifelse(row_names[i] == "1", 1, 2)
            c_idx = ifelse(col_names[j] == "1", 1, 2)
            full_ct[r_idx, c_idx] = ct[i, j]
          }
        }
        
        # Calculate marginal totals
        row_totals = rowSums(full_ct)
        col_totals = colSums(full_ct)
        grand_total = sum(full_ct)
        
        # Create LaTeX minipage with table
        table_latex = paste(c(
          "\\begin{minipage}[t]{0.19\\textwidth}",
          "\\centering",
          "\\footnotesize",
          paste0("\\textbf{", gsub("_", "\\_", var1_name, fixed = TRUE), 
                " $\\times$ ", gsub("_", "\\_", var2_name, fixed = TRUE), "}\\\\[0.1cm]"),
          "\\begin{tabular}{c|cc|c}",
          "\\toprule",
          paste0(" & \\multicolumn{2}{c|}{", gsub("_", "\\_", var2_name, fixed = TRUE), "} & \\\\"),
          "\\cmidrule{2-3}",
          paste0(gsub("_", "\\_", var1_name, fixed = TRUE), " & 1 & 0 & Total \\\\"),
          "\\midrule",
          paste0("1 & ", full_ct[1,1], " & ", full_ct[1,2], " & ", row_totals[1], " \\\\"),
          paste0("0 & ", full_ct[2,1], " & ", full_ct[2,2], " & ", row_totals[2], " \\\\"),
          "\\midrule",
          paste0("Total & ", col_totals[1], " & ", col_totals[2], " & ", grand_total, " \\\\"),
          "\\bottomrule",
          "\\end{tabular}",
          "\\end{minipage}"
        ), collapse = "\n")
        
        return(table_latex)
      }
      
      # Arrange tables with variable layout: 4 rows on first page, 6 rows on subsequent pages
      tables_per_row = 4
      first_page_rows = 4
      subsequent_page_rows = 6
      first_page_tables = first_page_rows * tables_per_row  # 20 tables
      subsequent_page_tables = subsequent_page_rows * tables_per_row  # 30 tables
      
      # Calculate page breaks
      table_idx = 1
      page_num = 1
      
      # Calculate total pages
      remaining_after_first = max(0, total_pairs - first_page_tables)
      total_pages = 1 + ceiling(remaining_after_first / subsequent_page_tables)
      
      while (table_idx <= total_pairs) {
        # Determine page parameters
        if (page_num == 1) {
          rows_this_page = first_page_rows
          tables_this_page = first_page_tables
        } else {
          rows_this_page = subsequent_page_rows
          tables_this_page = subsequent_page_tables
        }
        
        # Calculate actual tables to display on this page
        tables_remaining = total_pairs - table_idx + 1
        tables_to_show = min(tables_this_page, tables_remaining)
        
        # Add page break (except for first page)
        if (page_num > 1) {
          latex_content = c(latex_content, "\\newpage")
        }
        
        # Add page header
        latex_content = c(latex_content, 
                          paste0("\\subsection*{Cross-tabulations - Page ", page_num, 
                                  " of ", total_pages, "}"),
                          "\\vspace{0.3cm}")
        
        # Create tables row by row for this page
        for (row in 1:rows_this_page) {
          row_tables = c()
          
          for (col in 1:tables_per_row) {
            current_table_idx = table_idx + (row - 1) * tables_per_row + (col - 1)
            
            if (current_table_idx <= total_pairs && current_table_idx < table_idx + tables_to_show) {
              var1 = pairs$var1[current_table_idx]
              var2 = pairs$var2[current_table_idx]
              table_latex = create_crosstab(var1, var2, data)
              row_tables = c(row_tables, table_latex)
            } else {
              # Add empty minipage to maintain layout
              row_tables = c(row_tables, "\\begin{minipage}[t]{0.19\\textwidth} \\end{minipage}")
            }
          }
          
          # Add the row of tables
          if (length(row_tables) > 0) {
            latex_content = c(latex_content, 
                              paste(row_tables, collapse = "\\hfill"),
                              "\\\\[0.8cm]")
          }
        }
        
        # Move to next page
        table_idx = table_idx + tables_to_show
        page_num = page_num + 1
      }
      
      # Add summary statistics section
      latex_content = c(latex_content,
                        "\\newpage",
                        "\\section*{Dataset Summary}",
                        paste0("\\textbf{Number of variables:} ", p, "\\\\"),
                        paste0("\\textbf{Number of observations:} ", nrow(data), "\\\\"),
                        paste0("\\textbf{Total cross-tabulations:} ", total_pairs, "\\\\[0.5cm]"))
      
      # Add variable summary
      latex_content = c(latex_content, "\\textbf{Variable Summary:}\\\\[0.2cm]")
      
      for (var in vars) {
        var_summary = table(data[[var]], useNA = "ifany")
        prop_1 = ifelse("1" %in% names(var_summary), 
                        round(var_summary["1"] / sum(var_summary) * 100, 1), 0)
        
        latex_content = c(latex_content,
                          paste0("\\texttt{", gsub("_", "\\_", var, fixed = TRUE), 
                                  "}: ", prop_1, "\\% ones (", 
                                  ifelse("1" %in% names(var_summary), var_summary["1"], 0),
                                  "/", sum(var_summary), ")\\\\"))
      }
      
      # End document
      latex_content = c(latex_content, "", "\\end{document}")
      
      # Write to file
      writeLines(latex_content, output_file)
      
      # Return summary information
      remaining_after_first = max(0, total_pairs - first_page_tables)
      total_pages_final = 1 + ceiling(remaining_after_first / subsequent_page_tables)
      
      cat("✓ LaTeX report generated successfully!\n")
      cat("  File:", output_file, "\n")
      cat("  Variables:", p, "\n") 
      cat("  Observations:", nrow(data), "\n")
      cat("  Cross-tabulations:", total_pairs, "\n")
      cat("  Pages:", total_pages_final, "(", first_page_tables, "tables on page 1,", 
          subsequent_page_tables, "tables on subsequent pages)\n")
      cat("\nTo compile PDF:\n")
      cat("  pdflatex", output_file, "\n")
      
      # Return invisible list with metadata
      invisible(list(
        output_file = output_file,
        n_variables = p,
        n_observations = nrow(data),
        n_crosstabs = total_pairs,
        n_pages = total_pages_final,
        first_page_tables = first_page_tables,
        subsequent_page_tables = subsequent_page_tables
      ))
    }

  )
)
