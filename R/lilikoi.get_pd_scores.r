lilikoi.get_pd_scores <- function(metaboliteMeasurements, metabolitePathwayTable, pathwayMetaboliteDatabase = lilikoi::data.smpdb) {
  #' A lilikoi.get_pd_scores Function
  #'
  #' This function allows you to compute Pathway Deregulation Score deriving
  #' make sure that you have the below database for the metabolites and pathway list:
  #' meta_path.RData
  #'
  #' @param metaboliteMeasurements Metabolite levels
  #' @param metabolitePathwayTable This is the metabolitePathwayTable from lilikoi.metab_to_pathway function. This table includes the metabolite ids and their corresponding HMDB ids
  #' @param pathwayMetaboliteDatabase An external database of metabolites and their associated pathways.
  #' @examples
  #' filename <- system.file("extdata", "plasma_breast_cancer.csv", package = "lilikoi")
  #' metaboliteMeasurements <- read.csv(file = filename, check.names = FALSE, row.names = 1)[, 1:20]
  #' metaboliteNames <- colnames(metaboliteMeasurements)[-1]
  #' matches <- lilikoi.metab_to_pathway(metaboliteNames, "name")
  #' PDSmatrix <- lilikoi.get_pd_scores(metaboliteMeasurements, matches, lilikoi::data.smpdb[1:23,])
  #' @export
  #' @import dplyr

  # We use filter_ for SE vs NSE so R CMD CHECk passes
  metaboliteIds <- metabolitePathwayTable %>% filter_("pathway != \"NA\"") %>% select("Query", "HMDB")

  # Keep just the columns which we have matched
  metaboliteLevelsWithMatchingQuery <- metaboliteMeasurements[, t(metaboliteIds["Query"])]

  # Now swap those names with the HMDB id
  colnames(metaboliteLevelsWithMatchingQuery) <- t(metaboliteIds["HMDB"])

  # Finally transpose that set
  metaboliteLevelsWithHmdbIds <- t(metaboliteLevelsWithMatchingQuery)

  # uses mapping dbs
  isHealthySample <- metaboliteMeasurements$Label %>% as.numeric %>% -1 %>% as.logical

  pathwayList <- .lilikoi.get_pathways_with_their_metabolites(pathwayMetaboliteDatabase)

  metaboliteNames <- row.names(metaboliteLevelsWithHmdbIds)

  ## although pathifier is built for gene inputs, turns out it works for metabolites as well todo:
  ## figure out why this is so slow
  pds <- .quantify_pathways_deregulation(as.matrix(metaboliteLevelsWithHmdbIds), allgenes = metaboliteNames,
    syms = pathwayList, pathwaynames = names(pathwayList), normals = isHealthySample, attempts = 5,
    min_exp = 0, min_std = 0)

  pdsMatrix <- matrix(as.data.frame(pds$scores), nrow = length(names(pds$scores)), byrow = TRUE)

  colnames(pdsMatrix) <- colnames(metaboliteLevelsWithHmdbIds)
  rownames(pdsMatrix) <- names(pds$scores)
  mode(pdsMatrix) <- "numeric"

  pdsMatrix
}

.lilikoi.get_pathways_with_their_metabolites <- function(pathwayMetaboliteTable) {
  # The below complicated line simply turns our metabolite pathway table into the list needed by
  # pathifier
  lapply(split(pathwayMetaboliteTable, pathwayMetaboliteTable[, "pathway"]), function(x) x[, setdiff(names(x),
    "pathway")])
}

.getmeasuredgenesinpathway <- function(syms, allgenes) {
  l <- length(syms)
  pathways <- vector("list", l)
  for (i in 1:l) {
    n <- length(syms[[i]])
    isin <- matrix(FALSE, n)
    for (j in 1:n) {
      isin[j] <- (length(grep(paste("\\b", trim(syms[[i]][j]), "\\b", sep = ""), allgenes)) > 0)
    }
    pathways[[i]] <- unique(syms[[i]][isin])
  }
  pathways
}

.getpathway <- function(sym, allgenes, data) {
  l <- length(sym)
  x <- NULL
  isin <- rep(FALSE, l)
  for (i in 1:l) {
    ind <- unique(grep(paste("\\b", trim(sym[i]), "\\b", sep = ""), allgenes))
    n <- length(ind)
    if (n > 0) {
      if (n == 1)
        t <- data[ind, ] else t <- colMeans(data[ind, ])
      if (var(t) > 0) {
        isin[i] <- TRUE
        x <- c(x, t)
      }
    }
  }
  if (is.null(x)) {
    list(x = NULL, isin = isin)
  } else {
    list(x = matrix(x, nrow = ncol(data)), isin = isin)
  }
}

#' @import princurve
.score_pathway <- function(x, m, ranks, calcerr = FALSE, thresh = 5e-04, maxit = 200, start, logfile = "") {
  x <- x[, apply(x, 2, stats::sd) > 0.001]
  k <- dim(x)[2]
  if (k < 3) {
    c <- NULL
    cat(file = logfile, append = TRUE, "scoring failed (k=", k, ").\n")
  } else {
    d <- matrix(0, 1, m)
    if (start == "by pca") {
      start <- NULL
    } else if (start == "by ranks") {
      start <- stats::aggregate(x, by = list(ranks), FUN = mean)
      start <- as.matrix(start[, -1])
    }
    c <- princurve::principal_curve(x, start = start, thresh = thresh, maxit = maxit)
  }
  if (!is.null(c)) {
    d[c$ord[1]] <- 0
    for (j in 2:m) {
      d[c$ord[j]] <- d[c$ord[j - 1]] + stats::dist(c$s[c$ord[(j - 1):j], ])
    }
    d <- d/d[c$ord[m]]
    if (calcerr) {
      e <- matrix(0, 1, k)
      for (i in 1:k) {
        e[i] <- mean((c$s[, i] - x[, i])^2)
      }
    } else {
      e <- FALSE
    }
    list(score = d, error = e, thecurve = c)
  } else {
    cat(file = logfile, append = TRUE, "scoring failed.\n")
    NULL
  }
}

.samplings_stdev <- function(m, n, attempts, z, ranks, samplings, start, logfile = "") {
  dall <- array(dim = c(attempts, n))
  skip <- 0
  for (a in 1:attempts) {
    res <- .score_pathway(z[samplings[a, ], ], m, ranks[samplings[a, ]], start = start, logfile = logfile)
    if (!is.null(res)) {
      dall[a, samplings[a, ]] <- res$score
    } else {
      skip <- skip + 1
    }
  }
  if (skip < attempts/2) {
    mean(apply(dall, 2, stats::sd, na.rm = TRUE), na.rm = TRUE)
  } else {
    Inf
  }

}

.score_all_pathways_helper <- function(z, ranks, samplings, i, attempts, maximize_stability, logfile = "",
  start) {
  n <- dim(z)[1]
  k <- dim(z)[2]
  m <- dim(samplings)[2]
  mincheck <- 5
  kmin <- max(floor(0.8 * k), mincheck + 1)
  mindelta <- min(0.009, max(0.002, 1.5/k))
  sig <- matrix(0, 1, k)
  res <- .score_pathway(z, n, ranks, calcerr = TRUE, start = start, logfile = logfile)
  if (is.null(res)) {
    cat(file = logfile, append = TRUE, "pathway ", i, "> scoring failed 1.\n")
  } else {
    sig <- .samplings_stdev(m, n, attempts, z, ranks, samplings, start = start)
    if (sig > 10000) {
      cat(file = logfile, append = TRUE, "pathway ", i, "> scoring failed 2 (sig:", sig, ").\n")
      res <- NULL
    } else {
      origsig <- sig
      cat(file = logfile, append = TRUE, "pathway ", i, "> sig:", sig, "\n")
      isin <- 1:k
      if (maximize_stability) {
        testsig <- max(mincheck, floor(0.1 * k))
        newsig <- rep(0, testsig)
        while ((k >= kmin) & (sig > 0.05)) {
          se <- sort(res$error, index.return = TRUE, decreasing = TRUE)
          for (j in 1:testsig) {
          newsig[j] <- .samplings_stdev(m, n, attempts, z[, -se$ix[j]], ranks, samplings, start = start)
          }
          wj <- which.min(newsig)
          cat(file = logfile, append = TRUE, "pathway ", i, " k=", k, "(", ncol(res$thecurve$s),
          ") wj=", wj, ">new sig:", newsig[wj])
          if (sig - newsig[wj] < mindelta) {
          cat(file = logfile, append = TRUE, " x rejected\n")
          break
          }
          cat(file = logfile, append = TRUE, " | accepted!\n")
          sig <- newsig[wj]
          isin <- isin[-se$ix[wj]]
          z <- z[, -se$ix[wj]]
          k <- k - 1
          res <- .score_pathway(z, n, ranks, calcerr = TRUE, start = start, logfile = logfile)
          if (is.null(res)) {
          cat(file = logfile, append = TRUE, "pathway ", i, "> scoring failed 3.\n")
          break
          }
        }
      }
    }
  }
  if (is.null(res)) {
    NULL
  } else {
    list(score = res$score, thecurve = res$thecurve, z = z, isin = isin, sig = sig, origsig = origsig,
      k = k)
  }
}

#' This function and it's private methods are ported from:
#' \url{https://www.bioconductor.org/packages/release/bioc/html/pathifier.html}
#'
#' We port it instead of including the bioconductor package as the bioconductor version
#' is not compatible with the latest princurve on cran.
#' @import R.oo
.quantify_pathways_deregulation <- function(data, allgenes, syms, pathwaynames, normals = NULL, ranks = NULL,
  attempts = 100, maximize_stability = TRUE, logfile = "", samplings = NULL, min_exp = 4, min_std = 0.4) {
  cat(file = logfile, append = FALSE, "robust_score_bydist. min_exp=", min_exp, ", min_std=", min_std,
    "\n")
  data[data < min_exp] <- min_exp
  n <- ncol(data)
  if (is.null(normals)) {
    normals <- rep(TRUE, n)
    start <- "by pca"
  } else {
    start <- "by ranks"
  }
  if (is.null(ranks))
    ranks <- !normals
  ranks <- rank(ranks)
  if ((length(normals) != n) || (length(ranks) != n)) {
    stop("invalid dimentions")
  }
  l <- length(syms)
  nn <- sum(normals)
  m <- floor(0.8 * (n - nn)) + nn
  if (is.null(samplings)) {
    samplings <- matrix(0, attempts, m)
    w <- which(!normals)
    for (a in 1:attempts) {
      samplings[a, ] <- sort(c(w[sample(n - nn, m - nn)], which(normals)))
    }
  }
  s <- NULL
  ind <- NULL
  for (i in 1:l) {
    # if using pathifier for large number of pathways, you might want to use the doMC library to
    # parallelize your code, in that case replace the above for with: s <- foreach (i=1:l,
    # .options.multicore=list(preschedule=FALSE), .combine=rbind, .inorder=FALSE, .verbose=FALSE,
    # .errorhandling='stop', .multicombine=TRUE) %dopar% {
    pathway <- syms[[i]]
    pathwayindata <- .getpathway(pathway, allgenes, data)
    k1 <- sum(pathwayindata$isin)
    if (k1 < 3) {
      si <- NULL
      cat(file = logfile, append = TRUE, "skipping pathway ", i, " k1=", k1, "\n")
    } else {
      x <- pathwayindata$x
      pathway <- pathway[pathwayindata$isin]
      xm <- colMeans(x[normals, ])
      xs <- apply(x[normals, ], 2, stats::sd)
      xs[xs < min_std] <- min_std
      if (0 %in% xs) {
        si <- NULL
        cat(file = logfile, append = TRUE, "skipping pathway ", i, " (0 in xs)\n")
      } else {
        z <- (x - matrix(rep(xm, each = n), nrow = n))/(matrix(rep(xs, each = n), nrow = n))
        t <- stats::prcomp(z)
        k2 <- max(sum(t$sdev > 1.1), 4)
        k2 <- min(k2, k1, 0.75 * dim(x)[1], sum(t$sdev > 0.25))
        if (k2 < 3) {
          si <- NULL
          cat(file = logfile, append = TRUE, "skipping pathway ", i, " k2=", k2, "\n")
        } else {
          pca <- t$x[, 1:k2]
          res <- .score_all_pathways_helper(pca, ranks, samplings, i, attempts, maximize_stability,
          logfile, start = start)
          if (is.null(res)) {
          si <- NULL
          cat(file = logfile, append = TRUE, "skipping pathway ", i, "\n")
          } else {
          ind <- c(ind, i)
          si <- list(res$score, pathway, res$sig, res$origsig, res$k, res$thecurve$s, res$thecurve$ord,
            res$z, res$isin, xm, xs, t$center, t$rotation, k2)
          }
        }
      }
    }
    s <- rbind(s, si)  # if using doMC foreach above replace this line with simply: si
  }
  cat(file = logfile, append = TRUE, length(ind), "pathways processed with start=", start, "\n")
  rownames(s) <- pathwaynames[ind]
  list(scores = s[, 1], genesinpathway = s[, 2], newmeanstd = s[, 3], origmeanstd = s[, 4], pathwaysize = s[,
    5], curves = s[, 6], curves_order = s[, 7], z = s[, 8], compin = s[, 9], xm = s[, 10], xs = s[,
    11], center = s[, 12], rot = s[, 13], pctaken = s[, 14], samplings = samplings, sucess = ind,
    logfile = logfile)
}


