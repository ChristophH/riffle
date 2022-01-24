#' @useDynLib riffle
NULL

#' Non-parametric differential expression test for sparse non-negative data
#'
#' @param y A matrix of counts; must be (or inherit from) class dgCMatrix; genes are row,
#' cells are columns
#' @param group_labels The group labels (e.g. cluster identities); 
#' will be converted to factor
#' @param compare Specifies which groups to compare, see details; default is 'each_vs_rest'
#' @param R The number of random permutations used to derive the p-values; default is 99
#' @param log2FC_th Threshold to remove genes from testing; absolute log2FC must be at least
#' this large for a gene to be tested; default is \code{log2(1.2)}
#' @param mean_th Threshold to remove genes from testing; gene mean must be at least this
#' large for a gene to be tested; default is 0.05
#' @param cells_th Threshold to remove genes from testing; gene must be detected (non-zero count)
#' in at least this many cells in the group with higher mean; default is 5
#' @param only_pos Test only genes with positive fold change (mean in group 1 > mean in group2); 
#' default is FALSE
#' @param only_top_n Test only the this number of genes from both ends of the log2FC spectrum
#' after all of the above filters have been applied; useful to get only the top markers; 
#' only used if set to a numeric value; default is NULL
#' @param mean_type Which type of mean to use; if \code{'geometric'} (default) the geometric mean is
#' used; to avoid \code{log(0)} we use \code{log1p} to add 1 to all counts and log-transform, 
#' calculate the arithmetic mean, and then back-transform and subtract 1 using \code{exp1m}; if
#' this parameter is set to \code{'arithmetic'} the data is used as is
#' @param verbosity Integer controlling how many messages the function prints; 
#' 0 is silent, 1 (default) is not
#'
#' @return Data frame of results
#' 
#' @section Details:
#' This model-free test is applied to each gene (row) individually but is
#' optimized to make use of the efficient sparse data representation of
#' the input. A permutation null distribution us used to assess the 
#' significance of the observed difference in mean between two groups.
#' 
#' The observed difference in mean is compared against a distribution
#' obtained by random shuffling of the group labels. For each gene every 
#' random permutation yields a difference in mean and from the population of
#' these background differences we estimate a mean and standard
#' deviation for the null distribution. 
#' This mean and standard deviation are used to turn the observed
#' difference in mean into a z-score and then into a p-value. Finally,
#' all p-values (for the tested genes) are adjusted using the Benjamini & Hochberg
#' method (fdr). The log2FC values in the output are \code{log2(mean1 / mean2)}.
#' Empirical p-values are also calculated: \code{emp_pval = (b + 1) / (R + 1)}
#' where b is the number of times the absolute difference in mean from a random 
#' permutation is at least as large as the absolute value of the observed difference
#' in mean, R is the number of random permutations. This is an upper bound of
#' the real empirical p-value that would be obtained by enumerating all possible
#' group label permutations.
#' 
#' There are multiple ways the group comparisons can be specified based on the compare
#' parameter. The default, \code{'each_vs_rest'}, does multiple comparisons, one per 
#' group vs all remaining cells. \code{'all_vs_all'}, also does multiple comparisons, 
#' covering all groups pairs. If compare is set to a length two character vector, e.g.
#' \code{c('T-cells', 'B-cells')}, one comparison between those two groups is done.
#' To put multiple groups on either side of a single comparison, use a list of length two. 
#' E.g. \code{compare = list(c('cluster1', 'cluster5'), c('cluster3'))}.
#' 
#' @import Matrix
#' @importFrom matrixStats rowMeans2 rowSds
#' @importFrom stats p.adjust pnorm
#' 
#' @export
#'
#' @examples
#' \donttest{
#' clustering <- 1:ncol(pbmc) %% 2
#' vst_out <- sctransform::vst(pbmc, return_corrected_umi = TRUE)
#' de_res <- diff_mean_test(y = vst_out$umi_corrected, group_labels = clustering)
#' }
#'
diff_mean_test <- function(y, group_labels, 
                           compare = 'each_vs_rest', 
                           R = 99, log2FC_th = log2(1.2), 
                           mean_th = 0.05, cells_th = 5, only_pos = FALSE,
                           only_top_n = NULL,
                           mean_type = 'geometric', 
                           verbosity = 1) {
  if (is.na(match(x = mean_type, table = c('geometric', 'arithmetic')))) {
    stop('mean_type must be geometric or arithmetic')
  }
  if (!inherits(x = y, what = 'dgCMatrix')) {
    stop('y must be a dgCMatrix')
  }
  if (R < 13) {
    stop('R must be at least 13')
  }
  if (!is.null(only_top_n) & (!is.numeric(only_top_n) | length(only_top_n) > 1)) {
    stop('only_top_n must be NULL or a single numeric value')
  }
  group_labels <- droplevels(as.factor(group_labels))
  lab_tab <- table(group_labels)
  group_levels <- levels(group_labels)
  G <- length(group_levels)
  if (length(group_labels) != ncol(y)) {
    stop('length of group labels must be equal to the number of columns in y')
  }
  
  if (verbosity > 0) {
    message('Non-parametric DE test for count data')
    message(sprintf('Using %s mean and %d random permutations', mean_type, R))
    message('Input: ', nrow(y), ' genes, ', ncol(y), ' cells; ', G, ' groups')
  }
  
  # Set up the comparisons we want to do; each comparison is a list
  # name1, name2, labels grp1, labels grp2
  if (compare[1] == 'each_vs_rest' && G == 2) {
    compare <- group_levels
    if (verbosity > 0) {
      message('There are only two groups in the data. Changing compare argument from "each_vs_rest" to group levels')
    }
  }
  if (compare[1] == 'each_vs_rest') {
    comparisons <- lapply(group_levels, function(x) list(x, 'rest', x, setdiff(group_levels, x)))
  } else if (compare[1] == 'all_vs_all') {
    comparisons <- list()
    for (i in 1:(G-1)) {
      for (j in (i+1):G) {
        comparisons[[length(comparisons) + 1]] <- list(group_levels[i], group_levels[j], group_levels[i], group_levels[j])
      }
    }
  } else if (inherits(x = compare, what = 'character') &&
             length(compare) == 2 &&
             all(compare %in% group_levels))  {
    if (compare[1] == compare[2]) {
      stop('Group 1 and 2 need to be different - please check your compare argument')
    }
    comparisons <- list(list(compare[1], compare[2], compare[1], compare[2]))
  } else if (inherits(x = compare, what = 'list') &&
             length(compare) == 2 &&
             all(unlist(lapply(compare, inherits, what = 'character')))) {
    compare <- lapply(compare, unique)
    if (length(intersect(compare[[1]], compare[[2]])) > 0) {
      stop('Intersection between group 1 and 2 - please check your compare argument')
    }
    comparisons <- list(list('group1', 'group2', compare[[1]], compare[[2]]))
  } else {
    stop("Make sure the compare argument is 'each_vs_rest' or 'all_vs_all' or a length 2 
          character vector with both entries present in the group_labels argument or 
          a list of length 2 with each entry being a character vector of group labels")
  }
  
  # for all the genes, get the number of non-zero observations per group
  cells <- row_nonzero_count_grouped_dgcmatrix(matrix = y, group = group_labels)
  # if we want to use the geometric mean, it's fastest to convert all counts to
  # log1p upfront, then use expm1 of arithmetic mean later on
  if (mean_type == 'geometric') {
    y@x <- log(y@x + 1)
    means <- row_mean_grouped_dgcmatrix(matrix = y, group = group_labels, shuffle = FALSE)
  } else {
    means <- row_mean_grouped_dgcmatrix(matrix = y, group = group_labels, shuffle = FALSE)
  }
  
  # Run the test for each comparison
  res_lst <- lapply(comparisons, function(comp) {
    # we might only be using a subset of the input cells; set up here
    sel_columns1 <- group_labels %in% comp[[3]]
    sel_columns2 <- group_labels %in% comp[[4]]
    sel_columns <- sel_columns1 | sel_columns2
    comp_group_labels <- factor(sel_columns2[sel_columns])
    
    if (verbosity > 0) {
      message(sprintf('Comparing %s (group1, N = %d) to %s (group2, N = %d)',
                      comp[[1]], sum(sel_columns1), comp[[2]], sum(sel_columns2)))
    }
    if (sum(sel_columns1) == 0 || sum(sel_columns2) == 0) {
      return()
    }
    comp_cells <- do.call(cbind, lapply(comp[3:4], function(x) combine_counts(cells, x)))
    comp_means <- do.call(cbind, lapply(comp[3:4], function(x) combine_means(means, lab_tab, x, mean_type)))
    
    res <- data.frame(gene = rownames(means),
                      group1 = comp[[1]],
                      mean1 = comp_means[, 1],
                      cells1 = comp_cells[, 1],
                      group2 = comp[[2]],
                      mean2 = comp_means[, 2],
                      cells2 = comp_cells[, 2])
    res$mean_diff <- res$mean1 - res$mean2
    res$log2FC <- log2(res$mean1 / res$mean2)
    
    # remove genes according to the filters
    if (log2FC_th > 0 || mean_th > 0 || cells_th > 0 || only_pos || !is.null(only_top_n)) {
      sel1 <- abs(res$log2FC) >= log2FC_th
      sel2 <- res$mean1 >= mean_th | res$mean2 >= mean_th
      sel3 <- (res$log2FC >= 0 & res$cells1 >= cells_th) | (res$log2FC <= 0 & res$cells2 >= cells_th)
      if (only_pos) {
        sel4 <- res$log2FC > 0
      } else {
        sel4 <- TRUE
      }
      res <- res[sel1 & sel2 & sel3 & sel4, , drop = FALSE]
      if (!is.null(only_top_n)) {
        sel0 <- rank(-res$log2FC) <= only_top_n
        if (!only_pos) {
          sel0 <- sel0 | rank(res$log2FC) <= only_top_n
        }
        res <- res[sel0, , drop = FALSE]
      }
      if (verbosity > 0) {
        message(sprintf('Keeping %d genes after initial filtering', nrow(res)))
      }
      # handle the case where no genes remain after filtering
      if (nrow(res) == 0) {
        return(res)
      }
    }
    
    # now get the empirical null distribution for mean_diff
    y_ss <- y[rownames(res), sel_columns, drop = FALSE]
    if (mean_type == 'geometric') {
      mean_diff_rnd <- do.call(cbind, lapply(1:R, function(i) {
        means_r <- expm1(row_mean_grouped_dgcmatrix(matrix = y_ss, group = comp_group_labels, shuffle = TRUE))
        means_r[, 1, drop = FALSE] - means_r[, 2, drop = FALSE]
      }))
    } else {
      mean_diff_rnd <- do.call(cbind, lapply(1:R, function(i) {
        means_r <- row_mean_grouped_dgcmatrix(matrix = y_ss, group = comp_group_labels, shuffle = TRUE)
        means_r[, 1, drop = FALSE] - means_r[, 2, drop = FALSE]
      }))
    }
    
    # use null distribution to get empirical p-values
    # also approximate null with normal and derive z-scores and p-values
    res$emp_pval <- (rowSums((abs(mean_diff_rnd) - abs(res$mean_diff)) >= 0) + 1) / (R + 1)
    res$emp_pval_adj <- p.adjust(res$emp_pval, method = 'BH')
    #res$zscore <- (res$mean_diff - rowMeans2(mean_diff_rnd)) / rowSds(mean_diff_rnd)
    sds <- sqrt(rowSums(mean_diff_rnd^2)/(R-1))
    res$zscore <- (res$mean_diff - rowMeans2(mean_diff_rnd)) / sds
    res$pval <- 2 * pnorm(-abs(res$zscore))
    res$pval_adj <- p.adjust(res$pval, method = 'BH')
    
    if (length(comparisons) > 1) {
      rownames(res) <- NULL
    }
    return(res)
  })
  res <- Reduce(rbind, res_lst)
  if (length(compare) == 1 && compare == 'each_vs_rest' && !is.null(res)) {
    res$group1 <- factor(res$group1, levels = group_levels)
    res$group2 <- factor(res$group2)
  } 
  if (length(compare) == 1 && compare == 'all_vs_all' && !is.null(res)) {
    res$group1 <- factor(res$group1, levels = group_levels)
    res$group2 <- factor(res$group2, levels = group_levels)
  }
  return(res)
}

# helper functions

combine_counts <- function(group_counts, columns) {
  as.matrix(rowSums(group_counts[, columns, drop = FALSE]))
}

# combine per-group-mean to get the mean spanning multiple groups
# in an act of irrational premature optimization, we pass the 
# log-space mean when mean_type is geometric - need to make sure to 
# transform with exp1m before returning
combine_means <- function(means, n_items, columns, mean_type) {
  if (length(columns) == 1) {
    if (mean_type == 'arithmetic') {
      return(means[, columns, drop = FALSE])
    }
    if (mean_type == 'geometric') {
      return(expm1(means[, columns, drop = FALSE]))
    }
  }
  means <- means[, columns]
  n_items <- n_items[columns]
  tmp <- sweep(x = means, MARGIN = 2, STATS = n_items, FUN = '*')
  
  if (mean_type == 'arithmetic') {
    return(as.matrix(rowSums(tmp) / sum(n_items)))
  }
  if (mean_type == 'geometric') {
    return(as.matrix(expm1(rowSums(tmp) / sum(n_items))))
  }
}

#' Find differentially expressed genes that are conserved across samples
#'
#' @param y A matrix of counts; must be (or inherit from) class dgCMatrix; genes are rows,
#' cells are columns
#' @param group_labels The group labels (i.e. clusters or time points); 
#' will be converted to factor
#' @param sample_labels The sample labels; will be converted to factor
#' @param balanced Boolean, see details for explanation; default is TRUE
#' @param compare Specifies which groups to compare, see details; currently only 'each_vs_rest' 
#' (the default) is supported
#' @param pval_th P-value threshold used to call a gene differentially expressed when summarizing 
#' the tests per gene
#' @param ... Parameters passed to diff_mean_test
#' 
#' @return Data frame of results
#' 
#' @section Details:
#' This function calls diff_mean_test repeatedly and aggregates the results per group and gene.
#' 
#' If balanced is TRUE (the default), it is assumed that each sample spans multiple groups, 
#' as would be the case when merging or integrating samples from the same tissue followed by 
#' clustering. Here the group labels would be the clusters and cluster markers would have support
#' in each sample.
#' 
#' If balanced is FALSE, an unbalanced design is assumed where each sample contributes to one
#' group. An example is a time series experiment where some samples are taken from time point 
#' 1 while other samples are taken from time point 2. The time point would be the group label
#' and the goal would be to identify differentially expressed genes between time points that
#' are supported by many between-sample comparisons.
#' 
#' Output columns:
#' \describe{
#' \item{group1}{Group label of the frist group of cells}
#' \item{group2}{Group label of the second group of cells; currently fixed to 'rest'}
#' \item{gene}{Gene name (from rownames of input matrix)}
#' \item{n_tests}{The number of tests this gene participated in for this group}
#' \item{log2FC_min,median,max}{Summary statistics for log2FC across the tests}
#' \item{mean1,2_median}{Median of group mean across the tests}
#' \item{pval_max}{Maximum of p-values across tests}
#' \item{de_tests}{Number of tests that showed this gene having a log2FC going in the same
#' direction as log2FC_median and having a p-value <= pval_th}
#' }
#' 
#' The output is ordered by group1, -de_tests, -abs(log2FC_median), pval_max
#' 
#' @import Matrix
#' @importFrom dplyr n group_by summarise arrange
#' @importFrom rlang .data
#' @importFrom stats median
#' @importFrom magrittr %>%
#' 
#' @export
#'
#' @examples
#' \donttest{
#' clustering <- 1:ncol(pbmc) %% 2
#' sample_id <- 1:ncol(pbmc) %% 3
#' vst_out <- sctransform::vst(pbmc, return_corrected_umi = TRUE)
#' de_res <- diff_mean_test_conserved(y = vst_out$umi_corrected, 
#' group_labels = clustering, sample_labels = sample_id)
#' }
#'
diff_mean_test_conserved <- function(y, group_labels, sample_labels, balanced = TRUE, 
                                     compare = 'each_vs_rest', pval_th = 1e-4, ...) {
  if (!inherits(x = y, what = 'dgCMatrix')) {
    stop('y must be a dgCMatrix')
  }
  group_labels <- droplevels(as.factor(group_labels))
  sample_labels <- droplevels(as.factor(sample_labels))
  
  res <- NULL
  if (compare[1] == 'each_vs_rest') {
    if (balanced) {
      res_lst <- lapply(levels(sample_labels), function(sl) {
        sel <- sample_labels == sl
        res <- diff_mean_test(y = y[, sel, drop = FALSE], group_labels = group_labels[sel], 
                              compare = compare, ...)
        if (!is.null(res)) {
          res$sample <- sl
        }
        res
      })
    } else {
      # fix special case when there are only two groups
      if (length(levels(group_labels)) == 2) {
        group_labels_to_do <- levels(group_labels)[1]
        gl_rest <- levels(group_labels)[2]
      } else {
        group_labels_to_do <- levels(group_labels)
        gl_rest <- 'rest'
      }
      # for each group, compare each sample against all samples that are not in group individually
      res_lst <- lapply(group_labels_to_do, function(gl) {
        gl_sel <- group_labels == gl
        samples_in_group <- sample_labels[gl_sel]
        res_lst <- lapply(unique(samples_in_group), function(sl_in_group) {
          sl_in_group_sel <- sample_labels == sl_in_group
          other_samples_not_in_group <- sample_labels[!(gl_sel | sl_in_group_sel)]
          res_lst <- lapply(unique(other_samples_not_in_group), function(osl_not_in_group) {
            sel <- (gl_sel & sl_in_group_sel) | (!gl_sel & sample_labels == osl_not_in_group)
            tmp_group <- c(gl, gl_rest)[as.numeric(!gl_sel & sample_labels == osl_not_in_group) + 1]
            res <- diff_mean_test(y = y[, sel, drop = FALSE], group_labels = tmp_group[sel], 
                                  compare = c(gl, gl_rest), ...)
            if (nrow(res) > 0) {
              res$sample1 <- sl_in_group
              res$sample2 <- osl_not_in_group
            } else {
              res <- NULL
            }
            res
          })
          do.call(rbind, res_lst)
        })
        do.call(rbind, res_lst)
      })
    }
    res <- do.call(rbind, res_lst)
    levels(res$group1) <- levels(group_labels)
  }
  if (!is.null(res)) {
    res <- group_by(res, .data$group1, .data$group2, .data$gene) %>%
      summarise(n_tests = n(),
                log2FC_min = min(.data$log2FC), 
                log2FC_median = median(.data$log2FC),
                log2FC_max = max(.data$log2FC),
                mean1_median = median(.data$mean1),
                mean2_median = median(.data$mean2),
                pval_max = max(.data$pval),
                de_tests = sum((sign(.data$log2FC) == sign(.data$log2FC_median)) &
                                 (.data$pval <= pval_th)),
                .groups = 'drop') %>%
      arrange(.data$group1, -.data$de_tests, -abs(.data$log2FC_median), .data$pval_max)
  }
  return(res)
}
