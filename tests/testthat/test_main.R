context("diff_mean_test and diff_mean_test_conserved")

test_that('diff_mean_test works', {
  skip_on_cran()
  set.seed(42)
  
  grouping <- as.factor(sample(c('a','b','c'), size = ncol(riffle::pbmc), replace = TRUE))
  
  # default 
  de_res <- riffle::diff_mean_test(y = riffle::pbmc, group_labels = grouping)
  
  # one cell vs many
  sel <- c(which(grouping == 'a')[1], which(grouping == 'b'))
  de_res <- riffle::diff_mean_test(y = riffle::pbmc[, sel], group_labels = grouping[sel])
  
  # one cell vs few
  sel <- c(which(grouping == 'a')[1], which(grouping == 'b')[c(1, 2)])
  de_res <- riffle::diff_mean_test(y = riffle::pbmc[, sel], group_labels = grouping[sel])
  
  # one cell vs one
  sel <- c(which(grouping == 'a')[1], which(grouping == 'b')[1])
  de_res <- riffle::diff_mean_test(y = riffle::pbmc[, sel], group_labels = grouping[sel])
})

test_that('diff_mean_test_conserved unbalanced works', {
  skip_on_cran()
  set.seed(42)
  
  samples <- as.factor(sample(c('a1', 'a2','b1', 'b2'), size = ncol(riffle::pbmc), replace = TRUE))
  condition <- as.factor(gsub(pattern = '\\d$', replacement = '', x = samples))
  
  # default scenario
  # print(table(samples, condition))
  de_res <- riffle::diff_mean_test_conserved(y = riffle::pbmc, 
                                             group_labels = condition, 
                                             sample_labels = samples, 
                                             balanced = FALSE)
  
  # one sample has just one cell
  sel <- rep(x = TRUE, times = ncol(riffle::pbmc))
  sel[which(samples == 'a2')] <- FALSE
  sel[which(samples == 'a2')[1]] <- TRUE
  # print(table(samples[sel], condition[sel]))
  de_res <- riffle::diff_mean_test_conserved(y = riffle::pbmc[, sel], 
                                             group_labels = condition[sel], 
                                             sample_labels = samples[sel], 
                                             balanced = FALSE,
                                             R = 499,
                                             log2FC_th = 0, 
                                             mean_th = 0.01)
  
  # two samples have just one cell
  sel <- rep(x = TRUE, times = ncol(riffle::pbmc))
  sel[which(samples == 'a1')] <- FALSE
  sel[which(samples == 'a1')[1]] <- TRUE
  sel[which(samples == 'a2')] <- FALSE
  sel[which(samples == 'a2')[1]] <- TRUE
  # print(table(samples[sel], condition[sel]))
  de_res <- riffle::diff_mean_test_conserved(y = riffle::pbmc[, sel], 
                                             group_labels = condition[sel], 
                                             sample_labels = samples[sel], 
                                             balanced = FALSE,
                                             R = 499,
                                             log2FC_th = 0, 
                                             mean_th = 0.01)
  
  # three samples have just one cell
  sel <- rep(x = TRUE, times = ncol(riffle::pbmc))
  sel[which(samples == 'a1')] <- FALSE
  sel[which(samples == 'a1')[1]] <- TRUE
  sel[which(samples == 'a2')] <- FALSE
  sel[which(samples == 'a2')[1]] <- TRUE
  sel[which(samples == 'b1')] <- FALSE
  sel[which(samples == 'b1')[1]] <- TRUE
  # print(table(samples[sel], condition[sel]))
  de_res <- riffle::diff_mean_test_conserved(y = riffle::pbmc[, sel], 
                                             group_labels = condition[sel], 
                                             sample_labels = samples[sel], 
                                             balanced = FALSE,
                                             R = 499,
                                             log2FC_th = 0, 
                                             mean_th = 0.01)
  
  
})

