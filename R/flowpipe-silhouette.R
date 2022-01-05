### Find intensity cutoffs for each channel.

#' @export
baseline_cut <- function(x, length_out = 100)
{
  breaks <- seq(min(x), max(x), length.out = length_out)
  h <- hist(x, breaks, plot = FALSE)
  baseline <- max(4, mean(h$counts) / 10)
  bound <- range(which(h$counts > baseline))
  r <- x[x >= breaks[bound[1]] & x < breaks[bound[2]]]

  r
}


#' @export
bisect <- function(
  x,
  max_sample = 2000,
  return_silhouette = FALSE,
  use_baseline = TRUE,
  valley_length = 100,
  min_x = 0,
  plot_cutoff = FALSE
)
{
  if (length(x) > max_sample) y <- sample(x, max_sample)
  else y <- x

  if (use_baseline) z <- baseline_cut(y)
  else z <- y

  valley <- seq(min(z, na.rm = TRUE), max(z, na.rm = TRUE), length.out = valley_length)
  if (!is.null(min_x))
    valley <- valley[valley > min_x]

  d <- stats::dist(z %>% `attributes<-`(NULL)) # 'dist()' has problems w/ some attributes

  sil <- sapply(valley,
    function(v)
    {
      cluster <- 1 * (z > v) + 1
      if (length(unique(cluster)) < 2)
        return (-2)

      ss <- cluster::silhouette(cluster, d)

      return (mean(ss[, 3]))
    }, simplify = TRUE)

  valley <- valley[which.max(sil)]

  if (plot_cutoff) {
    ## I'll need to jazz this up so I can check the flow data w/ file name & channel.
    plot(stats::density(z))
    plinth::vline(sprintf("%.2f", valley), abline... = list(col = "red"), text... = list(y = cp_coords()$y))
  }

  if (!return_silhouette)
    return (valley)
  else
    return (c(cutoff = valley, silhouette = max(sil)))
}

## usage:
# cutoff <- bisect(c(rnorm(1000, 0), rnorm(1000, 5, 10), rnorm(1000, 15, 5)), plot_cutoff = TRUE)


#' @export
trisect <- function(
  x,
  max_sample = 2000,
  use_baseline = TRUE,
  plot_cutoff = FALSE,
  ...
)
{
  if (length(x) > max_sample) y <- sample(x, max_sample)
  else y <- x

  if (use_baseline) z <- baseline_cut(y)
  else z <- y

  d <- stats::dist(z %>% `attributes<-`(NULL)) # 'dist()' has problems w/ some attributes

  cutoff1 <- bisect(z, max_sample = max_sample, use_baseline = use_baseline, ...)

  t1 <- bisect(z[z > cutoff1], max_sample = Inf, use_baseline = FALSE, min_x = NULL)
  cluster <- .bincode(z, sort(c(-Inf, cutoff1, t1, +Inf), decreasing = FALSE))
  ss <- cluster::silhouette(cluster, d)
  s1 <- mean(ss[, 3])

  t2 <- bisect(z[z <= cutoff1], max_sample = Inf, use_baseline = FALSE, min_x = NULL)
  cluster <- .bincode(z, sort(c(-Inf, cutoff1, t2, +Inf), decreasing = FALSE))
  ss <- cluster::silhouette(cluster, d)
  s2 <- mean(ss[, 3])

  if (s2 > s1)
    cutoff2 <- t2
  else
    cutoff2 <- t1

  r <- sort(c(cutoff1, cutoff2), decreasing = TRUE)

  if (plot_cutoff) {
    ## I'll need to jazz this up so I can check the flow data w/ file name & channel.
    plot(stats::density(z))
    plinth::vline(sprintf("%.2f", r), abline... = list(col = "red"), text... = list(y = cp_coords()$y))
  }

  return (r)
}

## usage:
# cutoffs <- trisect(c(rnorm(1000, 0), rnorm(1000, 5, 10), rnorm(1000, 15, 5)), plot_cutoff = TRUE)


#' @export
quadrisect <- function(
  x,
  max_sample = 2000,
  use_baseline = TRUE,
  plot_cutoff = FALSE,
  ...
)
{
  if (length(x) > max_sample) y <- sample(x, max_sample)
  else y <- x

  if (use_baseline) z <- baseline_cut(y)
  else z <- y

  d <- stats::dist(z %>% `attributes<-`(NULL)) # 'dist()' has problems w/ some attributes

  cutoff2 <- bisect(z, max_sample = max_sample, use_baseline = use_baseline, ...)

  z1 <- z[z <= cutoff2]; z2 <- z[z > cutoff2]
  cutoff1 <- bisect(z1, max_sample = max_sample, use_baseline = FALSE, min_x = NULL)
  cutoff3 <- bisect(z2, max_sample = max_sample, use_baseline = FALSE, min_x = NULL)

  r <- sort(c(cutoff1, cutoff2, cutoff3), decreasing = FALSE)

  if (plot_cutoff) {
    ## I'll need to jazz this up so I can check the flow data w/ file name & channel.
    plot(stats::density(z))
    plinth::vline(sprintf("%.2f", r), abline... = list(col = "red"), text... = list(y = cp_coords()$y))
  }

  return (r)
}

## usage:
# cutoffs <- quadrisect(c(rnorm(1000, 0), rnorm(1000, 5, 10), rnorm(1000, 15, 5)), plot_cutoff = TRUE)


#' @export
multisect <- function(
  x,
  bins = 2,
  max_sample = 2000,
  use_baseline = TRUE,
  plot_cutoff = FALSE,
  random_seed = NULL,
  resample_on_error = 5,
  ...
)
{
  bins <- head(bins, 1)
  if (bins %nin% 2:4)
    stop("Incorrect bin count")

  if (!is.null(resample_on_error)) {
    if (is.logical(resample_on_error)) {
      if (resample_on_error) resample_on_error <- 5
      else resample_on_error <- 1
    }
  } else {
    resample_on_error <- 1
  }

  ## Make sure there are enough random seeds for all possible failures.
  if (is.null(random_seed)) {
    random_seed <- 666 + seq(from = 0, length.out = resample_on_error + 1)
  } else {
    length(random_seed) <- resample_on_error + 1
    random_seed[is.na(random_seed)] <- max(random_seed, na.rm = TRUE) + seq(sum(is.na(random_seed)))
  }

  i <- 0
  run_multisect <- function()
  {
    while (TRUE) {
      withRestarts({
        if (i >= resample_on_error) {
          cat("\nError: Too many restarts"); flush.console()

          r <- rep(NA_real_, bins - 1) %>%
            `attr<-`("random_seed", random_seed[i])

          break
        }

        set.seed(random_seed[i + 1])

        if (length(x) > max_sample) y <- sample(x, max_sample)
        else y <- x

        if (use_baseline) z <- baseline_cut(y)
        else z <- y

        d <- stats::dist(z %>% `attributes<-`(NULL)) # 'dist()' has problems w/ some attributes

        cutoff1 <- bisect(z, max_sample = max_sample, use_baseline = use_baseline, ...)
        cutoff2 <- cutoff3 <- NULL

        if (bins != 2) {
          t1 <- bisect(z[z > cutoff1], max_sample = Inf, use_baseline = FALSE, min_x = NULL)
          if (bins == 3) {
            cluster <- .bincode(z, sort(c(-Inf, cutoff1, t1, +Inf), decreasing = FALSE))
            ss <- cluster::silhouette(cluster, d) # Original R code (now in C++): cluster:::silhouette.default.R
            s1 <- mean(ss[, 3])
          }

          t2 <- bisect(z[z <= cutoff1], max_sample = Inf, use_baseline = FALSE, min_x = NULL)
          if (bins == 3) {
            cluster <- .bincode(z, sort(c(-Inf, cutoff1, t2, +Inf), decreasing = FALSE))
            ss <- cluster::silhouette(cluster, d)
            s2 <- mean(ss[, 3])

            if (s2 > s1) cutoff2 <- t2
            else cutoff2 <- t1
          } else {
            cutoff2 <- t1; cutoff3 <- t2
          }
        }

        r <- sort(c(cutoff1, cutoff2, cutoff3), decreasing = FALSE)
        attr(r, "random_seed") <- random_seed[i + 1]

        if (plot_cutoff) {
          plot(stats::density(z), main = attr(x, "main_title"), cex.main = 0.8)
          plinth::vline(sprintf("%.2f", r), abline... = list(col = "red"), text... = list(y = plinth::cp_coords()$y))
        }

        break
      },
        restart = function() {
          use_baseline <<- !use_baseline

          i <<- i + 1
        })
    }

    r
  }


  tryCatch({
    withCallingHandlers({
        r <- run_multisect()
      },
        error = function(e) {
          if (any(grepl("must be a finite number", e$message, fixed = TRUE))) {
            cat("\n  Warning: Possibly insufficient positive events. Resampling..."); flush.console()
          } else {
            cat(sprintf("\n  Warning: %s. Resampling...", e$message)); flush.console()
          }

          invokeRestart("restart")
        }
    )
  }, error = function(e) { message("\nError: ", e$message); flush.console() })

  return (r)
}

## usage:
# cutoffs <- multisect(c(rnorm(1000, 0), rnorm(1000, 5, 10), rnorm(1000, 15, 5)), bins = 4, plot_cutoff = TRUE)
