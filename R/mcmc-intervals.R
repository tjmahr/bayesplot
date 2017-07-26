#' Plot interval estimates from MCMC draws
#'
#' Plot central (quantile-based) interval estimates from MCMC draws. See the
#' \strong{Plot Descriptions} section, below, for details.
#'
#' @name MCMC-intervals
#' @family MCMC
#'
#' @template args-mcmc-x
#' @template args-pars
#' @template args-regex_pars
#' @template args-transformations
#' @param ... Currently unused.
#' @param prob The probability mass to include in the inner interval (for
#'   \code{mcmc_intervals}) or in the shaded region (for \code{mcmc_areas},
#'   \code{mcmc_joy}). The default is \code{0.5} for \code{mcmc_intervals} and
#'   \code{mcmc_areas} and \code{0.8} for \code{mcmc_joy}.
#' @param prob_outer The probability mass to include in the outer interval. The
#'   default is \code{0.9} for \code{mcmc_intervals}
#'   and \code{1} for \code{mcmc_areas} and \code{mcmc_joy}.
#' @param point_est The point estimate to show. Either \code{"median"} (the
#'   default), \code{"mean"}, or \code{"none"}.
#' @param rhat An optional numeric vector of \eqn{\hat{R}}{Rhat} estimates, with
#'   one element per parameter included in \code{x}. If \code{rhat} is provided,
#'   the intervals/areas and point estimates in the resulting plot are colored
#'   based on \eqn{\hat{R}}{Rhat} value. See \code{\link{rhat}} for methods for
#'   extracting \eqn{\hat{R}}{Rhat} estimates.
#' @param bw,adjust,kernel For \code{mcmc_areas} and \code{mcmc_joy}, optional
#'   arguments passed to \code{\link[stats]{density}} to override default kernel
#'   density estimation parameters.
#'
#' @template return-ggplot
#'
#' @section Plot Descriptions:
#' \describe{
#'   \item{\code{mcmc_intervals}}{
#'    Plots of uncertainty intervals computed from posterior draws with all
#'    chains merged.
#'   }
#'   \item{\code{mcmc_areas}}{
#'    Density plots computed from posterior draws with all chains merged,
#'    with uncertainty intervals shown as shaded areas under the curves.
#'   }
#'   \item{\code{mcmc_joy}}{
#'    Similar to \code{mcmc_areas} but the density plots are staggered and
#'    ordered by posterior median. (See the \pkg{ggjoy} package for more
#'    details.) The implementation of \code{mcmc_joy} in \pkg{bayesplot}
#'    was inspired by a blog post by TJ Mahr
#'    (\url{http://rpubs.com/tjmahr/joyplot}).
#'   }
#' }
#'
#' @examples
#' # some parameter draws to use for demonstration
#' x <- example_mcmc_draws(params = 6)
#' dim(x)
#' dimnames(x)
#'
#' color_scheme_set("brightblue")
#' mcmc_intervals(x)
#' mcmc_intervals(x, pars = c("beta[1]", "beta[2]"))
#' mcmc_areas(x, regex_pars = "beta\\[[1-3]", prob = 0.8) +
#'  ggplot2::labs(
#'    title = "Posterior distributions",
#'    subtitle = "with medians and 80% intervals"
#'  )
#'
#' color_scheme_set("red")
#' mcmc_areas(
#'    x,
#'    pars = c("alpha", "beta[4]"),
#'    prob = 2/3,
#'    prob_outer = 0.9,
#'    point_est = "mean"
#' )
#'
#' # color by rhat value
#' color_scheme_set("blue")
#' fake_rhat_values <- c(1, 1.07, 1.3, 1.01, 1.15, 1.005)
#' mcmc_intervals(x, rhat = fake_rhat_values)
#'
#' color_scheme_set("gray")
#' p <- mcmc_areas(x, pars = c("alpha", "beta[4]"), rhat = c(1, 1.1))
#' p + legend_move("bottom")
#' p + legend_move("none") # or p + legend_none()
#'
#' \donttest{
#' # apply transformations
#' mcmc_intervals(
#'   x,
#'   pars = c("beta[2]", "sigma"),
#'   transformations = list("sigma" = "log", "beta[2]" = function(x) x + 3)
#' )
#'
#' # apply same transformation to all selected parameters
#' mcmc_intervals(x, regex_pars = "beta", transformations = "exp")
#' }
#'
#' \dontrun{
#' # example using fitted model from rstanarm package
#' library(rstanarm)
#' fit <- stan_glm(
#'  mpg ~ 0 + wt + factor(cyl),
#'  data = mtcars,
#'  iter = 500
#' )
#' x <- as.matrix(fit)
#'
#' color_scheme_set("teal")
#' mcmc_intervals(x, point_est = "mean", prob = 0.8, prob_outer = 0.95)
#' mcmc_areas(x, regex_pars = "cyl", bw = "SJ",
#'            rhat = rhat(fit, regex_pars = "cyl"))
#' }
#'
#'
NULL

#' @rdname MCMC-intervals
#' @export
mcmc_intervals <- function(x,
                           pars = character(),
                           regex_pars = character(),
                           transformations = list(),
                           ...,
                           prob = 0.5,
                           prob_outer = 0.9,
                           point_est = c("median", "mean", "none"),
                           rhat = numeric()) {
  check_ignored_arguments(...)
  stopifnot(prob_outer >= prob)
  x <- prepare_mcmc_array(x, pars, regex_pars, transformations)
  .mcmc_intervals(
    x = merge_chains(x),
    prob_inner = prob,
    prob_outer = prob_outer,
    point_est = point_est,
    show_density = FALSE,
    rhat = rhat
  )
}

#' @rdname MCMC-intervals
#' @export
mcmc_areas <- function(x,
                       pars = character(),
                       regex_pars = character(),
                       transformations = list(),
                       ...,
                       prob = 0.5,
                       prob_outer = 1,
                       point_est = c("median", "mean", "none"),
                       rhat = numeric(),
                       bw = NULL,
                       adjust = NULL,
                       kernel = NULL) {
  check_ignored_arguments(...)
  stopifnot(prob_outer >= prob)
  x <- prepare_mcmc_array(x, pars, regex_pars, transformations)
  .mcmc_intervals(
    x = merge_chains(x),
    prob_inner = prob,
    prob_outer = prob_outer,
    point_est = point_est,
    rhat = rhat,
    show_density = TRUE,
    bw = bw,
    adjust = adjust,
    kernel = kernel
  )
}

#' @rdname MCMC-intervals
#' @export
#' @param height_scalar For \code{mcmc_joy}, a multiplicative factor used
#'   to increase or decrease the heights of the density plots.
#' @param order_pars For \code{mcmc_joy}, by default the parameters are plotted
#'   in order of their median value. To override this default and plot them in
#'   the order provided in \code{x} set \code{order_pars} to \code{FALSE}.
#'
#' @importFrom dplyr rename_ do_ left_join right_join
mcmc_joy <- function(x,
                     pars = character(),
                     regex_pars = character(),
                     transformations = list(),
                     ...,
                     prob = 0.8,
                     prob_outer = 1,
                     bw = "nrd0",
                     adjust = 1,
                     kernel = "gaussian",
                     height_scalar = 1.5,
                     order_pars = TRUE) {
  check_ignored_arguments(...)
  stopifnot(prob_outer >= prob, length(height_scalar) == 1)
  x <- prepare_mcmc_array(x, pars, regex_pars, transformations)
  x <- merge_chains(x)

  plotdata <- .mcmc_joy_melt_and_order(x, order_pars)
  dens_inner <- .mcmc_joy_dens_data(plotdata, prob, bw, adjust, kernel)
  dens_outer <- .mcmc_joy_dens_data(plotdata, prob_outer, bw, adjust, kernel)

  ggplot(dens_outer) +
    aes_(x = ~ value, y = ~ Parameter, height = ~ height * height_scalar) +
    hline_at(
      1:length(unique(dens_outer$Parameter)),
      size = 0.1,
      color = "gray"
    ) +
    ggjoy::geom_ridgeline(
      data = dens_inner,
      color = NA,
      fill = get_color("l")
    ) +
    ggjoy::geom_ridgeline(
      fill = NA,
      color = get_color("dh")
    ) +
    scale_y_discrete(expand = c(0.01, 0)) +
    dont_expand_x_axis() +
    yaxis_text(face = "bold") +
    yaxis_title(FALSE) +
    xaxis_title(FALSE)
}



# internal ----------------------------------------------------------------

# @param x A matrix (not a 3-D array) created by merge_chains()
# @param order_pars User's order_pars argument
# @return A data frame with columns 'Parameter' (factor) and 'value' (numeric).
#   If order_pars=TRUE then then the factor levels of 'Parameter' are arranged
#   by the median of the 'value' variable.
#
.mcmc_joy_melt_and_order <- function(x, order_pars = FALSE) {
  mx <- reshape2::melt(x)[, -1]
  colnames(mx) <- c("Parameter", "value")

  if (!order_pars)
    return(mx)

  mx %>%
    group_by_( ~ Parameter) %>%
    summarise_(Med = ~ median(value)) %>%
    ungroup() %>%
    mutate_(OrdParameter = ~ factor(Parameter, levels = Parameter[order(Med)])) %>%
    select_(.dots = list( ~ Parameter, ~ OrdParameter)) %>%
    dplyr::right_join(mx, by = "Parameter") %>%
    select_(.dots = list( ~ -Parameter)) %>%
    rename_(Parameter = ~ OrdParameter)
}

#
# @param x The data frame returned by .mcmc_joy_data
# @param prob The user's 'prob' or 'prob_outer' argument
# @param bw,adjust,kernel Optional args passed to stats::density
#
# @return A data frame with columns 'Parameter' (factor), (factor), 'value'
#   (numeric), 'height' (numeric) that can be passed to ggjoy::geom_ridgeline
.mcmc_joy_dens_data <- function(x, prob, bw, adjust, kernel) {
  x %>%
    group_by_(.dots = list(~ Parameter)) %>%
    do_(.dots = list(~ .compute_dens_i_joy(
      .$value,
      prob = prob,
      bw = bw,
      adjust = adjust,
      kernel = kernel
    ))) %>%
    ungroup() %>%
    mutate_(height = ~ dens / max(dens)) %>%
    select_(.dots = list(~ Parameter, ~ height, ~ value))
}

# @param x A matrix (not a 3-D array) created by merge_chains()
.mcmc_intervals <- function(x,
                            prob_inner = 0.5,
                            prob_outer = 0.95,
                            point_est = c("median", "mean", "none"),
                            rhat = numeric(),
                            show_density = FALSE,
                            bw = NULL,
                            adjust = NULL,
                            kernel = NULL) {
  n_param <- ncol(x)
  parnames <- colnames(x)

  probs <- c(0.5 - prob_outer / 2,
             0.5 - prob_inner / 2,
             0.5,
             0.5 + prob_inner / 2,
             0.5 + prob_outer / 2)

  quantiles <- t(apply(x, 2, quantile, probs = probs))
  y <- as.numeric(seq(n_param, 1, by = -1))
  x_lim <- c(min(quantiles[, 1]), max(quantiles[, 5]))
  x_range <- diff(x_lim)
  x_lim[1] <- x_lim[1] - 0.05 * x_range
  x_lim[2] <- x_lim[2] + 0.05 * x_range

  data <- data.frame(parnames, y, quantiles)
  colnames(data) <- c("parameter", "y", "ll", "l", "m", "h", "hh")
  point_est <- match.arg(point_est)
  no_point_est <- identical(point_est, "none")
  if (point_est == "mean")
    data$m <- unname(colMeans(x))

  color_by_rhat <- isTRUE(length(rhat) > 0)
  if (color_by_rhat) {
    rhat <- factor(factor_rhat(rhat), levels = c("high", "ok", "low"))
    if (length(rhat) != nrow(data))
      stop(
        "'rhat' has length ", length(rhat),
        " but 'x' has ", nrow(data), " parameters.",
        call. = FALSE
      )

    data$rhat <- rhat
  }

  graph <- ggplot(data)

  # faint vertical line at zero if zero is within x_lim
  if (0 > x_lim[1] && 0 < x_lim[2])
    graph <- graph + vline_0(color = "gray90", size = 0.5)

  if (show_density) {
    # density outline
    n_dens_pts <- 512
    y_dens <- matrix(0, nrow = n_dens_pts, ncol = n_param)
    x_dens <- matrix(0, nrow = n_dens_pts, ncol = n_param)
    for (i in seq_len(n_param)) {
      d_temp <- .compute_dens_i(
        x = x[, i],
        from = quantiles[i, 1],
        to = quantiles[i, 5],
        n = n_dens_pts,
        bw = bw,
        adjust = adjust,
        kernel = kernel
      )
      x_dens[, i] <- d_temp$x
      y_max <- max(d_temp$y)
      y_dens[, i] <- d_temp$y / y_max * 0.8 + y[i]
    }
    df_dens <- data.frame(
      x = as.vector(x_dens),
      y = as.vector(y_dens),
      name = rep(parnames, each = n_dens_pts)
    )
    if (color_by_rhat)
      df_dens$rhat <- rep(rhat, each = n_dens_pts)

    dens_args <- list(
      data = df_dens,
      mapping = aes_(
        x = ~ x,
        y = ~ y,
        group = ~ name,
        color = if (!color_by_rhat) NULL else ~ rhat
      ),
      lineend = "round"
    )
    if (!color_by_rhat)
      dens_args$color <- get_color("d")
    g_dens <- do.call("geom_line", dens_args)

    #shaded interval
    y_poly <- matrix(0, nrow = n_dens_pts + 2, ncol = n_param)
    x_poly <- matrix(0, nrow = n_dens_pts + 2, ncol = n_param)
    for (i in seq_len(n_param)) {
      d_temp <- .compute_dens_i(
        x = x[, i],
        from = quantiles[i, 2],
        to = quantiles[i, 4],
        n = n_dens_pts,
        bw = bw,
        adjust = adjust,
        kernel = kernel
      )
      x_poly[, i] <-
        c(d_temp$x[1], as.vector(d_temp$x), d_temp$x[n_dens_pts])
      y_max <- max(d_temp$y)
      y_poly[, i] <-
        as.vector(c(0, as.vector(d_temp$y) / y_max * 0.8, 0) + y[i])
    }
    df_poly <-
      data.frame(
        x = as.vector(x_poly),
        y = as.vector(y_poly),
        name = rep(parnames, each = n_dens_pts + 2)
      )
    if (color_by_rhat)
      df_poly$rhat <- rep(rhat, each = n_dens_pts + 2)
    g_poly <-
      geom_polygon(data = df_poly, aes_(
        x = ~ x,
        y = ~ y,
        group = ~ name,
        fill = if (color_by_rhat) ~ rhat else ~ y
      ))

    # point estimate
    df_dens$parameter <- df_dens$name
    pt_data <- dplyr::summarise_(
      # find y value at which to stop vertical pt est segment
      dplyr::group_by_(
        dplyr::left_join(df_dens, data[, c("parameter", "m")],
                         by = "parameter"),
        .dots = list(~ parameter)
      ),
      .dots = list(maxy = ~ y[which.min(abs(x - m))])
    )
    segment_args <- list(
      data = dplyr::left_join(data, pt_data, by = "parameter"),
      mapping = aes_(
        x = ~ m,
        xend = ~ m,
        y = ~ y,
        yend = ~ maxy,
        color = if (!color_by_rhat) NULL else ~ rhat
      ),
      size = 1.25
    )
    if (!color_by_rhat)
      segment_args$color <- get_color("m")
    g_point <- do.call("geom_segment", segment_args)

    # bottom line
    bottom_args <- list(
      mapping = aes_(
        x = ~ ll,
        xend = ~ hh,
        y = ~ y,
        yend = ~ y,
        color = if (!color_by_rhat) NULL else ~ rhat
      )
    )
    if (!color_by_rhat)
      bottom_args$color <- get_color("d")
    g_bottom <- do.call("geom_segment", bottom_args)

    graph <- graph + g_poly
    if (!no_point_est)
      graph <- graph + g_point
    graph <- graph + g_bottom + g_dens

    if (color_by_rhat) {
      graph <- graph +
        scale_fill_diagnostic("rhat") +
        scale_color_diagnostic("rhat")
    } else {
      graph <- graph + scale_fill_gradient(low = get_color("l"),
                                           high = get_color("l"),
                                           guide = "none")
    }

  } else { # No densities

    # outer interval
    graph <-
      graph + geom_segment(aes_(
        x = ~ ll,
        xend = ~ hh,
        y = ~ y,
        yend = ~ y
      ),
      colour = get_color("m"),
      lineend = "round")

    # inner interval
    segment_args <- list(
      mapping = aes_(
        x = ~ l,
        xend = ~ h,
        y = ~ y,
        yend = ~ y,
        color = if (!color_by_rhat) NULL else ~ rhat
      ),
      size = 2,
      lineend = "round",
      show.legend = FALSE
    )
    if (!color_by_rhat)
      segment_args$color <- get_color("d")
    graph <- graph + do.call("geom_segment", segment_args)

    # point estimate
    point_args <- list(
      mapping = aes_(
        x = ~ m,
        y = ~ y,
        color = if (!color_by_rhat) NULL else ~ rhat,
        fill = if (!color_by_rhat) NULL else ~ rhat
      ),
      size = 4,
      shape = 21
    )
    if (!color_by_rhat) {
      point_args$color <- get_color("dh")
      point_args$fill <- get_color("l")
    }

    if (!no_point_est)
      graph <- graph + do.call("geom_point", point_args)

    if (color_by_rhat)
      graph <- graph +
      scale_color_diagnostic("rhat") +
      scale_fill_diagnostic("rhat")
  }

  graph +
    scale_y_continuous(
      breaks = y,
      labels = parnames,
      limits = c(0.5, n_param + 1)
    ) +
    xlim(x_lim) +
    legend_move(ifelse(color_by_rhat, "top", "none")) +
    yaxis_text(face = "bold") +
    yaxis_title(FALSE) +
    xaxis_title(FALSE)
}


# compute density estimates for mcmc_areas
.compute_dens_i <- function(x, bw, adjust, kernel, n, from, to) {
  args <- c(
    # can't be null
    list(
      x = x,
      from = from,
      to = to,
      n = n
    ),
    # might be null
    bw = bw,
    adjust = adjust,
    kernel = kernel
  )
  do.call("density", args)
}

# compute density estimates for mcmc_joy
.compute_dens_i_joy <- function(x, prob, bw, adjust, kernel, n = 1e3) {
  alpha <- (1 - prob) / 2
  probs <- c(alpha, 1 - alpha)
  bounds <- range(quantile(x, probs))
  xdens <- stats::density.default(x, bw, adjust, kernel, n = n,
                                  from = bounds[1], to = bounds[2])
  data.frame(value = xdens$x, dens = xdens$y)
}

