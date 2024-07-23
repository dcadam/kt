## Various Utility functions 

custom_min_max_normalize <- function(x, a = 0, b = 1.57) {
  a + (x - min(x)) * (b - a) / (max(x) - min(x))
}

make_tbl <- function(input) {
  
  kk <- c(0.1, 0.5, 1.0, 10, 100)
  output <- tibble()
  for (i in 1:length(input)) {
    tmp <- input[[i]] |> 
      mutate(k = kk[i])
    
    output <- bind_rows(output, tmp)
  }
  
  return(output)
}


extract_est_mix <- function(list_in, n_lim = 300) {
  
  pb <- txtProgressBar(min = 1, max = length(list_in), style = 3)
  
  out <- vector(mode = "list", length = length(list_in))
  for (i in 1:length(list_in)) {
    
    tmp <- list_in[[i]]
    
    data_steps <- tibble(window = names(tmp[[1]]))
    
    
    
    data_out <- vector(mode = "list", length = nrow(data_steps))
    for (d in 1:nrow(data_steps)) {
      
      data_out[[d]] <- data.table()
      
      for (j in 1:length(tmp)) {
        
        tryCatch( #error catch
          {
            
            tmp_out <- tmp[[j]][[d]] |> 
              mutate(window = data_steps$window[[d]])
            
            
            data_out[[d]] <- rbindlist(list(data_out[[d]], tmp_out))
            
            
          },
          error = function(e) {}
        )
        
      }
      
      # if (d != 1) {
      #   
      #   data_out[[d]] <- data_out[[d]] |> 
      #     rename(t = t_start)
      #   
      # }
    }
    
    data_out <- map(data_out, function(x) {
      
      x |> 
        mutate(r1 = exp(r1),
               r2 = exp(r2)) |> 
        mutate(rc = (r1+r2) / 2) |> 
        mutate(kc = (k1+k2) / 2) |> 
        group_by(t) |> 
        mutate(n = n()) |> 
        filter(n >= n_lim) |> 
        reframe(
          Rt =  quantile(rt, 0.5, na.rm = TRUE),
          Rt_lower = quantile(rt, 0.025, na.rm = TRUE),
          Rt_upper = quantile(rt, 0.975, na.rm = TRUE),
          Rc =  quantile(rc, 0.5, na.rm = TRUE),
          Rc_lower = quantile(rc, 0.025, na.rm = TRUE),
          Rc_upper = quantile(rc, 0.975, na.rm = TRUE),
          Kc =  quantile(kc, 0.5, na.rm = TRUE),
          Kc_lower = quantile(kc, 0.05, na.rm = TRUE),
          Kc_upper = quantile(kc, 0.95, na.rm = TRUE)) |> 
        mutate(data_step = unique(x$window))
      
    })
    
    
    
    out[[i]] <- list_rbind(data_out)
    
    setTxtProgressBar(pb, i)
    
  }
  
  return(out)
  
  
}


extract_est <- function(list_in, n_lim = 300) {
  
  out <- vector(mode = "list", length = length(list_in))
  for (i in 1:length(list_in)) {
    
    tmp <- list_in[[i]]
    
    data_steps <- tibble(window = names(tmp)) |>
      separate(window, into = c("window", "sporadic"), sep = "_sporadic_")
    
    pb <- txtProgressBar(min = 1, max = nrow(data_steps), style = 3)
    
    data_out <- vector(mode = "list", length = nrow(data_steps))
    for (d in 1:nrow(data_steps)) {
      
      data <- tmp[[d]]
      
      ## combine all 300 epidemics and 25 subsamples
      rxk_epidemics <- map(data, function(x) { 
        x |>
          list_rbind()
      }) |>
        list_rbind()
      
      
      if (str_detect(string = data_steps$window[[d]], pattern = "daily")) {
        rxk_epidemics <- rxk_epidemics
      } else {
        rxk_epidemics <- rxk_epidemics |>
          rename(t = t_start)
      }
      
      
      
      rxk_qunatiles <- rxk_epidemics |>
        group_by(t, o_prop, ct_prob) |>
        mutate(n = n()) |> 
        filter(n >= n_lim) |> 
        reframe(
          Rc =  mean(r, na.rm = TRUE),
          Rc_lower = quantile(r, 0.025, na.rm = TRUE),
          Rc_upper = quantile(r, 0.975, na.rm = TRUE),
          kc =  quantile(k, 0.5, na.rm = TRUE),
          kc_lower = quantile(k, 0.05, na.rm = TRUE),
          kc_upper = quantile(k, 0.95, na.rm = TRUE)
        ) |>
        mutate(
          sporadic = data_steps$sporadic[[d]],
          data_step = data_steps$window[[d]]
        )
      
      
      
      
      data_out[[d]] <- rxk_qunatiles
      
      setTxtProgressBar(pb, d)
      
    }
    
    close(pb)
    
    out[[i]] <- list_rbind(data_out)
    
  }
  
  return(out)
  
}



sporadicAs <- function(sporadic_in, sporadic_as) {
  
  if (sporadic_as == 0) {
    
    sporadic_out <- sporadic_in
    
    return(sporadic_out)
    
  } 
  
  if (sporadic_as == 1) {
    
    sporadic_out <- map(sporadic_in, function(x) {
      
      sporadic <- x |> 
        mutate(t1 = as.numeric(source),
               t2 = as.numeric(case_id)) |> 
        filter(!t1 %in% t2 & obs_ct == 0) |> 
        pull(case_id)
      
      map_out <- x |> 
        mutate(obs_ct = case_when(case_id %in% sporadic ~ 1,
                                  TRUE ~ obs_ct))
      
      return(map_out)
      
      
    })
    
  }
  
  if (sporadic_as == "ex") {
    
    sporadic_out <- map(sporadic_in, function(x) {
      
      sporadic <- x |> 
        mutate(t1 = as.numeric(source),
               t2 = as.numeric(case_id)) |> 
        filter(!t1 %in% t2
               & obs_ct == 0) |> 
        pull(case_id)
      
      map_out <- as.data.table(x)[!case_id %in% sporadic] |> 
        as_tibble()
      
      return(map_out)
      
    })
    
  }
  
  return(sporadic_out)
  
}


FLXMRnegbin2 <- function (formula = . ~ ., 
                          theta = NULL, 
                          offset = NULL, 
                          method = "BFGS",
                          lower = -Inf, upper = Inf)
{
  .theta <- theta
  nbrefit <- function(x, y, w) {
    fit <- c(glm.fit(x, y, weights = w, offset = offset, 
                     family = MASS::negative.binomial(theta)), list(call = sys.call(), 
                                                                    offset = offset, control = eval(formals(glm.fit)$control), 
                                                                    method = "weighted.glm.fit"))
    fit$df.null <- sum(w) + fit$df.null - fit$df.residual - 
      fit$rank - is.null(.theta)
    fit$df.residual <- sum(w) - fit$rank - is.null(.theta)
    fit$x <- x
    fit
  }
  z <- methods::new("FLXMRglm", weighted = TRUE, formula = formula, 
                    name = "FLXMRglm: negative.binomial", offset = offset, 
                    family = "negative.binomial", refit = nbrefit)
  z@preproc.y <- function(x) {
    if (ncol(x) > 1L) 
      stop(paste("for the", family, "family y must be univariate"))
    x
  }
  z@defineComponent <- if (is.null(.theta)) {
    expression({
      predict <- function(x, ...) {
        dotarg <- list(...)
        if ("offset" %in% names(dotarg)) offset <- dotarg$offset
        p <- x %*% coef
        if (!is.null(offset)) p <- p + offset
        exp(p)
      }
      logLik <- function(x, y, ...) suppressWarnings(dnbinom(y, 
                                                             mu = predict(x, ...), size = theta, log = TRUE))
      methods::new("FLXcomponent", parameters = list(coef = coef, 
                                                     theta = theta), logLik = logLik, predict = predict, 
                   df = df)
    })
  }
  else {
    as.expression(substitute({
      predict <- function(x, ...) {
        dotarg <- list(...)
        if ("offset" %in% names(dotarg)) offset <- dotarg$offset
        p <- x %*% coef
        if (!is.null(offset)) p <- p + offset
        exp(p)
      }
      logLik <- function(x, y, ...) suppressWarnings(dnbinom(y, 
                                                             mu = predict(x, ...), size = theta, log = TRUE))
      methods::new("FLXcomponent", parameters = list(coef = coef), 
                   logLik = logLik, predict = predict, df = df)
    }, as.environment(list(theta = .theta))))
  }
  z@fit <- function(x, y, w, component) {
    if (is.null(component$theta)) {
      df <- ncol(x)
      theta <- if (is.null(.theta)) 
        1
      else .theta
      cf <- glm.fit(x, y, weights = w, family = MASS::negative.binomial(theta), 
                    offset = offset, start = component$coef)$coefficients
    }
    else {
      df <- ncol(x) + 1
      if (is.null(offset)) 
        offset <- 0
      nll <- function(par) {
        beta <- par[-df]
        theta <- exp(par[df])
        mu <- exp(drop(x %*% beta + offset))
        suppressWarnings(-sum(w * dnbinom(y, mu = mu, 
                                          size = theta, log = TRUE)))
      }
      gr <- function(par) {
        beta <- par[-df]
        theta <- exp(par[df])
        mu <- exp(drop(x %*% beta + offset))
        gr <- drop(y - mu * (y + theta)/(mu + theta))
        colSums(-w * cbind(gr * x, theta * (digamma(y + 
                                                      theta) - digamma(theta) + log(theta) + 1 - 
                                              log(mu + theta) - (y + theta)/(mu + theta))))
      }
      start <- c(component$coef, component$theta)
      if (length(start) < df) 
        start <- c(glm.fit(x, y, weights = w, family = MASS::negative.binomial(1), 
                           offset = offset)$coefficients, 0)
      opt <- optim(par = start, fn = nll, gr = gr, method = method, 
                   lower = lower,
                   upper = upper)
      cf <- opt$par[-df]
      theta <- exp(opt$par[df])
    }
    with(list(coef = cf, theta = theta, df = ncol(x) + is.null(.theta)), 
         eval(z@defineComponent))
  }
  z
}



generateOffspring<- function(case_data, pair_data, sporadic_logical) {
  
  output <- pair_data %>%
    group_by(infector_case) %>%
    count()
  
  if (sporadic_logical == TRUE) {
    
    sporadic <- case_data %>%
      filter(cluster_classification == "Sporadic local case" | 
               classification == "sporadic - unmatched" |
               classification == "sporadic - matched") %>%
      nrow()
    infectees <- pair_data %>%
      dplyr::select(infector_case, infectee_case) %>%
      gather() %>%
      filter(key == "infectee_case")
    infectors <- pair_data %>%
      dplyr::select(infector_case, infectee_case) %>%
      gather() %>%
      filter(key == "infector_case")
    duplicate <- infectors %>%
      left_join(., infectees, by = "value") %>%
      filter(key.y != "NA") %>%
      dplyr::select(value) %>%
      distinct()
    nterminal_infectees <- infectees %>%
      dplyr::select(value) %>%
      filter(!value %in% duplicate$value) %>%
      transmute(case_no = (value)) %>%
      nrow() + sporadic
    output <- enframe(c(output$n, rep(0, nterminal_infectees)))
    
    return(output)
  }
  
  else {
    
    infectees <- pair_data %>%
      dplyr::select(infector_case, infectee_case) %>%
      gather() %>%
      filter(key == "infectee_case")
    infectors <- pair_data %>%
      dplyr::select(infector_case, infectee_case) %>%
      gather() %>%
      filter(key == "infector_case")
    duplicate <- infectors %>%
      left_join(., infectees, by = "value") %>%
      filter(key.y != "NA") %>%
      dplyr::select(value) %>%
      distinct()
    nterminal_infectees <- infectees %>%
      dplyr::select(value) %>%
      filter(!value %in% duplicate$value) %>%
      transmute(case_no = (value)) %>%
      nrow()
    output <- enframe(c(output$n, rep(0, nterminal_infectees)))
    
    return(output)
    
    
  }
  
}


propresponsible <- function(R0, k, prop) {
  qm1 <- qnbinom(1 - prop, k + 1, mu = R0 * (k + 1) / k)
  remq <- 1 - prop - pnbinom(qm1 - 1, k + 1, mu = R0 * (k + 1) / k)
  remx <- remq / dnbinom(qm1, k + 1, mu = R0 * (k + 1) / k)
  q <- qm1 + 1
  1 - pnbinom(q - 1, k, mu = R0) - dnbinom(q, k, mu = R0) * remx
}


gen_date_pairs <- function(case_data, pairs_data) {
  a <- case_data %>%
    mutate(epi_week = epiweek(t_inf))
  b <- pairs_data %>%
    transmute(
      infector_case = as.character(infector_case),
      infectee_case = as.character(infectee_case)
    )
  output <- a %>%
    mutate(
      infector_case = as.character(case_no),
      infector_epi_date = t_inf
    ) %>%
    right_join(., b, by = "infector_case") %>%
    dplyr::select(infector_case, infector_epi_date, infectee_case)
  output <- a %>%
    mutate(
      infectee_case = as.character(case_no),
      infectee_epi_date = t_inf
    ) %>%
    right_join(., output, by = "infectee_case") %>%
    dplyr::select(infector_case, infectee_case, infector_epi_date, infectee_epi_date)
  output <- output %>%
    mutate(infector_epi_date = case_when(
      is.na(infector_epi_date) ~ infectee_epi_date,
      TRUE ~ infector_epi_date
    ))
  output <- output %>%
    mutate(infectee_epi_date = case_when(
      is.na(infectee_epi_date) ~ infector_epi_date,
      TRUE ~ infectee_epi_date
    ))
  output <- output %>%
    mutate(
      infector_epi_week = epiweek(infector_epi_date),
      infectee_epi_week = epiweek(infectee_epi_date)
    )
}


StatBinscatter <- ggplot2::ggproto(
  "StatBinscatter", 
  Stat,
  compute_group = function(data, scales, bins = 10) {
    bins     <- min(floor(nrow(data)/10), bins)
    x_bin    <- ggplot2::cut_number(data$x + 1e-12*runif(nrow(data)), bins)
    x_means  <- stats::ave(data$x, x_bin, FUN = mean)
    y_means  <- stats::ave(data$y, x_bin, FUN = mean)
    y_se     <- stats::ave(data$y, x_bin, FUN = sd)
    y_obs    <- stats::ave(data$y, x_bin, FUN = length)
    result   <- data.frame(x    = x_means, 
                           y    = y_means, 
                           ymax = y_means + 1.96*y_se/sqrt(y_obs),
                           ymin = y_means - 1.96*y_se/sqrt(y_obs))
    result   <- unique(result)
    return(result)
  },
  required_aes = c("x", "y")
)


stat_binscatter <- function(mapping = NULL, data = NULL, geom = "point",
                            position = "identity", na.rm = FALSE, show.legend = NA, 
                            inherit.aes = TRUE, ...) {
  ggplot2::layer(
    stat = StatBinscatter, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}



facet_wrap_custom <- function(..., scale_overrides = NULL) {
  # take advantage of the sanitizing that happens in facet_wrap
  facet_super <- facet_wrap(...)
  
  # sanitize scale overrides
  if(inherits(scale_overrides, "scale_override")) {
    scale_overrides <- list(scale_overrides)
  } else if(!is.list(scale_overrides) || 
            !all(vapply(scale_overrides, inherits, "scale_override", FUN.VALUE = logical(1)))) {
    stop("scale_overrides must be a scale_override object or a list of scale_override objects")
  }
  
  facet_super$params$scale_overrides <- scale_overrides
  
  ggproto(NULL, CustomFacetWrap,
          shrink = facet_super$shrink,
          params = facet_super$params
  )
}


scale_override <- function(which, scale) {
  if(!is.numeric(which) || (length(which) != 1) || (which %% 1 != 0)) {
    stop("which must be an integer of length 1")
  }
  
  if(is.null(scale$aesthetics) || !any(c("x", "y") %in% scale$aesthetics)) {
    stop("scale must be an x or y position scale")
  }
  
  structure(list(which = which, scale = scale), class = "scale_override")
}

CustomFacetWrap <- ggproto(
  "CustomFacetWrap", FacetWrap,
  init_scales = function(self, layout, x_scale = NULL, y_scale = NULL, params) {
    # make the initial x, y scales list
    scales <- ggproto_parent(FacetWrap, self)$init_scales(layout, x_scale, y_scale, params)
    
    if(is.null(params$scale_overrides)) return(scales)
    
    max_scale_x <- length(scales$x)
    max_scale_y <- length(scales$y)
    
    # ... do some modification of the scales$x and scales$y here based on params$scale_overrides
    for(scale_override in params$scale_overrides) {
      which <- scale_override$which
      scale <- scale_override$scale
      
      if("x" %in% scale$aesthetics) {
        if(!is.null(scales$x)) {
          if(which < 0 || which > max_scale_x) stop("Invalid index of x scale: ", which)
          scales$x[[which]] <- scale$clone()
        }
      } else if("y" %in% scale$aesthetics) {
        if(!is.null(scales$y)) {
          if(which < 0 || which > max_scale_y) stop("Invalid index of y scale: ", which)
          scales$y[[which]] <- scale$clone()
        }
      } else {
        stop("Invalid scale")
      }
    }
    
    # return scales
    scales
  }
)

breaks_fun <- function(x) {
  count <<- count + 1L
  switch(
    count,
    c(0.01, 0.1, 1),
    c(0.1, 1, 10),
    c(0.1, 1, 10),
    c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000),
    c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000)
  )
}




roundUpNice <- function(x, nice=c(1,2,3,4,5,6,8,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

