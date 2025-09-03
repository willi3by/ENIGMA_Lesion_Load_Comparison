pp_check.blrm <- function(modblrm, 
                          type = "dens_overlay", 
                          ndraws = 20,
                          group = NULL) {
  
  # Validate model class
  if (!any(class(modblrm) %in% c("blrm", "rmsb"))) {
    stop("rms object must be of class blrm", call. = FALSE)
  }
  
  # Validate type
  valid_types <- sub("^ppc_", "", as.character(bayesplot::available_ppc("")))
  if (!type %in% valid_types) {
    stop("Type '", type, "' is not a valid ppc type. ",
         "Valid types are:\n", paste(valid_types, collapse = " , "))
  }
  
  ppc_fun <- get(paste0("ppc_", type), asNamespace("bayesplot"))
  
  # Get data the old way (Option 1) for maximum compatibility
  newdata <- eval(modblrm$call$data)[all.vars(modblrm$sformula)]
  if (anyNA(newdata)) {
    warning("NA responses in sample")
    newdata <- newdata[complete.cases(newdata), ]
  }
  
  valid_vars  <- modblrm$Design$name 
  if ("group" %in% names(formals(ppc_fun))) {
    if (is.null(group)) stop("Argument 'group' is required for ppc type '", type, "'.")
    if (!group %in% valid_vars) stop("Variable '", group, "' not found in data.")
  }
  
  # Extract response
  y <- newdata[, all.vars(modblrm$sformula)[1]]
  y_length <- length(unique(y))
  
  # Binary case
  if (y_length == 2) {
    pred_binary <- predict(modblrm, newdata, fun = plogis, funint = FALSE, 
                           posterior.summary = "all")
    pred_binary <- pred_binary[seq_len(ndraws), , drop = FALSE]
    
    # Vectorized binomial sampling
    yrep <- matrix(
      rbinom(length(pred_binary), 1, pred_binary),
      nrow = nrow(pred_binary)
    )
    
  } else {
    # Ordinal / continuous case
    pred_ordinal <- predict(modblrm, newdata, type = "fitted.ind", 
                            posterior.summary = "all")
    pred_ordinal <- pred_ordinal[seq_len(ndraws), , , drop = FALSE]
    
    draws <- dim(pred_ordinal)[1]
    obs   <- dim(pred_ordinal)[2]
    cats  <- dim(pred_ordinal)[3]
    
    # Flatten to [draws*obs, cats]
    pred_long <- matrix(pred_ordinal, nrow = draws * obs, ncol = cats)
    
    # Bulk multinomial sampling
    samples <- max.col(t(apply(pred_long, 1, function(p) rmultinom(1, 1, p))),
                       ties.method = "first")
    
    # Map back to Y levels
    ylevels <- modblrm$ylevels
    yrep <- matrix(as.numeric(ylevels[samples]), nrow = draws, ncol = obs)
  }
  
  # Prepare args for ppc function
  ppc_args <- list(y, yrep)
  if (!is.null(group)) ppc_args$group <- newdata[[group]]
  
  do.call(ppc_fun, ppc_args)
}
