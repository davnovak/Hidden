## Hidden Markov Models

ResetStats <- function(H, ...) UseMethod('ResetStats', H)

#' Delete computed statistics in a HiddenMarkovModel object
#'
#' Clears fields \code{probs}, \code{most_likely_paths_probs}, \code{most_likely_paths_pointers} and \code{most_likely_path}.
#'
#' @param H \code{HiddenMarkovModel} object, as created by the constructor function \code{HiddenMarkovModel}
#'
#'
ResetStats.HiddenMarkovModel <- function(H) {

  prob_list        <- list(filtered = NULL, predicted = NULL, backward = NULL, smoothed = NULL)
  H$probs          <- list()
  for (i in 1:H$n_observations) {
    H$probs[[i]] <- prob_list
  }
  H$most_likely_paths_probs    <- NULL
  H$most_likely_paths_pointers <- NULL
  H$most_likely_path           <- NULL
  invisible(H)
}

#' Construct a HiddenMarkovModel object
#'
#' Constructs an object for simulation of a hidden Markov model.
#' This can be used to apply the filtration, prediction, smoothing and Viterbi algorithms.
#' In the future, a plotting capability might also be added.
#'
#' @param Tr transition matrix; a \code{K}-by-\code{K} numeric matrix where \code{K} is the number of possible hidden states.
#' @param Em emission matrix; a \code{K}-by-\code{N} numeric matrix where \code{N} is the number of possible observable states.
#' @param X0 initial probability distribution; a column vector of initial probabilities of hidden states at time zero.
#' @param hidden_states; character vector of names of hidden states in order. If not specified, they are named '\code{x1}', '\code{x2}', ... , '\code{xK}'.
#' @param observables; character vector of names of possible observable states in order. If not specified, they are named '\code{e1}', '\code{e2}', ... , '\code{eN}'.
#' @param name; internal name of the \code{HiddenMarkovModel} object. If not specified, it is named '\code{HMM_}' followed by the system date.
#'
HiddenMarkovModel <- function(Tr, Em, X0, hidden_states = NULL, observables = NULL, name = paste0('HMM_', as.character(Sys.Date()))) {

  if (!is.matrix(Tr) && substr(class(Tr), nchar(class(Tr)) - 5, nchar(class(Tr))) != 'Matrix') stop('Tr must be a numeric matrix')
  if (!is.matrix(Em) && substr(class(Em), nchar(class(Em)) - 5, nchar(class(Em))) != 'Matrix') stop('Em must be a numeric matrix')
  if (!is.matrix(X0) && substr(class(X0), nchar(class(X0)) - 5, nchar(class(X0))) != 'Matrix') stop('X0 must be a numeric matrix')
  if (nrow(X0) != ncol(Tr))                                    stop('Number of rows of X0 and number of columns of Tr must be equal (number of observables)')
  if (ncol(Tr) != nrow(Em))                                    stop('Number of columns of Tr and number of rows of Em must be equal (number of hidden states)')
  if (!is.null(hidden_states) && !is.character(hidden_states)) stop('hidden_states, if specified, must be a character vector')
  if (!is.null(observables) && !is.character(observables))     stop('observables, if specified, must be a character vector')

  H      <- new.env(hash = TRUE)
  H$name <- name
  H$Tr   <- Tr
  H$Em   <- Em
  H$X0   <- X0
  H$n_hidden_states <- ncol(Tr)
  H$n_observables   <- ncol(Em)
  H$hidden_states   <- if (is.null(hidden_states)) { paste0('x', 1:H$n_hidden_states) } else { hidden_states }
  H$observables     <- if (is.null(observables))   { paste0('e', 1:H$n_observables) }   else { observables }
  H$probs           <- list()
  colnames(H$Em)    <- H$observables
  colnames(H$Tr)    <- rownames(H$Tr) <- rownames(H$Em) <- rownames(H$X0) <- H$hidden_states

  structure(H, class = 'HiddenMarkovModel')
}

#' Print information about a HiddenMarkovModel object
#'
#' Prints information about the properties of a \code{HiddenMarkovModel} object.
#'
print.HiddenMarkovModel <- function(H) {

  cat('Hidden Markov Model object\n')
  cat('... name: ', H$name, '\n')
  cat('... hidden states names:\n\t',     paste(H$hidden_states, collapse = ', '), '\n')
  cat('... observable states names:\n\t', paste(H$observables, collapse = ', '), '\n')
  cat('... transition matrix:\n')
  print(H$Tr)
  cat('... emission matrix:\n')
  print(H$Em)
  cat('... initial probability distribution:\n')
  print(H$X0)
  if (!is.null(H$observations)) {
    cat('... observations (', H$n_observations, '):\n\t', paste(H$observations, collapse = ' -> '), '\n', sep = '')
  } else {
    cat('... observations: not set')
  }
}

#' Plot a HiddenMarkovModel object diagram
#'
#' This method is not implemented.
#' I might or might not implement it.
#'
plot.HiddenMarkovModel <- function(H, ...) {

  message('Plotting not implemented yet')
  invisible(H)
}

PrintStats <- function(H, ...) UseMethod('PrintStats', H)

#' Report HiddenMarkovModel object statistics
#'
#' Prints computed statistics stored in a \code{HiddenMarkovModel} object: probabilities obtained by forward pass, backward pass, smoothing and the Viterbi algorithm.
#'
PrintStats.HiddenMarkovModel <- function(H) {

  cat('Model name:', H$name, '\n')
  cat('... filtered probabilities:\n')
  p <- do.call(cbind, lapply(H$probs, function(p) p$filtered))
  if (!is.null(dim(p))) colnames(p) <- paste('t=', 1:ncol(p), sep = '')
  print(p)
  cat('... predicted probabilities:\n')
  p <- do.call(cbind, lapply(H$probs, function(p) p$predicted))
  if (!is.null(dim(p))) colnames(p) <- paste('t=', 1:ncol(p), sep = '')
  print(p)
  cat('... backward probabilities:\n')
  p <- do.call(cbind, lapply(H$probs, function(p) p$backward))
  if (!is.null(dim(p))) colnames(p) <- paste('t=', 1:ncol(p), sep = '')
  print(p)
  cat('... smoothed probabilities:\n')
  p <- do.call(cbind, lapply(H$probs, function(p) p$smoothed))
  if (!is.null(dim(p))) colnames(p) <- paste('t=', 1:ncol(p), sep = '')
  print(p)
  if (!is.null(H$most_likely_path)) {
    cat('... most likely hidden states path:\n\t', paste(H$most_likely_path, collapse = ' -> '), sep = '')
  }
}

SetObservations <- function(H, ...) {

  UseMethod('SetObservations', H)
}

#' Add observations to a HiddenMarkovModel object
#'
#' Adds a vector of observations at progressing time points to a \code{HiddenMarkovModel} object.
#' Please beware: this erases any stats you have computed on your model so far.
#'
#' @param H \code{HiddenMarkovModel} object, as created by the constructor function \code{HiddenMarkovModel}.
#' @param values character or integer vector; either a vector of observable state names (taken from the names of observables specified for \code{H}) or a vector of their indices with respect to the order in which they are listed (use \code{print(H)} to check).
#'
SetObservations <- function(H, values) {

  idcs.na <- which(is.na(values))
  na <- (length(idcs.na) > 0)
  v <- if (na) { values[-idcs.na] } else { values }
  if (is.character(values) && any(!values %in% H$observables)) stop('values must only contain names/indices of observables or NAs')
  if (na && (max(which(!is.na(values))) > min(idcs.na)))       stop('No names/indices of observables may occur after the first NA in values')
  if (is.numeric(values) && any(!v %in% 1:H$n_observables))    stop('If indices of observables are given in values, they must be from <1,N> where N is the number of observables')

  H$observations   <- if (is.character(values)) { values } else { H$observables[values] }
  H$n_observations <- length(H$observations)
  ResetStats(H)
  invisible(H)
}

ForwardPass <- function(H = parent.frame()$H, ...) UseMethod('ForwardPass', H)

#' Forward pass in a hidden Markov model
#'
#' This is an internal function of the package \code{hidden}.
#' It implements the recursive forward algorithm for computing posterior probability distributions of hidden states in hidden Markov models.
#' It is called by the methods \code{Filter.HiddenMarkovModel} and \code{Smooth.HiddenMarkovModel}.
#'
ForwardPass.HiddenMarkovModel <- function(H             = parent.frame()$H,
                                          idx           = parent.frame()$idx - 1,
                                          verbose       = parent.frame()$verbose) {

  if (idx > 1) {
    ForwardPass()
  } else {
    H$Fwd <- H$X0
  }
  H$Fwd <- H$Tr %*% H$Fwd
  if (verbose) {
    message(paste0('Step ', idx))
    message(paste0('Prob distribution after transition: ', paste(round(H$Fwd, digits = 4), collapse = ', ')))
  }
  obs <- H$observations[idx]
  H$prediction <- is.na(obs)
  if (!H$prediction) {
    H$Fwd <- H$Fwd * H$Em[, obs]
    H$Fwd <- H$Fwd / sum(H$Fwd)
    if (verbose) {
      message(paste0('Observation ', obs))
      message(paste0('Prob distribution after observation: ', paste(round(H$Fwd, digits = 4), collapse = ', ')))
    }
  } else {
    if (verbose) {
      message('Observation UNKNOWN')
    }
  }
  slot_name <- if (H$prediction) { 'predicted' } else { 'filtered' }
  H$probs[[idx]][[slot_name]] <- H$Fwd
  invisible(H)
}

BackwardPass <- function(H = parent.frame()$H, ...) UseMethod('BackwardPass', H)

#' Backward pass in a hidden Markov model
#'
#' This is an internal function of the package \code{hidden}.
#' It implements the recursive backward algorithm for computing posterior probability distributions of hidden states in hidden Markov models.
#' It is called by the method \code{Smooth.HiddenMarkovModel}.
#'
BackwardPass.HiddenMarkovModel <- function(H             = parent.frame()$H,
                                           idx           = parent.frame()$idx - 1,
                                           verbose       = parent.frame()$verbose) {

  if (idx > 0) {
    BackwardPass()
  } else {
    H$Bwd <- Matrix(rep(1, H$n_hidden_states), nr = H$n_hidden_states)
    H$probs[[H$n_observations]][['backward']] <- H$Bwd
    return(invisible(H))
  }
  obs <- rev(H$observations)[H$n_observations - idx + 1]
  if (is.na(obs)) stop(paste0('Backward pass failed: observation at position ', idx, ' is NA'))
  H$Bwd <- H$Bwd * H$Em[, obs]
  H$Bwd <- H$Bwd / sum(H$Bwd)
  if (verbose) {
    message(paste0('Step ', idx))
    message(paste0('Observation ', obs))
    message(paste0('Prob distribution after observation: ', paste(round(H$Bwd, digits = 4), collapse = ', ')))
  }
  H$Bwd <- H$Tr %*% H$Bwd
  if (verbose) {
    message(paste0('Prob distribution after transition: ', paste(round(H$Bwd, digits = 4), collapse = ', ')))
  }
  H$probs[[H$n_observations - idx]][['backward']] <- H$Bwd
  invisible(H)
}

Filter <- function(H, ...) UseMethod('Filter', H)

#' Apply filtering to a HiddenMarkovModel object
#'
#' Applies the filtering algorithm by a recursive forward pass through a hidden Markov model.
#' This lets us compute the posterior probability distribution of the most recent state given all evidence up to a current point in time.
#' A vector of observations needs to be included with the model for this algorithm to work (use \code{SetObservations.HiddenMarkovModel}).
#' It is acceptable for the vector of observations to include \code{NA}s at the end (then, we apply prediction instead of filtration at some point).
#'
#' After the algorithm finishes, see \code{YourHiddenMarkovModelObject$probs[[YourTimePointIndex]]$filtered} to retrieve the resulting probability distribution.
#'
#' @param H \code{HiddenMarkovModel} object, as created by the constructor function \code{HiddenMarkovModel}.
#' @param idx index of the last time step for which to run filtration (filtration probabilities for the ones preceding it will be stored automatically).
#' @param verbose logical; whether to report progression of the forward pass.
#'
Filter.HiddenMarkovModel <- function(H, idx, verbose = TRUE) {

  if (is.null(H$observations)) stop('No observations available in the model')
  if (!is.numeric(idx) || length(idx) != 1 || (!idx %in% 1:H$n_observations)) stop('idx must be a numeric index of time point where filtration/prediction stops')
  if (verbose) {
    message('## Forward pass')
    message('Initial probs: ', paste(H$X0, collapse = ', '))
  }
  ForwardPass(H, idx, verbose)
  H$Fwd        <- NULL
  H$prediction <- NULL
  gc()

  invisible(H)
}

Smooth <- function(H, ...) UseMethod('Smooth', H)

#' Apply smoothing to a HiddenMarkovModel object
#'
#' Applies the smoothing algorithm by a recursive forward pass and backward passthrough a hidden Markov model.
#' This lets us compute the posterior probability distribution of a past state given all evidence up to a current point in time.
#' A vector of observations needs to be included with the model for this algorithm to work (use \code{SetObservations.HiddenMarkovModel}).
#' It is not acceptable for the vector of observations to include any \code{NA}s.
#'
#' After the algorithm finishes, see \code{YourHiddenMarkovModelObject$probs[[YourTimePointIndex]]$smoothed} to retrieve the resulting probability distribution.
#' If interested in the probability distributions computed in the forward pass, see the slot \code{filtered} instead.
#' To see the probability distributions computed in the backward pass, see the slot \code{backward}.
#'
#' @param H \code{HiddenMarkovModel} object, as created by the constructor function \code{HiddenMarkovModel}.
#' @param idcs a single index or a vector of indices of time points for which smoothing is to be done.
#' @param verbose logical; whether to report progression of the forward pass and the backward pass.
#'
Smooth.HiddenMarkovModel <- function(H, idcs, verbose = TRUE) {

  if (is.null(H$observations) || any(is.na(H$observations)))   stop('A vector of observations with no NAs is required')
  if (!is.numeric(idcs) || any(!idcs %in% 1:H$n_observations)) stop('idcs must be a single index or a vector of indices of time points for which smoothing is to be done')
  if (verbose) {
    message('## Forward pass')
    message('Initial probs: ', paste(H$X0, collapse = ', '))
  }
  ForwardPass(H, H$n_observations, verbose)
  H$Fwd        <- NULL
  H$prediction <- NULL
  gc()
  if (verbose) {
    message('## Backward pass')
  }
  depth <- H$n_observations - min(idcs)
  BackwardPass(H, depth, verbose)
  H$Bwd        <- NULL
  gc()
  for (idx in idcs) {
    prob <- H$probs[[idx]][['filtered']] * H$probs[[idx]][['backward']]
    prob <- prob / sum(prob)
    if (verbose) message(paste0('Smoothed probs at time t=', idx, ': ', paste(prob, collapse = ', ')))
    H$probs[[idx]][['smoothed']] <- prob
  }
  invisible(H)
}

MostLikelyPath <- Viterbi <- function(H, ...) UseMethod('Viterbi', H)

#' Apply the Viterbi algorithm to a HiddenMarkovModel object
#'
#' Applies the Viterbi algorithm using dynamic programming on a hidden Markov model.
#' This lets us compute the most likely progression (path) of hidden states given a sequence of observable states.
#' A vector of observations needs to be included with the model for this algorithm to work (use \code{SetObservations.HiddenMarkovModel}).
#' It is not acceptable for the vector of observations to include any \code{NA}s.
#'
#' After the algorithm finishes, see \code{YourHiddenMarkovModelObject$most_likely_path} to retrieve the most likely sequence of hidden states.
#' If interested in the dynamic programming matrix with path probabilities, see the slot \code{most_likely_paths_probs} instead.
#' To view the matrix with pointers for backtracking in the Viterbi algorithm, see the slot \code{most_likely_path_pointers}.
#'
#' @param H \code{HiddenMarkovModel} object, as created by the constructor function \code{HiddenMarkovModel}.
#' @param idx.end an index of the time points where path terminates.
#' @param verbose logical; whether to report progression of the Viterbi algorithm.
#'
Viterbi.HiddenMarkovModel <- function(H, idx.end = H$n_observations) {

  if (is.null(H$observations))               stop('No known observations available for model')
  if (any(is.na(H$observations[1:idx.end]))) stop('Observations on the path may not contain NAs')

  Ncol <- H$n_observations + 1
  T.max <- T.argmax <- Matrix(0, nrow = H$n_hidden_states, ncol = Ncol)
  rownames(T.max) <- rownames(T.argmax) <- H$hidden_states
  colnames(T.max) <- colnames(T.argmax) <- paste('t=', 0:H$n_observations, sep = '')

  T.max[, 1] <- H$Em %*% H$X0
  for (idx.time in 2:(idx.end + 1)) {
    for (idx.state in 1:H$n_hidden_states) {
      probs.transition             <- H$Tr[, idx.state]                             # probabilities of transitions given current state
      probs.current_distribution   <- T.max[, idx.time - 1]                         # maximum probabilities up to here
      probs.emission               <- H$Em[idx.state, H$observations[idx.time - 1]] # fixed probability of emission
      routes                       <- probs.transition * probs.current_distribution * probs.emission
      T.max[idx.state, idx.time]    <- max(routes)
      T.argmax[idx.state, idx.time] <- which.max(routes)
    }
  }

  idcs.likeliest <- states.likeliest <- vector(mode = 'list', length = Ncol)
  idcs.likeliest[[Ncol]]   <- which.max(T.max[, Ncol]) # most probable end state index
  states.likeliest[[Ncol]] <- rownames(T.max)[idcs.likeliest[[Ncol]]]
  for (idx.time in Ncol:2) {
    idcs.likeliest[[idx.time - 1]]   <- T.argmax[idcs.likeliest[[idx.time]], idx.time]
    states.likeliest[[idx.time - 1]] <- rownames(T.max)[idcs.likeliest[[idx.time]]]
  }

  H$most_likely_paths_probs    <- T.max
  H$most_likely_paths_pointers <- T.argmax
  H$most_likely_path           <- unlist(states.likeliest)
  gc()
  invisible(H)
}

