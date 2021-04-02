#' @export
makeContrasts <- function(..., contrasts=NULL, levels)
#  Construct matrix of custom contrasts
#  Gordon Smyth
#  30 June 2003.  Last modified 21 August 2006.
{
  e <- substitute(list(...))
  if(is.factor(levels)) levels <- levels(levels)
  if(!is.character(levels)) levels <- colnames(levels)
  #if(any(levels != make.names(levels))) stop("Parameter names must by syntactically valid names in R")
  n <- length(levels)
  if(n < 1) stop("No levels to construct contrasts from")
  indicator <- function(i,n) {
    out <- rep(0,n)
    out[i] <- 1
    out
  }
  levelsenv <- new.env()
  for (i in 1:n) assign(levels[i], indicator(i,n), pos=levelsenv)

#  Contrasts given as character vector
  if(!is.null(contrasts)) {
    if(length(e)>1) stop("Can't specify both ... and contrasts")
    e <- as.character(contrasts)
    ne <- length(e)
    cm <- matrix(0,n,ne, dimnames=list(Levels=levels,Contrasts=e))
    if(ne==0) return(cm)
    for (j in 1:ne) {
      ej <- parse(text=e[j])
      cm[,j] <- eval(ej, envir=levelsenv)
    }
    return(cm)
  }

#  Contrasts given as list of expressions
  ne <- length(e)
  enames <- names(e)[2:ne]
  easchar <- as.character(e)[2:ne]
  if(is.null(enames))
    cn <- easchar
  else
    cn <- ifelse(enames=="",easchar,enames)
  cm <- matrix(0,n,ne-1, dimnames=list(Levels=levels,Contrasts=cn))
  if(ne < 2) return(cm)
  for (j in 1:(ne-1)) {
    ej <- e[[j+1]]
    if(is.character(ej)) ej <- parse(text=ej)
    ej <- eval(ej, envir=levelsenv)
#    Character variable
    if(!is.numeric(ej)) {
      colnames(cm)[j] <- as.character(ej)[1]
      if(is.character(ej)) ej <- parse(text=ej)
      ej <- eval(ej, envir=levelsenv)
    }
    cm[,j] <- ej
  }
  cm
}
