## MARKOV.R:
## Functions for dealing with transition matrices representing Markov diagrams.
##
## The input format for 'matrix' is tttdlist$ttarray, where: tttdlist
## contains: (1) a TTARRAY: the vocalisation matrix proper, in which
## all rows sum to 1, and (2) TDARRAY: an absolute (i.e. based on
## real-valued time intervals rather than sample counts) turn duration
## vector which corresponds to the static distribution for TTARRAY
## (i.e. TTARRAY^n[], the matrix representation of the
## Chapman-Kolmogorov equation, as n approaches INF) See vocalgraphs.R
## $ttarray format
##
### This program is free software; you can redistribute it and/or
### modify it under the terms of the GNU General Public License
### as published by the Free Software Foundation; either version 2
### of the License, or (at your option) any later version.
###
### (c) 2017 S. Luz (luzs@acm.org)

pauseTypes <- c('Pause', 'SwitchingPause', 
                'GrpPause', 'GrpSwitchingPause') 
categories <- c('Vocalisation', 'SwitchingVocalisation', 'Pause', 'SwitchingPause',
                'GrpPause', 'GrpSwitchingPause',
                'GrpVocalisation') 
catDescriptions <- array(dim=length(categories),dimnames=list(categories),
                         data=c('Vocalisation', 'Switching Vocalisation', 'Pause', 'Switching pause',
                           'Group pause', 'Group switching Pause',
                           'Group vocalisation'))

##' Compute the stationary distribution for a Markov diagram
##'
##' Return static matrix (i.e. the stationary distribution) for the
##' Markov process represented by the given adjacency matrix. In the
##' particular case of vocaldia's, each column should roughly
##' correspond to the amount of time a speaker held the floor for).
##' Of course, not all Markov chains converge, an example being:
##' \preformatted{
##'            1
##'      /----->-------\
##'     A               B
##'      \----<--------/
##'            1
##'
##' which gives
##'
##' .      | 0  1 |             | 0x0+1x1  0x1+1x0|   | 1  0 |
##' .  M = | 1  0 |  and  M^2 = | 1x0+0x1  1x1+1x0| = | 0  1 |
##'
##' }
## whose matrix representation is:
## \deqn{M= \left( \begin{array}{cc}
##   0 & 1 \\
##   1 & 0 \end{array} \right)\]
## and
## \[ M^2 = \left( \begin{array}{cc}
##    0\times 0+1 \times 1  0 \times 1+1 \times 0 \\
##    1\times 0+0\times 1  1\times 1+1\times 0 \end{array} \right)
##    =
##  \left( \begin{array}{cc}
##   1 & 0 \\
##   0 & 1 \end{array} \right)\]
##  }
##' @title staticMatrix Iterate until transition probabilities converge (or give up).
##' @param matrix an adjecency matrix of trnasition probabilities
##' @param limit maximum number of iterations until we give up on
##'     convergence
##' @param digits the number of decimal places to compare
##' @param history if TRUE, keep track of all matrix products
##' @return a matrixseries object; that is, a list where each element
##'     is either the initial matrix or the product of the two
##'     preceding matrices
##' @examples
##' data(vocdia)
##' x2 <- staticMatrix(vocmatrix$ttarray, digits=4, history=TRUE)
##' ## original matrix
##' round(x2[[1]],3)
##' ## stationary matrix (M^139)
##' round(x2[[length(x2)]],3)
##' @export
staticMatrix <- function (matrix, limit=1000, digits=4, history=F) {
  exp <- 2
  ma <- matrix
  mb <- ma %*% matrix
  if (history)
      mseries <- list(ma, mb)
  while (exp < limit && !all(round(ma,digits=digits) == round(mb, digits=digits))) {
    exp <- exp + 1
    ma <- mb
    mb <- ma %*% matrix
    if (history)
      mseries[[exp]] <- mb
  }
  ## when history isn't specified, we only keep first and last iterations
  if (!history)
      mseries <- list(matrix, mb)
  if (exp == limit)
    print(paste("staticMatrix: exp limit reached before convergence at matrix^",exp, sep=""))
  else
    print(paste("staticMatrix: values converged for matrix^",exp, sep=""))
  class(mseries) <- 'matrixseries'
  return(mseries)
}

## 
##' Matrix exponentials
##'
##' A (sort of) exponential function for matrix multiplication (to be
##' used with \code{\link{staticMatrix}}).
##' @title matrixExp: raise matrix to exp.
##' @param matrix a matrix
##' @param exp the power to which matrix will be raised
##' @param mmatrix a placeholder. 
##' @return matrix^exp
##' @examples
##' data(vocdia)
##' matrixExp(vocmatrix$ttarray, 3)
##' @export
matrixExp <- function(matrix, exp, mmatrix=matrix) {
  if (exp < 1){
    warning("matrixExp only accepts positive values")
  }
  i <- 1
  #m <-  matrix
  while (i  < exp){
    i <- i + 1
    mmatrix <- mmatrix %*% matrix
  }
  return(mmatrix)
}

##' Visualise  convergence properties of vocalisation graphs
##'
##' A 'toy' for visualisation of convergence properties of
##' vocalisation graphs. Plot the convergence paths of each
##' Vocalisation event (i.e. each row-column transition probability,
##' grouped by colour according to the inciding node)
##' @title plotConvergence: plots Markov diagram convergence.
##' @param x an object of class matrixseries; a list where the
##'     \eqn{i^{th}} element corresponds to \eqn{M^i}.
##' @param par graphic parameters alist
##' @param interact if TRUE, pauses the drawing after each node.
##' @param ... extra graphics parameters for plot.
##' @return the matrixseries
##' @examples
##' data(vocdia)
##' plot(staticMatrix(vocmatrix$ttarray, digits=4, history=TRUE))
##' @export
plot.matrixseries <- function(x, ...,
                            par=list(), interact=F)
{
  mseries <- x
  op <- par(no.readonly = TRUE); on.exit(par(op))
  par(par)
  #mseries <- staticMatrix(matrix, limit=limit, digits=digits, history=T)
  matrix <- startmatrix(mseries)
  convpoint <- length(mseries)
  mm <- mseries[[convpoint]]
  taillen <- round(convpoint/4)
  limit <- convpoint + taillen
  mseries <- c(mseries,
               lapply(1:taillen,
                      function(y){mm <<- mm %*% matrix}))
  plot(sapply(1:limit, function(x){mseries[[x]][1,1]}), axes=F,
       type='l', ylim=c(0,1), xlab='iteration',
       ylab='amount of speech (steady-state value)',
       ...)
  axis(1)
  axis(2)
  offs = limit/20
  y = limit
  names <- rownames(matrix)
  for (j in 1:ncol(matrix)) {
    y <- y - offs
    print(paste("Plotting convergence for col ",j))
    text(y, mseries[[limit]][1,j]+.02, labels=names[j])
    for (i in 1:nrow(matrix)) {
      if (i == 1 && j == 1)
        next;
      cseries <- sapply(1:limit, function(x){mseries[[x]][i,j]})
      lines(cseries, col=j)
    }
    if (interact)
        locator(1)
  }
  mseries
}

##' Access initital matrix in a \code{matrixseries}
##'
##' Access initital matrix in a \code{matrixseries}
##' @title startmatrix: return the first matrix of a converging series.
##' @param mseries a matrixseries object
##' @return the initial matrix.
##' @examples
##' data(vocdia)
##' x2 <- staticMatrix(vocmatrix$ttarray, digits=4, history=TRUE)
##' ## original matrix
##'  startmatrix(x2)
##' @rdname startmatrix
##' @export startmatrix
startmatrix <- function(mseries) UseMethod('startmatrix')

##' @rdname startmatrix
##' @method startmatrix default
##' @S3method startmatrix default
##' @export startmatrix
startmatrix.default <- function(mseries){
    warning(paste("startmatrix() does not know how to handle object of class ", 
                  class(mseries))) 
}

##' @rdname startmatrix
##' @method startmatrix matrixseries
##' @S3method startmatrix matrixseries
##' @export startmatrix
startmatrix.matrixseries <- function(mseries){
    mseries[[1]]
}

##' Anonymise a vocalisation diagram
##'
##' "anonymise" a \code{vocaldia} turn taking probability matrix by
##' replacing speaker names by variables \eqn{s_1,...,s_n$ s.t. $s_1} is
##' the speaker who spoke the least and \eqn{s_n} the one who did the most
##' talking.
##' @title anonymise: anonymise a vocalisation diagram
##' @param vd a vocalisation diagram (vocaldia object)
##' @return a new vocaldia with speaker names replaced by variables
##'     \eqn{s_1,...,s_n} s.t. \eqn{s_1} is the speaker who spoke the least
##'     and \eqn{s_n} the one who did the most talking.
##' @rdname anonymise
##' @examples
##' data(vocdia)
##' x2 <- getSampledVocalMatrix(subset(atddia, id=='Abbott_Maddock_01'),
##'                             individual=TRUE, nodecolumn='speaker')
##' anonymise(x2)
##' @export anonymise
anonymise <- function(vd) UseMethod('anonymise')

##' @rdname anonymise
##' @method anonymise vocaldia
##' @S3method anonymise vocaldia
anonymise.vocaldia <- function(vd){
  excluded <- c(categories, "Grp", "Floor")
  ## get array of speakers sorted decreasingly by accumulated turn duration 
  ordspk <- sort(vd$tdarray[!names(vd$tdarray) %in% excluded], decreasing=T)
  ## get array indices ordered as ordspk
  idx <- pmatch(names(ordspk), names(vd$tdarray))
  ## get a variable per speaker
  spkvars <- paste(rep("s",length(ordspk)), LETTERS[1:length(ordspk)], sep='')
  ## replace speaker names on accumulated turn duration table
  names(vd$tdarray)[idx] <- spkvars
  ## now do the same for transition turn transiton table row names... 
  idx <- pmatch(names(ordspk), dimnames(vd$ttarray)[[1]])
  dimnames(vd$ttarray)[[1]][idx] <- spkvars
  ## now do the same to col names... 
  idx <- pmatch(names(ordspk), dimnames(vd$ttarray)[[2]])
  dimnames(vd$ttarray)[[2]][idx] <- spkvars
  return(vd)
}

##' @rdname anonymise
##' @method anonymise default
##' @S3method anonymise default
##' @export anonymise
anonymise.default <- function(vd){
     warning(paste("anonymise() does not know how to handle object of class ", 
                  class(vd), '. Try passing a vocaldia.')) 
}



##' Assign types to the pauses (Floor events) in a sequence
##'
##' Identify the pauses in a vector as one of the pauses in
##' \code{pauseTypes}
##' @title identifyPauses: label pauses according to type.
##' @param vocvector a character vector containing a sequence of
##'     vocalisation events
##' @return A vector with all Floor events replaced by the appropriate
##'     pause type identifier.
##' @examples
##' data(vocdia)
##' identifyPauses(atddia$speaker[1:60])
##' @export
identifyPauses <- function(vocvector){
  vocvector <- as.character(vocvector)
  indices <- which(vocvector=='Floor')
  laindex <- length(vocvector)
  ##fvec <- vector(length=length(indices))
  vocvector <- as.character(vocvector)
  for (i in indices){
    if (i == 1 || i == laindex){
      vocvector[i] <- 'Pause'
      next
    }
    if (vocvector[i-1] == 'Grp' || vocvector[i-1] == 'GrpVocalisation'){
      if (vocvector[i+1] == 'Grp' || vocvector[i+1] == 'GrpVocalisation')
        vocvector[i] <- 'GrpPause'
      else
        vocvector[i] <- 'GrpSwitchingPause'
      next
    }
    if (vocvector[i-1] == vocvector[i+1])
      vocvector[i] <- 'Pause'
    else
      vocvector[i] <- 'SwitchingPause'
  }
  vocvector
}

##' Identify switching vocalisations
##' 
##' SwitchingVocalisation is a vocalisation that signals a immediate
##' speaker transition; that is, a transition from speaker to
##' speaker (as opposed to speaker to Grp or speaker to Pause).
##'
##' E.g (speakers A, B, C):
##' \preformatted{
##' AAAAAAAABBBBBBBCCCCCBBBBBPauseBBBBSwitchingPauseAAAAAGrpVocalisation
##'        ^      ^    ^    ^        ^                  ^
##'        |      |    |    |        |                  |
##'        |      |    |    ---------------- Non-SwitchingVocalisation's
##'        |      |    |
##'        ---------------------> SwitchingVocalisation's
##' }
##'
##' @title identifyVocalisations: replace appropriate vocalisation
##'     types
##' @param vocvector a character vector containing a sequence of
##'     vocalisation events
##' @param idswitchvoc if TRUE distinguise between
##'     SwitchingVocalisation and Vocalisation.
##' @return A vector with all events replaced by the appropriate type
##'     identifier.
##' @examples
##' data(vocdia)
##' identifyVocalisations(atddia$speaker[1:60])
##' @export
identifyVocalisations <- function(vocvector, idswitchvoc=T){
  vocvector <- as.character(vocvector)
  vocvector <- identifyGrpVocalisations(vocvector)
  vi <- which(!(vocvector %in% c(pauseTypes,categories,'Floor','Grp')))
  if (idswitchvoc){
    ## find indices of SwitchingVocalisations
    vsi <- vi[which(sapply(1:(length(vi)-1),
                           function(i){
                             vi[i+1]==vi[i]+1 &&
                               vocvector[vi[i]]!=vocvector[vi[i+1]]
                           }))]
    vocvector[vi] <- 'Vocalisation'
    vocvector[vsi] <- 'SwitchingVocalisation'
  }
  else
    vocvector[vi] <- 'Vocalisation'
  vocvector
}


##' Identify group vocalisations
##' 
##' Standardise identifier for group vocalisations
##' @title identifyGrpVocalisations: replace appropriate vocalisation
##'     types
##' @param vocvector a character vector containing a sequence of
##'     vocalisation events
##' @return A vector with all events replaced by the appropriate type
##'     identifier.
##' @examples
##' data(vocdia)
##' identifyGrpVocalisations(atddia$speaker[1:60])
##' @export
identifyGrpVocalisations <- function(vocvector){
  vocvector <- as.character(vocvector)
  grpindices <- which(vocvector=='Grp')
  vocvector[grpindices] <- 'GrpVocalisation'
  vocvector
}

## probs and entropy

##' Conditional (transition ) probability
##'
##' Retrieve \eqn{p(a|b)}, probability of a transition from b to a in an
##' adjacency matrix
##' @title getPofAgivenB: transtion probability.
##' @param a target node
##' @param b source node
##' @param ttarray adjacency matrix
##' @return a transition probability.
##' @export
getPofAgivenB <- function(a, b, ttarray){
  if (! all(c(a,b) %in% names(ttarray[1,]))) ## one of the nodes doesn't exist; return 0
    0
  else
  ttarray[b,a]
}
##' Compute the entropy of a distribution.
##'
##' Compute the entropy of a distribution.
##' @title getEntropy: safely return the Shannon entropy of a distribution.
##' @param distribution a probability distribution.
##' @return a numeric value.
##' @export
getEntropy <- function (distribution){
  PtimesLOG2 <- distribution * log((1 / distribution), 2)
  ## define "0 log 0 = 0" (i.e. get rid of NaN's)
  sum(PtimesLOG2[!is.na(PtimesLOG2)])
}


## external format conversion
##############################

##' Plot a vocalisation diagram
##'
##' Plot a vocalisation diagram
##' @title plot.vocaldia
##' @param x a vocalisation diagram
##' @param package the package to be used for ploting (igraph
##'     (default) or Rgraphviz)
##' @param ... arguments for the layout algorithm
##' @return \code{NULL}
##' @examples
##' data(vocdia)
##' require('igraph')
##' plot(getSampledVocalMatrix(subset(atddia, id=='Abbott_Maddock_01'),
##'                           individual=TRUE, nodecolumn='speaker'))
##' @export
plot.vocaldia <- function(x, ..., package='igraph'){
    vd <- x
    if (requireNamespace("igraph", quietly = TRUE)){
    ##        require('igraph')
    g <- igraph.vocaldia(vd)
    plot(g, layout=igraph::layout.fruchterman.reingold(g), ...)
    return(g)
    }
    ##if (requireNamespace("Rgraphviz", quietly = TRUE)){
    ##    cat('Rgraphviz support under construction. PLease use igraph instead')
    ## if (package=='Rgraphviz'){
    ##     require('Rgraphviz')
    ##     g <- graphNEL.vocaldia(vd)
    ##     w <- as.character(round(unlist(edgeWeights(g)), digits=3))
    ##     w <- w[setdiff(seq(along=w), removedEdges(g))]
    ##     names(w) <- edgeNames(g)
    ##     ea <- at <- list()
    ##     ea$label <- w
    ##     at$edge$fontzise=30
    ##     plot(g, edgeAttrs=ea, attrs=at, ...)
    ##     return(g)
    ##}
    warning(paste('Package ',package, ' not supported. Try igraph.')) 
}
    

##' Create an igraph vocalisation diagram
##'
##' Create a vocalisation diagram
##' @title igraph.vocaldia: Create an igraph vocalisation diagram
##' @param vd a vocalisation diagram
##' @param ... arguments for the layout algorithm
##' @return an igraph
##' @examples
##' data(vocdia)
##' igraph.vocaldia(getSampledVocalMatrix(subset(atddia, id=='Abbott_Maddock_01'),
##'                   individual=TRUE, nodecolumn='speaker'))
##' @export
igraph.vocaldia <- function(vd, ...){
    g <- igraph::graph.adjacency(vd$ttarray, weighted=T)
    igraph::V(g)$label <- names(vd$ttarray[1,])
    igraph::E(g)$label <- round(igraph::E(g)$weight,digits=3)
    igraph::V(g)$size <- 25*exp(vd$tdarray)
    g$layout <- igraph::layout.kamada.kawai(g)
    g
}

## ##' Create a graphNEL vocalisation diagram
## ##'
## ##' Create a vocalisation diagram
## ##' @title graphNEL.vocaldia: Create a graphNEL vocalisation diagram
## ##' @param vd a vocalisation diagram
## ##' @param ... arguments for the layout algorithm
## ##' @return a graphNEL
## ##' @examples
## ##' data(vocdia)
## ##' graphNEL.vocaldia(getSampledVocalMatrix(subset(atddia, id=='Abbott_Maddock_01'),
## ##'                     individual=TRUE, nodecolumn='speaker'))
## ##' @export
## graphNEL.vocaldia <- function(vd, ...){
##     as(unclass(vd$ttarray), 'graphNEL')
## }
## ##

##' Write vocalisation diagram to file in dot (graphviz) notation
##'
##' Write a vocalisation diagram
##' @title write.vocaldia 
##' @param vd a vocalisation diagram
##' @param file name of file to which dot diagram will be written.
##' @param ... arguments passed on to toDotNotation. If "", write to STDOUT.
##' @return \code{NULL}
##' @examples
##' data(vocdia)
##' write.vocaldia(getSampledVocalMatrix(subset(atddia, id=='Abbott_Maddock_01'),
##'                             individual=TRUE, nodecolumn='speaker'),
##'                             file=tempfile(pattern='vocaldia-', fileext='.dot') )
##' @export
write.vocaldia <- function(vd,  file="", ...){
    o <- toDotNotation(vd,  ...)
    o <- paste('## Generated automatically by: ', format(sys.call()), o)
    if (file!="")
        cat("Writing ", file, '\n')
        cat(o, file=file)
}

##' Create vocalisation diagram to file in dot (graphviz) notation
##'
##' Create a vocalisation diagram in dot notation
##' @title toDotNotation: conver vocaldia to graphviz dot notation 
##' @param vd a vocalisation diagram
##' @param individual if TRUE write individual node names
##' @param varsizenode if true set varsizenode in dot
##' @param shape node shape
##' @param fontsize font size
##' @param rankdir direction of ranking (LR, RF etc)
##' @param nodeattribs attributes for node 
##' @param comment comments
##' @return character data containing the diagram in dot format.
##' @examples
##' data(vocdia)
##' toDotNotation(getSampledVocalMatrix(subset(atddia, id=='Abbott_Maddock_01'),
##'                                     individual=TRUE, nodecolumn='speaker'))
##' @seealso graphviz manual
##' @export
toDotNotation <- function(vd, individual=T, varsizenode=T, shape='circle',
                          fontsize=16, rankdir='LR', nodeattribs='fixedsize=true;',
                          comment="")
{
  head <- paste("## diagram generated by vocalgraphs.r\n## ", comment,
                "\ndigraph finite_state_machine {\n",
                'shape=',shape,';',
                'fontsize=', fontsize,';',
                'rankdir=', rankdir, ';',
                nodeattribs)
  links <- ""
  nodes <-  dimnames(vd$ttarray)[[1]]
  for  (i in nodes){
    if (vd$tdarray[i] == 0) next
    width <- log(1000*vd$tdarray[i],base=5)
    width <- if (width < .4 ) .4 else width ## minimal acceptable width
    if (individual){
      nodelabel <- i
    }
    else if (width < .6) {
      nodelabel <- catDescriptions[i]
    }
    else {
      nodelabel <- catDescriptions[i]
    }
    head <- paste(head, "  ", i, 
                  "[",
                  (if (varsizenode) paste("width =", width,', ') else ""),
                  (if (varsizenode)
                   sprintf("label = \"%s \\n%.3f\"",
                           nodelabel,
                           vd$tdarray[i])
                   else
                   paste("label = \"",nodelabel,"\", ")
                   ),
                   if (width < .6) "fontsize=8", 
                  "];\n")
    for (j in nodes) {
      if (vd$ttarray[i,j] == 0) next 
      links <- paste(links, "  ", i, "->", j,
                     "[ label =", sprintf("%.3f", vd$ttarray[i,j]),
                     "];\n")
    }
  }
  o <- paste(head, links, "}\n")
  return(o);
}
