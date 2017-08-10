## VOCALGRAPHS.R:
## Functions for creating adjacency matrices of vocalisation graphs from dataframes.
## 
### This program is free software; you can redistribute it and/or
### modify it under the terms of the GNU General Public License
### as published by the Free Software Foundation; either version 2
### of the License, or (at your option) any later version.
###
### (c) 2017 S. Luz (luzs@acm.org)
##source('markov.R')

##' @docType package
##' @bibliography /home/luzs/lib/tex/bib/luz.bib

##' Create and manipulate vocalisation matrices and diagrams 
##'
##' \code{vocaldia} provides a number of functions for generating
##' vocalisation diagrams (vocaldias) from data frames containing, minimally, a
##' column for start time of a vocalisation event (speech, silence,
##' group-talk etc), a column for end time, and a column for the event
##' identifier.
##'
##' Functions \code{\link{getSampledVocalMatrix}} and
##' \code{\link{getTurnTakingProbMatrix}} generate alternative
##' versions of adjacency matrices for
##' vocaldias. \code{\link{staticMatrix}} generates steady state
##' diagrams from a vocaldia. \code{\link{printARFFfile}} generates a
##' 'flat' representation of vocaldias for classifier training and
##' evaluation.
##'
##' @author Saturnino Luz \email{luzs@@acm.org}
##' @references
##'   S. Luz. Automatic identification of experts and performance
##'   prediction in the multimodal math data corpus through analysis
##'   of speech interaction. In \emph{Proceedings of the 15th ACM on
##'   International Conference on Multimodal Interaction, ICMI'13},
##'   pages 575--582, New York, NY, USA, 2013. ACM.
##'
##'   S. Luz. The non-verbal structure of patient case discussions in
##'   multidisciplinary medical team meetings. \emph{ACM Transactions on
##'   Information Systems}, 30(3):17:1--17:24, 2012
##'
##'   Dabbs, J. M. J. and Ruback, B. Dimensions of group process: Amount and
##'   structure of vocal interaction. \emph{Advances in Experimental Social
##'   Psychology} 20, 123-169, 1987.
##'
##'   Jaffe , J. and Feldstein, S. Rhythms of
##'   dialogue. ser. \emph{Personality and Psychopathology}. Academic
##'   Press, New York, 1976.
##' @import graphics
"_PACKAGE"

##' A sample vocalisation matrix
##'
##' A \code{vocaldia} object containing a 3-speaker dialogue 
##'
##' @format A list containing 2 arrayes
##' \describe{
##'   \item{ttarray}{The vocaldia adjacency matrix}
##'   \item{tdarray}{The proportional durations (stationary probabilities) of each event (node)}
##' }
##' @source This dataset was generated from the Multomodal Leatrning
##'     Analytics dataset, for the eponymous ICMI'13 Grand
##'     Challenge. The use these vocaldias were put to is described in
##'     Luz (2013). The full dataset and code is available
##'     at https://gitlab.scss.tcd.ie/saturnino.luz/icmi-mla-challenge
##' @references
##'   S. Luz. Automatic identification of experts and performance
##'   prediction in the multimodal math data corpus through analysis
##'   of speech interaction. In \emph{Proceedings of the 15th ACM on
##'   International Conference on Multimodal Interaction, ICMI'13},
##'   pages 575--582, New York, NY, USA, 2013. ACM.
"vocmatrix"

##' A sample Medical Team Meeting dialogue encoded as a vocaldia
##'
##' A dataset containing 38 dialogues (17 control patients, and 21 AD
##' patients) and 7869 vocalisation events.
##'
##' @format A data frame with 7869 rows and 7 variables:
##' \describe{
##'   \item{id}{The dialogue indentifier}
##'   \item{begin}{The start time of a speech turn or silence interval}
##'   \item{end}{The end time of a speech turn or silence interval}
##'   \item{speaker}{An identifier for the speaker of the turn, or Floor for silence.}
##'   \item{role}{The speaker's role (patient, interviewer, other, or Floor}
##'   \item{trans}{The transcription of the turn (blanked out for anonymity)}
##'   \item{dx}{The diagnosis (ad or nonad}
##'  }
##' @source This dataset was generated from the Carolina Conversations
##'     Collection, and used in the work described in De La Fuente,
##'     Albert and Luz:
##'     "Detecting cognitive decline through dialogue processing",
##'     2017. For the full data set, please contact the Medical
##'     University of South Carolina (MUSC)
##'     http://carolinaconversations.musc.edu/
"atddia"

##' Generate a probabilistic vocalisation diagram through 'sampling'.
##'
##' A vocalisation diagram (vocaldia) is a representation of a
##' dialogue as a Markov process whose cell <m,n> contains the
##' transition probability from node n to node m). 
##' @title getSampledVocalCountMatrix: generate vocalisation diagrams
##' @param df a data frame consisting, minimally, of a column for
##'     vocalisation/pause start times, a column for end times, and a
##'     column identifying the speaker, speaker role or 'Floor' (for
##'     silences).
##' @param ... general parameter to be passed to
##'     \code{\link{getSampledVocalCountMatrix}}
##' @return a vocaldia object, consisting of a vocalisation matrix
##'     (vocmatrix) where cell <m,n> contains the transition
##'     probability from node n to node m, and a table of prior
##'     probabilities (stationary distribution) per node.
##' @author Saturnino Luz  \email{luzs@@acm.org}
##' @references
##'   S. Luz. Automatic identification of experts and performance
##'   prediction in the multimodal math data corpus through analysis
##'   of speech interaction. In \emph{Proceedings of the 15th ACM on
##'   International Conference on Multimodal Interaction, ICMI'13},
##'   pages 575--582, New York, NY, USA, 2013. ACM.
##' @examples
##' data(vocdia) 
##' getSampledVocalMatrix(subset(atddia,
##'          id=='Abbott_Maddock_01'),nodecolumn='speaker', individual=TRUE)
##' @seealso \code{\link{getSampledVocalCountMatrix}}
##' @export
getSampledVocalMatrix <- function (df, ...) {
  vcm <- getSampledVocalCountMatrix(df, ...)
  tdarray <- vcm$tdarray/sum(vcm$tdarray)   ## normalise turn-distribution vector
  rsum <-  apply(vcm$ttarray, 1, sum)                   ## sum rows
  ttarray <- apply(vcm$ttarray, 2, function(z){z/rsum}) ## normalise each column

  class(ttarray) <- 'vocalmatrix'
  class(tdarray) <- 'vocduration'
  vd <- list(ttarray=ttarray, tdarray=tdarray)
  class(vd) <- 'vocaldia'
  
  return(vd)
}

##' Generate a count vocalisation diagram through 'sampling'.
##'
##' A vocalisation diagram (vocaldia) is a representation of a
##' dialogue as a Markov process whose cell <m,n> contains the
##' transition probability from node n to node m). This function for
##' 'cases' (an identifier for a case or a vector of identifiers
##' identifying a set of cases) in data frame 'df', obtained by
##' sampling the timeline every 'rate'-th second (see
##' getSampledVocalCountMatrix).
##' @title getSampledVocalCountMatrix: generate vocalisation diagrams
##' @param cdf a data frame consisting, minimally, of a column for
##'     vocalisation/pause start times, a column for end times, and a
##'     column identifying the speaker, speaker role or 'Floor' (for
##'     silences).
##' @param rate the rate at which to sample the vocalisation events
##'     (in seconds)
##' @param individual whether to include individual speakers or group
##'     them into a single Vocalisation node
##' @param noPauseTypes if TRUE, ignore distinctions between pauses
##'     (SwitchingPause, GrpSwitchingPause, etc)
##' @param begin the name of the column containing the start time of
##'     the vocalisation event in a row.
##' @param end the name of the column containing the end time of the
##'     vocalisation event in the same row.
##' @param nodecolumn the name of the column containing the node
##'     (speaker) name (e.g. 'speaker', 'role').
##' @return a vocaldia object, consisting of a vocalisation matrix
##'     (vocmatrix) where cell <m,n> contains the counts of
##'     transitions from node n to node m, and a table of prior
##'     probabilities (stationary distribution) per node.
##' @seealso (Luz, 2013)
##' @examples
##' data(vocdia) 
##' getSampledVocalCountMatrix(subset(atddia,
##'      id=='Abbott_Maddock_01'), nodecolumn='role')
##' @export
getSampledVocalCountMatrix <- function (cdf, rate=1, individual=FALSE, noPauseTypes=FALSE,
                                        begin='begin', end='end',
                                        nodecolumn='role') {
    ##cdf <- subset(df, id %in% ids)
    if(!individual)
        cdf[[nodecolumn]] <- as.factor(identifyVocalisations(cdf[[nodecolumn]]))
    if(!noPauseTypes)
        cdf[[nodecolumn]] <- as.factor(identifyPauses(cdf[[nodecolumn]]))
    ## no need for df from this point on; one could split this into 2
    ## separate functions

    btime <- cdf[[begin]][1]
    etime <- cdf[[end]][length(cdf[[end]])]
    nodes <- as.character(unique(cdf[[nodecolumn]]))
    ## absolute duration (static distribution); ttarray^n (i.e. the
    ## matrix representation of the Chapman-Kolmogorov equation) should
    ## converge to tdarray as n approaches infinity
    tdarray <- array(data=0,
                     dim=c(length(nodes)),
                     dimnames=list(nodes))
    tdarray <-  sapply(nodes,function(y){sum(cdf[[end]][cdf[[nodecolumn]]==y]
                                             -cdf[[begin]][cdf[[nodecolumn]]==y])})
    ttarray <- array(data=0,
                     dim=c(length(nodes),length(nodes)),
                     dimnames=list(nodes, nodes))
    pspk <- character(length=0)
    for (t in seq(btime, etime, by=rate) ) {
        sp <- as.character(cdf[[nodecolumn]][cdf[[begin]] <= t & t < cdf[[end]]])
        if ( length(sp) > 1 ) {
            warning(paste("Invalid overlapping vocalisation at ", t))
            sp <- sp[1]
        }
        if ( length(sp) == 0 || length(pspk) == 0) {
            pspk <- sp
            next
        }
        ttarray[pspk, sp] <- ttarray[pspk, sp] + 1
        pspk <- sp
    }
    ## prevent zero-transition rows (and div by zero in getSampledVocalMatrix)
    if (sum(ttarray[pspk,])==0) 
        ttarray[pspk,pspk] <- 1
    return(list(ttarray=ttarray, tdarray=tdarray))
}

##' Identify turn types
##'
##' Return one of {Vocalisation, GrpVocalisation, ...} or identifier.
##' @title getTurnType: return type of turn
##' @param df  a data frame consisting, minimally, of a column for
##'     vocalisation/pause start times, a column for end times, and a
##'     column identifying the speaker, speaker role or 'Floor' (for
##'     silences).
##' @param i the identifier (index number) whose type will be returned
##' @param individual if TRUE, return the identifier, a Pause or Grp
##' @param nodecolumn the name of the column containing the node
##'     (speaker) name (e.g. 'speaker', 'role').
##' @param noPauseTypes if TRUE, ignore distinctions between pauses
##'     (SwitchingPause, GrpSwitchingPause, etc)
##' @return a string containing the turn type or identifier.
##' @examples
##' data(vocdia)
##' atddia[1:10,]
##' getTurnType(atddia, 3, nodecolumn='role') ## a vocalisation
##' getTurnType(atddia, 4, nodecolumn='role') ## a pause
##' @export
getTurnType <- function(df, i, individual=FALSE, nodecolumn='speaker',
                        noPauseTypes=F) {
  ##paste(i,", ")
  currentspeaker <- as.character(df[[nodecolumn]][i])
  nextspeaker <- 'none'
  if (length(df$id) > i && df$id[i] == df$id[i+1] ){
    nextspeaker <- as.character(df[[nodecolumn]][i+1])
  }
  prevspeaker <- if (i > 1 && df$id[i] == df$id[i-1]) as.character(df[[nodecolumn]][i-1]) else 'none'
  if (currentspeaker == 'Floor') {
    if (noPauseTypes) return('Pause') else return(getPauseType(prevspeaker, nextspeaker))
  }
  else if (currentspeaker == 'Grp'){
    if (individual) return(currentspeaker) else return('GrpVocalisation') 
  }
  else {
    if (individual) return(currentspeaker) else return('Vocalisation')
  }
}

##' Identify the type of pause between vocalisations.
##'
##' The type of pause a 'Floor' (silence) event represents can be:
##' 'Pause', 'SwitchingPause', 'GrpPause', or 'GrpSwitchingPause'. See
##' (Luz, 2013) for details.
##' @title getPauseType: name pause type between two vocalisation events.
##' @param prevspeaker speaker of the vocalisation immediately before Floor
##' @param nextspeaker speaker of the vocalisation immediately after Floor
##' @return the pause type.
##' @seealso \code{\link{namePauses}}
##' @examples
##' getPauseType('a', 'b')
##'  ## [1]  "SwitchingPause"
##' getPauseType('a', 'Grp')
##'  ## [1]  "SwitchingPause"
##' getPauseType('Grp', 'Grp')
##'  ## [1]  "GrpPause"
##' getPauseType('Grp', 'a')
##'  ## [1]  "GrpSwitchingPause"
##' getPauseType('a', 'a')
##'  ##[1] "Pause"
##' @export
getPauseType <- function(prevspeaker, nextspeaker){
  if (prevspeaker == 'Floor' ){
    warning("OOPS! Two consecutive Floor turns found. A coding error??")
    return('Pause')
  }
  if ( prevspeaker == 'none' || nextspeaker == 'none') return('Pause')
  if ( prevspeaker == 'Grp'){
    if (nextspeaker == 'Grp') return('GrpPause') else return('GrpSwitchingPause')
  }
  if (prevspeaker != nextspeaker) return('SwitchingPause') else return ('Pause')
}

##' Replace identified pause pause types in data frame.
##'
##' replace all 'Floor' speakers in df by 'Pause', 'SwitchingPause'
##' etc, and return a new data fame containing pause types in place of
##' 'Floor' (see markov.R, identifyPauses() for a better
##' implementation)
##' @title namePauses: name pause types.
##' @param df a data frame consisting, minimally, of a column for
##'     vocalisation/pause start times, a column for end times, and a
##'     column identifying the speaker, speaker role or 'Floor' (for
##'     silences).
##' @param nodecolumn the name of the column containing the node
##'     (speaker) name (e.g. 'speaker', 'role').
##' @return a data.frame with pauses in nodecolumn replaced by different pause types.
##' @seealso \code{\link{identifyPauses}} for a better implementation
##' @examples
##' data(vocdia)
##' x <- subset(atddia, id=='Abbott_Maddock_01')
##' x[1:15,1:6]
##' namePauses(x)[1:15,1:6]
##' @export
namePauses <- function(df, nodecolumn='role') {
  nspeakers <- sapply(as.integer(row.names(df)),
                      function (X) {getTurnType(df, X, individual=T)})
  df[[nodecolumn]] <- nspeakers
  return(df)
}


## non-sampled versions of vocalisation graphs (no-self transitions allowed)
###########################################################################

##' Convert a data frame into a vocalisation diagram using counts rather than sampling.
##'
##' Unlike \code{\link{getSampledVocalMatrix}}, this function is based
##' on transition counts rather than sampled intervals. As a result,
##' where in this version self transitions will always be set to 0
##' (since a vocalisation by a speaker is never followed by another
##' vocalisation by the same speaker) in the sampled version self
##' transitons will usually dominate the distribution, since the
##' speaker who is speaking now is very likely to be the one who were
##' speaking one second ago.
##' @title getTurnTakingProbMatrix: create a vocaldia from a
##'     data.frame.
##' @param df a data frame consisting, minimally, of a column for
##'     vocalisation/pause start times, a column for end times, and a
##'     column identifying the speaker, speaker role or 'Floor' (for
##'     silences).
##' @param individual whether to include individual speakers or group
##'     them into a single Vocalisation node
##' @param ... other parameters to be passed to
##'     \code{\link{getTurnTakingMatrix}}.
##' @return a vocaldia object, consisting of a vocalisation matrix
##'     (vocmatrix) where cell \eqn{(m,n)} contains the probabilities \eqn{P(n|m)}
##'     transitions to node \eqn{n} from node \eqn{m}, and a table of prior
##'     probabilities (stationary distribution) per node.
##' @seealso (Luz, 2013) and  \code{\link{getTurnTakingMatrix}}.
##'
##'   S. Luz. Automatic identification of experts and performance
##'   prediction in the multimodal math data corpus through analysis
##'   of speech interaction. In \emph{Proceedings of the 15th ACM on
##'   International Conference on Multimodal Interaction, ICMI'13},
##'   pages 575--582, New York, NY, USA, 2013. ACM.
##' @examples
##' x <- subset(atddia, id=='Abbott_Maddock_01')
##' getTurnTakingProbMatrix(x)
##' getTurnTakingProbMatrix(x, individual=TRUE)
##' @export
getTurnTakingProbMatrix <- function(df, individual=FALSE, ...)
{
  t <- getTurnTakingMatrix(df, individual=individual, ...)
  if (individual){
    nodes <-  dimnames(t$ttarray)[[1]]
  } else {
    nodes <- categories
  }
  for (i in nodes){
    s <- sum(t$ttarray[i,])
    t$ttarray[i,] <- if (s > 0) t$ttarray[i,]/sum(t$ttarray[i,]) else 0
  }
  t$tdarray <- t$tdarray/sum(t$tdarray)
  class(t$ttarray) <- 'vocalmatrix'
  class(t$tdarray) <- 'vocduration'
  class(t) <- 'vocaldia'
  return(t)
}

##' Generate a vocalisation diagram with absolute vocalisation durations.
##'
##' A vocalisation diagram (vocaldia) is a representation of a
##' dialogue as a Markov process whose cell <m,n> contains the
##' transition probability from node n to node m). Unlike
##' \code{\link{getSampledVocalCountMatrix}} this function
##' accummulates event durations directly, therefore resulting in no
##' self-transitions (in general).
##' @title getSampledVocalCountMatrix: generate vocalisation diagrams
##' @param df a data frame consisting, minimally, of a column for
##'     vocalisation/pause start times, a column for end times, and a
##'     column identifying the speaker, speaker role or 'Floor' (for
##'     silences).
##' @param begin the name of the column containing the start time of
##'     the vocalisation event in a row.
##' @param end the name of the column containing the end time of the
##'     vocalisation event in the same row.
##' @param nodecolumn the name of the column containing the node
##'     (speaker) name (e.g. 'speaker', 'role').
##' @param individual whether to include individual speakers or group
##'     them into a single Vocalisation node
##' @param noPauseTypes if TRUE, ignore distinctions between pauses
##'     (SwitchingPause, GrpSwitchingPause, etc)
##' @return a vocaldia object, consisting of a vocalisation matrix
##'     (vocmatrix) where cell <m,n> contains the counts of
##'     transitions from node n to node m, and a table of absolute
##'     durations of vocalisation events.
##' @seealso (Luz, 2013) and \code{\link{getTurnTakingMatrix}}.
##' @references
##'   S. Luz. Automatic identification of experts and performance
##'     prediction in the multimodal math data corpus through analysis
##'     of speech interaction. In \emph{Proceedings of the 15th ACM on
##'     International Conference on Multimodal Interaction, ICMI'13},
##'     pages 575--582, New York, NY, USA, 2013. ACM.
##' @examples
##' x <- subset(atddia, id=='Abbott_Maddock_01')
##' getTurnTakingMatrix(x)
##' getTurnTakingMatrix(x, individual=TRUE)
##' @export
getTurnTakingMatrix <- function(df, begin='begin', end='end',
                                nodecolumn='role', 
                                individual=FALSE, noPauseTypes=FALSE)
{
    ##df <- subset(df, df$id %in% ids)
  if (individual){
      nodes <- as.character(unique(df[[nodecolumn]]))
      nodes <- nodes[nodes != 'Floor']
      if (!noPauseTypes){
          nodes <- c(nodes, pauseTypes)
      }
      else {
          nodes <- c(nodes, 'Floor')
      }
  } else {
      nodes <- categories
  }
  ## ttarray (turn transition array) is a category x category matrix
  ## s.t. cell (i,j) stores the number of transitions from turn
  ## category i to j
  ttarray <- array(data=0,
                   dim=c(length(nodes),length(nodes)),
                   dimnames=list(nodes, nodes))
  ## tdarray (turn duration array) stores the total duration of each turn category
  tdarray <- array(data=0,
                   dim=c(length(nodes)),
                   dimnames=list(nodes))
  prevttype <- 'none'
  prevcase <- 'none'
  for (i in 1:length(df[[nodecolumn]])) {
    ##if ( !df$id[i] %in% ids ) next;
      ttype <- getTurnType(df, i, individual=individual,
                           nodecolumn=nodecolumn, noPauseTypes=noPauseTypes)
    tdarray[ttype] <-  tdarray[ttype] + (df$end[i]-df$begin[i]);
    if (prevcase == df$id[i]) { ##  no transition is recorded across cases
      ttarray[prevttype,ttype] <- ttarray[prevttype,ttype] + 1
    }
    else {
      prevcase <- df$id[i]
    }
    prevttype <- ttype
  }
  return(list(ttarray=ttarray, tdarray=tdarray))
}

##' Generate ARFF files from vocalisation diagrams
##'
##' Use this function to generate turn-taking diragrams in ARFF format for
##  processing with, for instance, the WEKA machine learning toolkit.
##' @title printARFFfile: Create arff files by creating and flattening vocaldias
##' @param df df a data frame consisting, minimally, of a column for
##'     vocalisation/pause start times, a column for end times, and a
##'     column identifying the speaker, speaker role or 'Floor' (for
##'     silences).
##' @param ids Ids of dialogues to generate (as defined in column named idcolumn)
##' @param idcolumn the name of the column containing the dialogue id
##' @param noPauseTypes if TRUE, ignore distinctions between pauses
##'     (SwitchingPause, GrpSwitchingPause, etc)
##' @param sampled if >0 use \code{\link{getSampledVocalMatrix}} with rate=sampled
##' @param individual whether to include individual speakers or group
##'     them into a single Vocalisation node
##' @param nodecolumn the name of the column containing the node
##'     (speaker) name (e.g. 'speaker', 'role').
##' @param classcolumn the name of the column containing the target class (or value).
##' @param file name of ARFF file to be generated, or "" (print to console).
##' @return NULL
##' @seealso
##'     \code{\link{getSampledVocalCountMatrix}},
##'     \code{\link{getTurnTakingProbMatrix}}.
##' @references
##'   S. Luz. Automatic identification of experts and performance
##'   prediction in the multimodal math data corpus through analysis
##'   of speech interaction. In \emph{Proceedings of the 15th ACM on
##'   International Conference on Multimodal Interaction, ICMI'13},
##'   pages 575--582, New York, NY, USA, 2013. ACM.
##' @examples
##' data(vocdia)
##' atdarff <- tempfile(pattern='vocaldia-', fileext='arff')
##' printARFFfile(atddia, individual=TRUE, classcolumn='dx',
##'               file=atdarff, noPauseTypes=FALSE)
##' library("foreign")
##' x1 <- read.arff(atdarff)
##' x1[1:3,]
##' ## remove empty columns
##' x1[,c(unlist(apply(x1[1:(ncol(x1)-1)],2,sum)!=0), TRUE)]
##' @export
printARFFfile <- function(df,
                          ids=c(),
                          idcolumn='id',
                          noPauseTypes=F,
                          sampled=0,
                          individual=TRUE,
                          nodecolumn="role",
                          classcolumn='dx',
                          file="")
{
  call <- format(sys.call())
  ## exclude dubious annotation
  #print(formals())
  head <- paste(c("% file automatically generated by vocalgraphs.R\n",
                  "% ", call,
                  "\n\n@RELATION mdtm\n\n"), collapse="")
  if (length(ids) == 0) {
     ids <- levels(df$id)
  }
  nodes <- c()
  tvector <- c(list())
  vi <- 0
  for  (i in ids) {
      if (sampled>0)
          t <- getSampledVocalMatrix(df[df[[idcolumn]]==i,],
                                     individual=individual,
                                     noPauseTypes=noPauseTypes,
                                     nodecolumn=nodecolumn, rate=sampled)
      else
          t <- getTurnTakingProbMatrix(df[df[[idcolumn]]==i,],
                                       nodecolumn=nodecolumn,
                                       individual=individual,
                                       noPauseTypes=noPauseTypes)
      nodes <- union(nodes, names(t$tdarray))
      tvector <- c(tvector, list(list(ttpm=t, id=i)))
  }
  ## attributes consits of all node names (the values of which
  ## represent the probability of that node) and all combinations of
  ## nodes (values representing transition probabilities)
  attributes <- c(paste(nodes),
                  sapply(nodes, function(x){paste(nodes,'-',x,sep='')}))
  targetclass <- unique(df[[classcolumn]])
  mla <- max(nchar(attributes))
  body <- c(sapply(attributes,
                   function(x){pad <- mla - nchar(x);
                       paste(c("@ATTRIBUTE ",x,rep(" ",pad)," REAL\n"),collapse="")}),
            paste(c("@ATTRIBUTE ", classcolumn ," {",
                    paste(targetclass,collapse=","),"}\n"),collapse=""),
            "\n\n@DATA\n")
  ##
  for (tv in tvector) {
      t <- tv$ttpm
      i <- tv$id
      ## tda: sparse array containing all spkr turn durations,
      ## whereas t$tdarray contain only those that occur in case i
      tda <- array(data=0, dim=c(length(nodes)), dimnames=list(nodes))
      tda[names(t$tdarray)] <- t$tdarray
      tta <- array(data=0, dim=c(length(nodes),length(nodes)),
                   dimnames=list(nodes,nodes))
      line <- paste(tda, collapse=",")
      n <- dimnames(t$ttarray)[[1]]
    for  (j in n){
      for  (k in n){
        if (t$ttarray[j,k] == 0) next
        tta[j,k] <- t$ttarray[j,k]
      }
    }
    line <- paste(c(line,
                    paste(tta, collapse=","),
                    as.character(df[[classcolumn]][df$id==i][1])),
                  collapse=",")
    body <- c(body, line,"\n")
  }
  o <- paste(c(head,body),collapse="")
  cat(o, file=file)
}

## loadDataset <- function(file='~/lib/projects/DTP/HRI-Dementia/data/ccc/dialog_ALL.csv'){
##     d <- read.csv(file, as.is=T)
##     ## fill in missing roles (possibly inconsistently; check)
##     d$role[d$role==''] <- 'Other'
##     d$begin <- as.numeric(d$begin)
##     d$end <- as.numeric(d$end)
##     numsil <- sum(sapply(1:(nrow(d)-1), function(X)d$end[X]<d$begin[X+1]))
##     dout <- matrix(nrow=nrow(d)+numsil-1, ncol=ncol(d))
##     colnames(dout) <- names(d)
##     j <- 0
##     for (i in 1:(nrow(d)-1)){
##         dout[j <- j+1,] <- as.character(d[i,])
##         if (d$end[i] < d$begin[i+1]) {
##             dout[j <- j+1,] <- c(d$id[i], d$end[i], d$begin[i+1],
##                                  'Floor', 'Floor', d$dx[i], '___')
##         }
##     }
##     daux <- data.frame(dout[,1], as.numeric(dout[,2]), as.numeric(dout[,3]),
##                        dout[,4], dout[,5], dout[,6])
##     ##daux <- data.frame(dout[,1:(ncol(d)-1)], stringsAsFactors=F)
##     ## keep transcriptions as strings rather than converting them to factors
##     dout <- cbind(daux,trans=dout[,ncol(d)], stringsAsFactors=F)
##     names(dout) <- colnames(d)
##     ## annotated silences (inconsistent; check manually)
##     sil <- c(grep('^[^a-zA-Z]+$', dout$trans), which(dout$speaker==''))
##     dout$role[sil] <- dout$speaker[sil] <- 'Floor'
##     dout
## }
