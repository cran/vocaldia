% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/markov.R
\name{toDotNotation}
\alias{toDotNotation}
\title{toDotNotation: conver vocaldia to graphviz dot notation}
\usage{
toDotNotation(
  vd,
  individual = T,
  varsizenode = T,
  shape = "circle",
  fontsize = 16,
  rankdir = "LR",
  nodeattribs = "fixedsize=true;",
  comment = ""
)
}
\arguments{
\item{vd}{a vocalisation diagram}

\item{individual}{if TRUE write individual node names}

\item{varsizenode}{if true set varsizenode in dot}

\item{shape}{node shape}

\item{fontsize}{font size}

\item{rankdir}{direction of ranking (LR, RF etc)}

\item{nodeattribs}{attributes for node}

\item{comment}{comments}
}
\value{
character data containing the diagram in dot format.
}
\description{
Create vocalisation diagram to file in dot (graphviz) notation
}
\details{
Create a vocalisation diagram in dot notation
}
\examples{
data(vocdia)
toDotNotation(getSampledVocalMatrix(subset(atddia,
                                           id=='Abbott_Maddock_01'),
                             individual=TRUE, nodecolumn='speaker'))
}
\seealso{
graphviz manual
}
