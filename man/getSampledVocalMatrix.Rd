% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vocalgraphs.R
\name{getSampledVocalMatrix}
\alias{getSampledVocalMatrix}
\title{getSampledVocalCountMatrix: generate vocalisation diagrams}
\usage{
getSampledVocalMatrix(df, ...)
}
\arguments{
\item{df}{a data frame consisting, minimally, of a column for
vocalisation/pause start times, a column for end times, and a
column identifying the speaker, speaker role or 'Floor' (for
silences).}

\item{...}{general parameter to be passed to
\code{\link{getSampledVocalCountMatrix}}}
}
\value{
a vocaldia object, consisting of a vocalisation matrix
    (vocmatrix) where cell <m,n> contains the transition
    probability from node n to node m, and a table of prior
    probabilities (stationary distribution) per node.
}
\description{
Generate a probabilistic vocalisation diagram through 'sampling'.
}
\details{
A vocalisation diagram (vocaldia) is a representation of a
dialogue as a Markov process whose cell <m,n> contains the
transition probability from node n to node m).
}
\examples{
data(vocdia) 
getSampledVocalMatrix(subset(atddia,
         id=='Abbott_Maddock_01'),nodecolumn='speaker', individual=TRUE)
}
\references{
S. Luz. Automatic identification of experts and performance
  prediction in the multimodal math data corpus through analysis
  of speech interaction. In \emph{Proceedings of the 15th ACM on
  International Conference on Multimodal Interaction, ICMI'13},
  pages 575--582, New York, NY, USA, 2013. ACM.
}
\seealso{
\code{\link{getSampledVocalCountMatrix}}
}
\author{
Saturnino Luz  \email{luzs@acm.org}
}
