% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vocalgraphs.R
\name{getTurnTakingMatrix}
\alias{getTurnTakingMatrix}
\title{getSampledVocalCountMatrix: generate vocalisation diagrams}
\usage{
getTurnTakingMatrix(
  df,
  begin = "begin",
  end = "end",
  nodecolumn = "role",
  individual = FALSE,
  noPauseTypes = FALSE
)
}
\arguments{
\item{df}{a data frame consisting, minimally, of a column for
vocalisation/pause start times, a column for end times, and a
column identifying the speaker, speaker role or 'Floor' (for
silences).}

\item{begin}{the name of the column containing the start time of
the vocalisation event in a row.}

\item{end}{the name of the column containing the end time of the
vocalisation event in the same row.}

\item{nodecolumn}{the name of the column containing the node
(speaker) name (e.g. 'speaker', 'role').}

\item{individual}{whether to include individual speakers or group
them into a single Vocalisation node}

\item{noPauseTypes}{if TRUE, ignore distinctions between pauses
(SwitchingPause, GrpSwitchingPause, etc)}
}
\value{
a vocaldia object, consisting of a vocalisation matrix
    (vocmatrix) where cell <m,n> contains the counts of
    transitions from node n to node m, and a table of absolute
    durations of vocalisation events.
}
\description{
Generate a vocalisation diagram with absolute vocalisation durations.
}
\details{
A vocalisation diagram (vocaldia) is a representation of a
dialogue as a Markov process whose cell <m,n> contains the
transition probability from node n to node m). Unlike
\code{\link{getSampledVocalCountMatrix}} this function
accummulates event durations directly, therefore resulting in no
self-transitions (in general).
}
\examples{
x <- subset(atddia, id=='Abbott_Maddock_01')
getTurnTakingMatrix(x)
getTurnTakingMatrix(x, individual=TRUE)
}
\references{
S. Luz. Automatic identification of experts and performance
    prediction in the multimodal math data corpus through analysis
    of speech interaction. In \emph{Proceedings of the 15th ACM on
    International Conference on Multimodal Interaction, ICMI'13},
    pages 575--582, New York, NY, USA, 2013. ACM.
}
\seealso{
(Luz, 2013) and \code{\link{getTurnTakingMatrix}}.
}
