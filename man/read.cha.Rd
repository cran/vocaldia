% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dbank.R
\name{read.cha}
\alias{read.cha}
\title{read.cha read CHA transcription file (format used by DementiaBank)}
\usage{
read.cha(file, sildir = NULL, silsuffix = "c.mp3.csv")
}
\arguments{
\item{file}{.cha file to reas}

\item{sildir}{silences directory}

\item{silsuffix}{silence files suffix}
}
\value{
a list containing the PID, a dataframe containing the speaker IDs and demogrephics, and a dataframe containing the speaker IDs, transcribed utterances, start and en times, speech rates etc.
}
\description{
read.cha: read CHA transcription file (format used by DementiaBank)
}
\author{
luzs
}
