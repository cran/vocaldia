# vocaldia

<!-- badges: start -->
<!-- badges: end -->

This package contains functions that create and manipulate vocalisation
diagrams. Vocalisation diagrams date back to early work in psychiatry
(Jaffe and Feldstein, 1970) and social psychology (Dabbs and Ruback,
1987) but have only recently been employed as a data representation
method for machine learning (Luz, 2013; Luz and Kane, 2009).

This provides a number of functions for generating vocalisation
diagrams (vocaldias) from data frames containing, minimally, a column
for start time of a vocalisation event (speech, silence, group-talk
etc), a column for end time, and a column for the event identifier. It
also contains some basic functions for reading and processing files
from DementiaBank (.cha transcripts and audio files). 

Functions `getSampledVocalMatrix` and `getTurnTakingProbMatrix`
generate alternative versions of adjacency matrices for
vocaldias. `staticMatrix` generates steady state diagrams from a
vocaldia. `printARFFfile` generates a 'flat' representation of
vocaldias for classifier training and evaluation.

## Installation

You can install the released version of vocaldia from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("vocaldia")
```

## Example

The following examples illustrate the use of vocaldia to create and
visualise vocalisation graphs and their properties. 


``` r
library(vocaldia)
## load some data
data(vocdia)

## select a dialogue
x <- subset(atddia, id=='Abbott_Maddock_01')

## show a probability matrix 
getTurnTakingProbMatrix(x)

## if you have igraph installed, visualise a vocal matrix
require('igraph')
subset(atddia, id=='Abbott_Maddock_01') %>% 
    getSampledVocalMatrix(individual=TRUE, nodecolumn='speaker') 
    %>% igraph.vocaldia %>% plot

## plot steady state of the Markov diagram
plot(staticMatrix(vocmatrix$ttarray, digits=4, history=TRUE))
```

See the following publication for further examples of use of this
package:

Luz S, De La Fuente Garcia S, Albert P. A Method for Analysis of
Patient Speech in Dialogue for Dementia Detection. In Resources and
ProcessIng of linguistic, para-linguistic and extra-linguistic Data
from people with various forms of cognitive impairment. Paris, France:
ELRA. 2018. p. 35-42 (https://arxiv.org/abs/1811.09919)
