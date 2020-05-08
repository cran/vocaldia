## dbank.R:
## Functions for reading and processing DementiaBank files
## 
### This program is free software; you can redistribute it and/or
### modify it under the terms of the GNU General Public License
### as published by the Free Software Foundation; either version 2
### of the License, or (at your option) any later version.
###
### (c) 2019 S. Luz (luzs@acm.org)

## 
##' Build a data frame createwith vocalisation statistics
##'
##' @title makeVocalStatsDataset: create a dataset of vocalisation statistics (1 row per patient) 
##' @param dir a string or vector containing the location (directory path) of the DementiaBank transcript files (.cha files)
##' @param sildir directory where silence csv files are stored 
##' @param silsuffix the suffix of the silence profile files 'c.mp3.csv'. The format of such files should be the format used by Audacity label files, i.e. 'start time, end time, label' (without header), where 'label' should be 'silence' 
##' @param srdir  directory where speech rate csv (1 value per utterance) files are stored 
##' @param srsuffix the suffix of the speech rate files (default: sre)
##' @param sprate compute speech rate? (default: TRUE; this has no effect at the moment)
##' @param sprate compute speech rate? (not in use yet)
##' @return a session's vocalisation feature stats
##' @rdname makeVocalStatsDataset
##' @importFrom stats median sd
##' @examples
##' \dontrun{
##' makeVocalStatsDataset(dir=c('ADReSS-IS2020-data/train/transcription/cc/',
##'                             'ADReSS-IS2020-data/train/transcription/cd/'),
##'                      sildir='ADReSS/silence/',
##'                      srdir='ADReSS/speech_rate/',
##'                      silsuffix='.wav-sil.csv')
##' }
##' @export
makeVocalStatsDataset <- function(dir=c("data/Pitt/Dementia/cookie",
                                        "data/Pitt/Control/cookie"),
                                  sildir=NULL,
                                  silsuffix='c.mp3.csv',
                                  srdir='data/Pitt/speech_rate/',
                                  srsuffix='sra',
                                  sprate=T)
{
    files <- c()
    for (d in dir) {
        files <- c(files, dir(d, pattern='.+\\.cha$', full.names=T))
    }
    l <- length(files)
    d <- data.frame(id=character(l), age=numeric(l),
                    gender=factor(rep('female', l),levels=c('male','female')),
                    mmse=numeric(l),
                    dur.mean=numeric(l), dur.sd=numeric(l), dur.median=numeric(l),
                    srate.mean=numeric(l), srate.max=numeric(l), srate.min=numeric(l),
                    srate.sd=numeric(l), number.utt=numeric(l),
                    sildur.mean=numeric(l), sildur.sd=numeric(l), sildur.median=numeric(l),
                    dur.max=numeric(l), sildur.max=numeric(l),
                    dur.min=numeric(l), sildur.min=numeric(l),
                    stringsAsFactors=F
                    )
    dx <- character(l)
    for (i in 1:l){
        f <- files[i]
        t <- makeSessionDataSet(f, sildir=sildir, silsuffix=silsuffix, srdir=srdir, srsuffix=srsuffix, sprate=sprate)
        d$id[i] <- gsub("(.cha|.*/)", "", f)
        d$age[i] <- as.numeric(t$ids$age[t$ids$PAR=='PAR'])
        d$gender[i] <- as.character(t$ids$gender[t$ids$PAR=='PAR'])
        d$mmse[i] <- as.numeric(t$ids$mmse[t$ids$PAR=='PAR'])
        ## get diagnostic (categorisation target)
        dx[i] <- as.character(t$ids$Dx[t$ids$PAR=='PAR'])
        ## limit analysis to patient speech
        t$trans <- t$trans[t$trans$speaker=='PAR',]
        dur <- t$trans$end - t$trans$begin
        d$dur.mean[i] <- mean(dur, na.rm=T)
        d$dur.sd[i] <- sd(dur, na.rm=T)
        d$dur.median[i] <- median(dur, na.rm=T)
        d$dur.max[i] <- max(dur, na.rm=T)
        d$dur.min[i] <- min(dur, na.rm=T)
        d$srate.mean[i] <- 1/mean(1/t$trans$speechrate) ## use hamonic mean for ratios 
        d$srate.max[i] <- max(t$trans$speechrate)
        d$srate.min[i] <- min(t$trans$speechrate)
        d$srate.sd[i] <- sd(t$trans$speechrate)
        d$number.utt[i] <- length(t$trans$utterance)
        sildur <- t$sil$end - t$sil$begin
        if (length(sildur) == 0 || sildur == 0) {sildur <- c(0,0)}
        #print(sildur)
        d$sildur.mean[i] <-mean(sildur, na.rm=T)
        d$sildur.sd[i] <- sd(sildur, na.rm=T)
        if (is.na(d$sildur.sd[i])) {d$sildur.sd[i] <- 0}
        d$sildur.median[i] <-median(sildur, na.rm=T)
        d$sildur.max[i] <- max(sildur, na.rm=T)
        d$sildur.min[i] <- min(sildur, na.rm=T)        
    }
    cbind(d, dx=dx)
}

##' makeSessionDataSet: create a data frame for a session (e.g. cookie scene description)
##'
##' @title makeSessionDataSet: create a data frame for a session (e.g. cookie scene description) based on .cha transcription files
##' @param f CHA file to read
##' @param sildir directory where silence profiles are stored
##' @param silsuffix suffix for silence files
##' @param srdir  directory where speech rate csv (1 value per utterance) files are stored 
##' @param srsuffix the suffix of the speech rate files (default: sre)
##' @param sprate estimate speech rate? (default: TRUE)
##' @return a speech session data frame
##' @author luzs
##' @export
makeSessionDataSet <- function(f, sildir=NULL, silsuffix='c.mp3.csv',
                               srdir='../data/ADReSS/speech_rate/',
                               srsuffix='sra',
                               sprate=T)
{
    cat('reading ', f, '\n') 
    t <- read.cha(f, sildir=sildir, silsuffix=silsuffix)
    if (sprate){
        srfile <- gsub('cha$', srsuffix, f)
        srfile <- gsub('.+/', srdir, srfile)
        t <- appendSpeechRate(t,file=srfile)
        t
    }
}

##' appendSpeechRate: append pre-generated speech rate data (see audioproc.R)
##'
##' @title appendSpeechRate: append pre-generated speech rate data to given dataframe t
##' @param t a table read through read.cha
##' @param file speech rate file 
##' @return dataframe t bound to speech rates per utterance
##' @author luzs
##' @importFrom utils read.csv
appendSpeechRate <- function(t,
                             file=NULL)
{
    speechrate <-  read.csv(file, col.names='speechrate')
    t$trans <- cbind(t$trans, speechrate)
    t
}


##' getSyllablesAndSilences: process Praat's grid for syllable nuclei, based on De Jong's approach 
##'
##' @title getSyllablesAndSilences: process Praat's grid for syllable nuclei
##' @param txtgrid Path to Praat grid file generated by praat-syllable-syllable-nuclei-v2
##' @return list of syllables and silences 
##' @author luzs
##' @references
##'  De Jong, N. H. and Wempe, T. (2009). Praat script to detect syllable nuclei
##'  and measure speech rate automatically. Behavior Research Methods,
##'  41(2):385–390, May.
##' @export
getSyllablesAndSilences <- function(txtgrid){
    grid <- readLines(txtgrid)
    sttier <- ''
    sttime <- 0
    sylstart <- c()
    silstart <- c()
    silend <- c()
    siltype <- c()
    silstarted <- FALSE
    for (l in grid) {
        if (length(grep('name = "syllables"', l, value=F))>0) {
            sttier <- 'syllables'
            next
        }
        if (length(grep('name = "silences"', l, value=F))>0){
            sttier <- 'silences'
            next
        }
        if (sttier == 'syllables') {
            r <- regexec('.*number *= *([0-9\\.]+)',l)
            if (r[[1]][1] > 0)
                sylstart <- c(sylstart, as.numeric(regmatches(l,r)[[1]][2]))
            next
        }
        if(sttier == 'silences'){
            r <- regmatches(l,regexec('.*(xmin|xmax) *= *([0-9\\.]+)',l))
            if (is.na(r[[1]][2]) || length(r[[1]][2])==0 ){
                r <- regmatches(l,regexec('.*(text) *= "*(.+)"',l))
                if (!is.na(r[[1]][2]) && length(r[[1]][2])>0 ){
                    siltype <- c(siltype, r[[1]][3])
                    cat('r-- ',r[[1]][3],'\n')
                }
                next
            }
            cat(l, '--', as.character(r),'\n')
            if (r[[1]][2] == 'xmin')
                silstart <- c(silstart, as.numeric(r[[1]][3]))
            else if (r[[1]][2] == 'xmax'){
                silend <- c(silend, as.numeric(r[[1]][3]))
                silstarted <- FALSE
            }
        }
    }
    list(sylstart=sylstart, silstart=silstart[-1], silend=silend[-1],
         siltype=siltype)
}
##' read.cha: read CHA transcription file (format used by DementiaBank)
##' 
##' @title read.cha read CHA transcription file (format used by DementiaBank)
##' @param file .cha file to reas
##' @param sildir silences directory
##' @param silsuffix silence files suffix
##' @return a list containing the PID, a dataframe containing the speaker IDs and demogrephics, and a dataframe containing the speaker IDs, transcribed utterances, start and en times, speech rates etc.  
##' @author luzs
##' @export
read.cha <- function(file, sildir=NULL, silsuffix='c.mp3.csv'){   
    text <- readLines(file)
    ## get rid of spurious formatting
    for(X in rev(grep('^\t.*', text)))
    {
        text[X-1] <- paste(c(text[X-1], sub('^\t', ' ', text[X])), collapse=' ')
    }
    text <- text[grep('^\t.*', text, invert=T)]
    cha <- list(pid=getPID(text),
                     ids=getIDs(text),
                trans=getTranscript(text),
                sil=getSilences(file, sildir=sildir, silsuffix=silsuffix)
                )
    cha
}

## Read description file
read.meta <- function(file="../data/Pitt/data.csv"){
    md <- read.csv(file, skip=2)
    md
}

##' getPIDs get study-wide unique patient IDs from CHA content
##' 
##' @title getIDs get  study-wide unique patient IDs from CHA content 
##' @param text a string vector containing the lines of a CHA file
##' @return a vector with participants IDs
##' @author luzs
##' @export
getPID <- function(text){
    for (line in text){
        pid <- regmatches(line,regexec('^@PID:\\W+(.+)\\W*',line))[[1]][2]
        if (!is.na(pid))
            break
    }
    if (is.na(pid))
        warning('No PID found!')
    pid
}

##' getIDs get speaker IDs from CHA content
##' 
##' @title getIDs get speaker role IDs (PAR, INV) and info from CHA content 
##' @param text a string vector containing the lines of a CHA file
##' @return a vector with participants IDs
##' @author luzs
##' @export
getIDs <- function(text){
    ids <- c() ##data.frame()
    for (line in text){
        id <- regmatches(line,regexec('^@ID:\\W+(.+)\\W*',line))[[1]][2]
        if (!is.na(id)){
            ids <- rbind(ids, strsplit(id, ';?\\|')[[1]])
        }
    }
    colnames(ids) <- c('language','UPMC','PAR','age','gender','Dx','Id','participant','mmse','ign')
    data.frame(ids,stringsAsFactors=F)
}

##' getTranscript
##' 
##' @title getTranscript: get transcription lines from .cha content 
##' @param text a string vector containing the lines of a CHA file
##' @return a list of transcriptions (participant and interviewer utterances)
##' @author luzs
##' @export
getTranscript <- function(text){
    trans <- c()
    buffer <- ''
    for (line in text){        
        tmp <- regmatches(line,
                          regexec('^\\*(PAR|INV):\\W+(.+)\\W*([0-9]+)_([0-9]+)',
                                  line))[[1]]
        if (length(tmp)>0){
            trans <- rbind(trans,
                           data.frame(speaker=tmp[2], utterance=I(tmp[3]),
                                      begin=as.numeric(tmp[4]),
                                      end=as.numeric(tmp[5])))
        }
    }
    trans
}

##' getSilences read silences file 
##'
##' @title getSilences read silences file 
##' @param file  CSV formatted silences file
##' @param sildir dir where silence files are
##' @param silsuffix ## suffix for silence files
##' @return silences dataframe
##' @author luzs
##' @export
getSilences <- function(file,
                        sildir=NULL,
                        silsuffix='c.mp3.csv' 
                        )
{
    sf <- sub('\\.cha', silsuffix,file)
    if (!is.null(sildir)){
        sf <- sub('.*/', sildir, sf)
    }
    else
        return(NULL)
    s <- read.csv(sf,sep='\t',header=F,
                  col.names=c('begin','end', 'desc'))
    s$begin <- s$begin * 1000
    s$end <- s$end * 1000
    s
}
