
asc <- function(x) as.character(x)
asi <- function(x) as.integer(x)
asd <- function(x) as.double(x)
asic  <- function(x) as.integer(as.character(x))
asdc  <- function(x) as.double(as.character(x))

round.mean <- function(x) round(mean(x))


se <- function(x)
     {
       y <- x[!is.na(x)] # remove the missing values, if any
       sqrt(var(as.vector(y))/length(y))
}

ci.width <- function (scores){
	stderr <- se(scores)
	len <- length(scores)
	return( qt(.975, df=len-1) * stderr )
}

ci <- function (scores){
	m <- mean(scores)
	w <- ci.width(scores)
	upper <- m + w
	lower <- m - w
	return(data.frame(lower=lower,upper=upper))
}



strlen <- function(str) {
	length(strsplit(as.character(str), split="")[[1]]) 
}

cbind.asym <- function(x, y) {
  x.len <- length(x)
  y.len <- length(y)
  if(x.len < y.len) {
    y <- c(y, rep(NA, y.len-x.len))
  } else {
    x <- c(y, rep(NA, x.len-y.len))
  }
  result <- cbind(x,y)
  return(result)
}

constructRoiInfo <- function(strings, wordinfo)
{
	cnts  <- c()
	means <- c()
	lens  <- c()
	for(string in strings) 
	{
		len <- strlen(string)
		words <- strsplit(string, split=" ")[[1]]
		wordCnt <- length(words)
		usedWords <- subset(wordinfo, word %in% words)
		if(length(usedWords[,1]) != wordCnt) {
			freqSum <- NA
		} else {
			freqSum <- sum(usedWords$freq.lemma)
		}
		meanFreq <- freqSum/wordCnt
		cnts <- append(cnts, wordCnt)
		means <- append(means, meanFreq)
		lens <- append(lens, len)
	}
	return(data.frame(wordCnt=cnts, meanFreq=means, lenChr=lens))
}


cmap <- function(d, from, to) {
	colnames(d) <- map(colnames(d), from, to)
	return(d)
}

fmap <- function(d, to, convert=NULL) {
   condVec <- as.character(d$condition)
   conds <- sort(unique(condVec))
   if(is.null(convert)) {
	convert = factor
   }	
   return(convert(map(condVec, conds, to)))
}

map <- function(vec, from, to, verbose=F) {
	newVec <- vec
	for( i in 1:length(from) ) {
                if(verbose){
                  print(from[i])
                }
		newVec[vec == from[i]] <- to[i]
	}
	return(newVec)
}

nmap <- function(vec, fromto, cast=I) cast(map(vec, names(fromto), fromto))

individual.rep <- function(vec, n) {
  ret <- c()
  if(length(n) == 1){
	n <- rep(n, length(vec))
  } 
  for( i in 1:length(vec)) {
    ret <- append(ret, rep(vec[i], n[i]))
  }
  return(ret)
}

in.range <- function(x,rng) {
  stopifnot(is.integer(x)|is.double(x))
  return(x>=rng[1] & x<=rng[2])
}



no.na <- function(x){x[is.na(x)] <- 0; return(x);}
no.zero <- function(x){x[x==0] <- NA; return(x);}
mean.se <- function(x) paste(round(mean(x)), ' (', round(se(x)), ')', sep="")
mean.2se <- function(x) paste(round(mean(x)), ' (', round(2*se(x)), ')', sep="")

mean.2se.n <- function(x) paste(round(mean(x)), ' (', round(2*se(x)),', n=', length(x), ')', sep="")


TOST <- function(mean1, mean2, theta, n1, n2, sigma) {
  d <- (mean2 - mean1) 
  t1 <- (d - theta)/(sigma * (sqrt((1/n1) + (1/n2))))
  t2 <- (d + theta)/(sigma * (sqrt((1/n1) + (1/n2))))
  tcrit <- qt(0.95, (n1 + n2 - 2))
  if ((t1 < -tcrit) && (t2 > tcrit)) {
    print(t1) 
    print(t2) 
    print(tcrit) 
    print(c("Equivalent"))
   } else{ 
    print(c("Failed to show equivalence"))
   }
}



mean.se.simple <- function(RT, subject, condition) {

 aggregate(list(RT=RT), list(condition=condition), function(x) c(M=mean(x), SE=se(x), N=length(x)))
  M <- tapply(RT, condition, mean)
  SE <- tapply(RT, condition, se)
  N <- tapply(RT, condition, length)
}

mean.se.cousineau <- function(RT, subject, condition, conditions.cnt=0) {
  library(plyr, warn.conflicts=TRUE)
  library(reshape, warn.conflicts=TRUE)
  if(conditions.cnt) correction <- conditions.cnt/(conditions.cnt-1)
  else correction <- 1
  d <- data.frame(RT, subject, condition)
  d$GM <- mean(tapply(RT, asc(subject), mean))
  d <- ddply(d, .(asc(subject)), transform, RT.w = RT - mean(RT) + GM)  
  temp<-melt(d, id.var=c("subject","condition"), measure.var="RT.w")
  (M.id.w <- cast(temp, condition  ~ ., 
          function(x) { cur.var <- var(x)*correction;
  			    c(M=mean(x), SE=sqrt(cur.var/length(x)), N=length(x) ) 
		})) 
}

mean.se.cousineau.proportion <- function (DV, subject, condition, conditions.cnt=0) {
    library(plyr)
    library(reshape)
    if(conditions.cnt) correction <- conditions.cnt/(conditions.cnt-1)
    else correction <- 1
    d <- data.frame(DV, subject, condition)
    d$GM <- mean(tapply(DV, asc(subject), mean))
    d <- ddply(d, .(asc(subject)), transform, DV.w = DV - mean(DV) + 
        GM)
    temp <- melt(d, id.var = c("subject", "condition"), measure.var = "DV.w")
    (M.id.w <- cast(temp, condition ~ ., function(x) { 
	cur.var <- mean(x)*(1-mean(x))*correction
	c(M = mean(x), SE = sqrt(cur.var/length(x)), N = length(x))
    }))
}

