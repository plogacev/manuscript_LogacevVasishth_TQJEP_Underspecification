
read.swets.data <- function(fname) {
  d <- read.csv(fname, sep=";", as.is=T)
  id.cols <- c("subj","item","trial","X1","response.yes","resp.RT","X2","qtype","attachment","correct.response.code")
  colnames(d)[1:10] <- id.cols
  
  d$qtype <- nmap(d$qtype, c("1"="RC questions","2"="superficial","3"="occasional")) 
  d$attachment <- nmap(d$attachment, c("1"="ambiguous","2"="N1","3"="N2")) 
  d$response.yes <- nmap(d$response.yes, c("1"=1,"2"=0), asi) 
  
  d$questionNP = NA
  d$questionNP[d$qtype=="RC questions"] <- with(subset(d, qtype=="RC questions"), nmap(correct.response.code, c("1"="N1","2"="N2")))
  d$correct.response.yes = NA
  d$correct.response.yes[d$qtype!="RC questions"] <- nmap(d$correct.response.code[d$qtype!="RC questions"], c("1"=1,"2"=0), as.integer)
  d$correct.response.yes[d$qtype=="RC questions" & d$attachment!="ambiguous"] <- with(subset(d, qtype=="RC questions" & attachment!="ambiguous"), questionNP == attachment)
  d$response.n1 = NA
  d$response.n1 = d$questionNP=="N1" & d$response.yes==1 | d$questionNP=="N2" & d$response.yes==0
  
  critical.region <- c("on the balcony", "far too often", "quite a bit", "at the reception", "last summer", "constantly",
                       "on the bicycle", "most evenings", "the other day", "in thought", "at the party", "in public",
                       "after the accident", "last summer", "on the broken glass", "for the cause", "last year",
                       "in the ocean", "for lying", "by falling off a horse", "at night", "after the tragedy",
                       "in a lot of trouble", "to herself", "to an ice-cream cone", "a note", "to the party",
                       "an inordinate amount", "a lot", "very thoroughly", "with a low income", "in the mirror",
                       "very much", "a natural beauty", "in the office", "all the time")
  
  library(plyr)
  crit.lens <- sapply(strsplit(critical.region, split=" "), function(x) sapply(x, function(x) strlen(x)))
  crit.lens <- ldply(crit.lens, function(lens) {
    res <- rep(0, 5)
    res[1:length(lens)] <- lens
    names(res) <- paste("pc.len", 1:5, sep="")
    res
  })
  crit.lens <- cbind(item=1:36, crit.lens)
  d <- merge(d, crit.lens)
  
  library(reshape)
  regions <- c("the", "n1", "of", "the.1", "n2", "who", "sub.verb", "reflexive", "region9",  "matrix.verb", "w11","w12",
               "w13", "w14", "w15", "w16", "w17", "w18")
  x <- melt(d, measure.vars=regions)
  colnames(x) <- nmap(colnames(x), c('variable'='rid', 'value'='RT'))
  x$rid <- ordered(x$rid, levels=regions)
  x$pc.cnt <- with(x, (pc.len1!=0)+(pc.len2!=0)+(pc.len3!=0)+(pc.len4!=0)+(pc.len5!=0))# +pc.len3>0+pc.len4>0+pc.len5>0
  x$pc.avglen <- with(x, (pc.len1+pc.len2+pc.len3+pc.len4+pc.len5)/pc.cnt)
  
  x
}
