
# load Swets' data
source("~/Documents/MyPapers/models_of_underspecification/Data_Swets_et_al/RScriptSwetsData.R")
d <- read.swets.data("~/Documents/MyPapers/models_of_underspecification/Data_Swets_et_al/Data_All.csv")

d = subset(d, qtype%in%c('RC questions'))
d = d[with(d, order(subj,trial,item,rid)),]
d.critical <- subset(d, rid=="reflexive") 
d.critical2 <- subset(d, rid=="region9") 
d.critical$RT = d.critical$RT + d.critical2$RT

# find participants with more than 35% correct answers in the unambiguous conditions
subject.acc = with(d.critical, tapply(response.yes==correct.response.yes, list(subj, attachment), mean ))
nmap(asc(d.critical$subj),  subject.acc[,'N1'], as.double) -> d.critical$acc.N1
nmap(asc(d.critical$subj),  subject.acc[,'N2'], as.double) -> d.critical$acc.N2

acc.min = .5
d.qrc = subset(d.critical, acc.N1 > acc.min & acc.N2 > acc.min)
d.excluded = unique( subset(d.critical, !(acc.N1 > acc.min & acc.N2 > acc.min))[,c('subj','qtype','acc.N1','acc.N2')] )

excluded.subjects.rc.N1 = unique(subset(d.excluded, acc.N1 > acc.min)$subj)
excluded.subjects.rc.N2 = unique(subset(d.excluded, acc.N2 > acc.min)$subj)

#d.qrc$attachment = ordered(d.qrc$attachment, levels=c("N1","N2","ambiguous"))
d.qrc$attachment.n2 = with(d.qrc, ifelse(attachment=="N2", 1, 0))
d.qrc$response.correct = NA
d.qrc$response.correct[d.qrc$attachment!="ambiguous"] = with(subset(d.qrc, attachment!="ambiguous"), ifelse(response.yes==correct.response.yes, 1, 0))
d.qrc$response.yes.str = nmap(d.qrc$response.yes, c('0'="`no'", '1'="`yes'"))
d.qrc$response.correct.str = nmap(asc(d.qrc$response.correct), c('1'='correct response', '0'='incorrect response'))
d.qrc$crit.region.cnt = d.qrc$pc.cnt + 1

# exclude outlier
d.qrc = subset(d.qrc, resp.RT < 15000 )

# change the variable coding for stan
d <- subset(d.qrc, qtype=="RC questions")[,c('item','subj','trial')]
d$iv_questN1 <- as.integer(d.qrc$questionNP=="N1")
d$crit_region_cnt <- d.qrc$crit.region.cnt
d$item <- as.integer(as.factor(d.qrc$item))
d$subj <- as.integer(as.factor(d.qrc$subj))
d$response_yes <- d.qrc$response.yes
d$response_RT <- d.qrc$resp.RT
d$reading_time <- d.qrc$RT
d$iv_cond <- ifelse(d.qrc$attachment=="ambiguous", 0, ifelse(d.qrc$attachment=="N1", 1, 2))

d.qrc.stan <- d

save(d.qrc, d.qrc.stan, file="SwetsDataForModelling.rda")

