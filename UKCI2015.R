# UKCI2015.R
# The 15th UK Workshop on Computational Intelligence UKCI2015, Exeter, 7th-9th Sept 2015
# 1. Draft Paper Submission Deadline: 8th June 2015, 
# 2. Notification of paper acceptance: 6th July 2015: Accepted!
# 3. Camera ready deadline: 4th August 2015. Done!
#               Up to eight pages.
#               source("http://bioconductor.org/biocLite.R")
#               biocLite("GOstats")

library(GOstats)
library(GO.db)
library(org.Sc.sgd.db)
library(sand)
library(housekeep)
library(caret)
library(kernlab)
library(dplyr)
library(tidyr)
library(VennDiagram)
library(ggplot2)
library(xtable)
library(bioassayR)
library(clusterProfiler)
library(mygene)
library(DOSE)
library(igraph)

# load in the STITCH data, proteins linked to Alzheimers and their numerous interactions
plist <- c('APP','PSEN1','PSEN2','APOE','BACE1') 
BACE1 <- file.path('C://R-files//UKCI2015', 'BACE1.txt') %>% read.delim(na.strings='')
PSEN1 <- file.path('C://R-files//UKCI2015', 'PSEN1.txt') %>% read.delim(na.strings='')
PSEN2 <- file.path('C://R-files//UKCI2015', 'PSEN2.txt') %>% read.delim(na.strings='')
APOE  <- file.path('C://R-files//UKCI2015', 'APOE.txt') %>% read.delim(na.strings='')
APP   <- file.path('C://R-files//UKCI2015', 'APP.txt') %>% read.delim(na.strings='') 

# Keep P1 and P2 (interacting proteins) with confidence score, discard the other 13 variables.
PSEN1<-PSEN1[,c("X.node1","node2","combined_score")]
PSEN2<-PSEN2[,c("X.node1","node2","combined_score")]
APOE<-APOE[,c("X.node1","node2","combined_score")]
APP<-APP[,c("X.node1","node2","combined_score")]
BACE1<-BACE1[,c("X.node1","node2","combined_score")]

# combine interactions and remove duplicates
total <- rbind(PSEN1,PSEN2,APOE,APP,BACE1) 
total <- unique(total)

# SIDER2 data (side-effects)
#db_maps   <- file.path('C://R-files//drugbank', 'label_mapping.tsv.gz') %>% read.delim(na.strings='',header = FALSE) 
db_med_effects  <- file.path('C://R-files//drugbank', 'meddra_adverse_effects.tsv.gz') %>% read.delim(na.strings='',header = FALSE) 
#db_adverse_effects <- file.path('C://R-files//drugbank', 'adverse_effects_raw.tsv.gz') %>% read.delim(na.strings='',header = FALSE) 
#db_indications <- file.path('C://R-files//drugbank', 'indications_raw.tsv.gz') %>% read.delim(na.strings='',header = FALSE) 
#db_fp      <- file.path('C://R-files//drugbank', 'meddra_freq_parsed.tsv.gz') %>% read.delim(na.strings='',header = FALSE) 

# load in data from the PharmGKB website
#GKB_offsides <- file.path('C://R-files//himmelstein', 'offsides.tsv') %>% read.delim(na.strings='')
#GKB_drug2drug <- file.path('C://R-files//himmelstein', 'twosides.tsv') %>% read.delim(na.strings='')
#GKB_boxwarn <- file.path('C://R-files//himmelstein', 'box_warnings.tsv') %>% read.delim(na.strings='')

# clean up unneeded variables
#rm(APOE,APP,BACE1,PSEN1,PSEN2)

#newdata <- sider.df[order(sider.df$pubchem_cid),]  # SORT DATA BY 

# read in Himmelstein SIDER2 data
#sider.df <- file.path('C://R-files//drugbank', 'sider2-processed.txt') %>% read.delim(na.strings='')

# This is correct for number of drugs per side-effect
results<-table(db_med_effects$V4)
results<-sort(results,decreasing=TRUE)
plot(results,xlab="Number of side-effects",ylab="Number of drugs",cex.main=1)
abline(v=(seq(0,1000,50)), col="darkgray", lty="dotted")
abline(h=(seq(0,1500,200)), col="darkgray", lty="dotted")


# This is correct for number of side -effects per drug.
results<-table(db_med_effects$V5)
results<-sort(results,decreasing=TRUE)
plot(results,xlab="Number of drugs",ylab="Number of side-effects",cex.main=1)
abline(v=(seq(0,4500,200)), col="darkgray", lty="dotted")
abline(h=(seq(0,2500,200)), col="darkgray", lty="dotted")

# The three common drugs used to treat alzheimers: 'donepezil','galantamine','rivastigmine'
# get the side effects for each one. NOTE: spellings important!!
donepezil<-db_med_effects$V5[grep('donepezil',db_med_effects$V4)]
galantamine<-db_med_effects$V5[grep('galantamine',db_med_effects$V4)]
rivastigmine<-db_med_effects$V5[grep('rivastigmine',db_med_effects$V4)]

# Now get only unique side-effects for each drug as duplications occur.
donepezil<-unique(donepezil)
galantamine<-unique(galantamine)
rivastigmine<-unique(rivastigmine)

# Now create the Venn diagram, saving it to a PDF file
plot.new()
venn.plot <- venn.diagram(list(donepezil,galantamine,rivastigmine), 
                          NULL, 
                          fill=c("red", "blue","green"), 
                          alpha=c(0.5,0.5,0.5), 
                          cex = 2, 
                          cat.fontface=2, 
                          margins =c(10,10),
                          cat.cex=2,
                          #main = "Venn Diagram showing shared side effects for donepezil,galantamine,rivastigmine",
                          category.names=c("donepezil", "galantamine","rivastigmine"))
grid.draw(venn.plot)

#---------------- Now determine combined side effects of the three drugs ------------
# need to convert from factors to strings
donepezil <- as.character(donepezil)
galantamine <- as.character(galantamine)
rivastigmine <- as.character(rivastigmine)

universe <- unique(c(donepezil,galantamine,rivastigmine))
GroupA <-universe %in% donepezil
GroupB <-universe %in% galantamine
GroupC <-universe %in% rivastigmine

## Side effects that are in GroupA and in GroupB but not in GroupD, 97 common side effects as per Venn.
allSE <- universe[GroupA & GroupB & GroupC]

# clean up memory
rm(GroupA,GroupB,GroupC,universe)

#--------- Now determine which drugs have similar side effects to our universe of side-effects -----
# But how many SE's would be deemed 'enough' for a drug to be classed as similar????

candDrug <- db_med_effects$V5[grep('donepezil',db_med_effects$V4)]
candSE <- db_med_effects$V5[grep('donepezil',db_med_effects$V5)]

# 215,852 observations, make a dataframe with drug name and number of side effects in common from the universe
# Give a count of SE commonality for each drug in list

NoDrugs <- unique(db_med_effects$V4) # We have 996 drugs

Names <- letters[1:5];
Dates<- 1:length(NoDrugs);
allDrugs<- data.frame(drugnames=Dates, NoSideEffects = vector(mode="numeric", length=length(Dates)), coverage=vector(mode="numeric",length=length(Dates))); 
allDrugNames <- lapply(NoDrugs, as.character)

# TAKES FIVE MINUTES TO CALCULATE
for (i in 1:length(NoDrugs)){
  #i=50
  allDrugs[i,1]<-as.character(allDrugNames[i])
  TempDrug <- db_med_effects$V5[grep(allDrugNames[i],db_med_effects$V4)]
  #TempSE <- db_med_effects$V5[grep(allDrugNames[i],db_med_effects$V5)]
  allDrugs[i,2] <- length(unique(TempDrug)); # side effects found for each drug
  TempDrug <- unique(TempDrug)
  coverage <- allSE[TempDrug]
  coverage <- coverage[!is.na(coverage)] # get rid of NA when we dont get a matching side-effect
  coverage <- (length(coverage)/length(allSE))*100; # percentage coverage of joint SE of three alzheiemers drugs

  allDrugs[i,3] <- coverage
    
  }
  
summary(allDrugs[,3])

#----------- get list of drugs with at least 10% shared side effects with our three Alzheimers drugs

candidates <- subset(allDrugs, allDrugs$coverage >= 10.0)
candidates <- candidates[order(-candidates$coverage),]  # sort descending so we get top ten for latex table.

# create a LaTex table for my document using the XTABLE package.

tli.table <- xtable(candidates)
digits(tli.table)[c(2,6)] <- 0
print(tli.table,floating=FALSE)

# ----- Identify the targets of these drugs, build protein networks
#   protein networks are calculate here ------------------
curr_drugs <- read.delim('C:\\R-files\\drugbank\\similardrugs1.csv', header=TRUE,sep=',') 

disease <- as.matrix(curr_drugs[,c(2,4)])    # V2=drug V4 =targets
g <- graph.edgelist(disease,directed=FALSE)

nodesize=degree(g)*2

ad <- get.adjacency(g)
nodecolor=character(ncol(ad))  # create a character for every column in adjaceny matrix,
x <- 1:ncol(ad)

nodelabel<-V(g)$name
d1<-as.character(disease[,1])
d2<-as.character(disease[,2])

for ( i in 1:ncol(ad)){
  z<-(sapply(nodelabel[i],grep,d2))
  if(is.integer(z)){
    nodecolor[i]<-"pink"}
  y<-(sapply(nodelabel[i],grep,d1)) ## Disease
  if(is.integer(y)){   
    nodecolor[i]<-"lightblue"}                              
}

#--- TKPLOT allows you drag nodes around and create a better graph --------
tkplot(g,layout = layout.fruchterman.reingold,vertex.label = nodelabel,
       vertex.label.color= "black",vertex.size=nodesize, vertex.color=nodecolor,
       edge.arrow.size=0, edge.curved=FALSE)
#-------------------------------------------------------------------------


# BUILD NETWORKS FOR PSEN1; PSEN2; APOE; APP AND BACE1
# AND WORK OUT STATISTICS

# PSEN1 PSEN1<-PSEN1[,c("X.node1","node2")]
# PSEN2
# APOE
# APP
# BACE1
PSEN1<-PSEN1[,c("X.node1","node2","combined_score")]
PSEN2<-PSEN2[,c("X.node1","node2")]
APOE<-APOE[,c("X.node1","node2","combined_score")]
APP<-APP[,c("X.node1","node2","combined_score")]
BACE1<-BACE1[,c("X.node1","node2","combined_score")]


psen1 <-graph.data.frame(PSEN1)
psen2 <-graph.data.frame(PSEN2)
apoe <-graph.data.frame(APOE)
app <- graph.data.frame(APP)
bace1 <-graph.data.frame(BACE1)


ecount(bace1)
vcount(bace1)
clusters(bace1)
diameter(bace1,weights=NA)
transitivity(bace1)
reciprocity(bace1)
is.connected(bace1)
average.path.length(bace1)

###==============================================
done <- read.delim('C:\\R-files\\UKCI2015\\donepezil.txt', header=TRUE,sep='\t') 
done <- done[,1:2]

done <-graph.data.frame(done)
ecount(done)
vcount(done)
diameter(done,weights=NA)
transitivity(done)
reciprocity(done)
is.connected(done)
average.path.length(done)
######################################################################################

gala <- read.delim('C:\\R-files\\UKCI2015\\galantamine.txt', header=TRUE,sep='\t') 
gala <- gala[,1:2]

gala <-graph.data.frame(gala)
ecount(gala)
vcount(gala)
diameter(gala,weights=NA)
transitivity(gala)
reciprocity(gala)
is.connected(gala)
average.path.length(gala)

######################################################################################
riva <- read.delim('C:\\R-files\\UKCI2015\\rivastigmine.txt', header=TRUE,sep='\t') 
riva <- riva[,1:2]

riva <-graph.data.frame(riva)
ecount(riva)
vcount(riva)
diameter(riva,weights=NA)
transitivity(riva)
reciprocity(riva)
is.connected(riva)

######################################################################################


#--------------------- ONTOLOGY INTEGRATION -------------------------------
# Gene ontology (GO) and disease ontology (DO) can be used to enrich the statistical information 
# gained from the graphs by integrating them with biological knowledge.

# alzheimers DOID:10652
# parkinsons DOID:14330
# schizophrenia DOID:5419
# psychosis DOID:8646
# depression DOID:1596
# pain disorder DOID:0060164


# ropinrole -> parkinsons 
# tramadol -> analgesic
# pregabalin -> anticonvulsant
# citalopram -> antidepressants
# paroxetine -> antidepressants
# aripiprazole -> schizophrenia
# olanzapine -> antipsychotics

# Get the semantic similarity between our diseases using doSim() function.
hgncList<-c("APP","PSEN1","PSEN2","APOE");
listDO<-c("DOID:10652","DOID:14330","DOID:5419","DOID:8646","DOID:1596","DOID:0060164","DOID:1826","DOID:750","DOID:14262","DOID:13938")
listDO1<-c("DOID:14330","DOID:5419","DOID:8646","DOID:1596","DOID:0060164")
alzDO="DOID:10652"
s<-doSim(listDO,listDO,measure="Wang")
simplot(s,color.low="white",color.high="red",labs="true",digits=2,labs.size=5,
        font.size=14,xlab="",ylab="")

# ----- change row and column names from DOID numbers to disease names -----
DoNames<-c("alzheimers","parkinsons","schizophrenia","psychosis","depression","pain_disorder","epilepsy","ulcer","candidiasis","amenorrhea")
colnames(s)<-DoNames
rownames(s)<-DoNames

# -------------- gene ontology --------------------
# first some Gene ID conversion from the clusterprofiler function bitr()

# ------------ doneprezil interacting proteins enrichment
dGene<-done[,1]
gene <- bitr(dGene,fromType = "SYMBOL",toType="ENTREZID",annoDb="org.Hs.eg.db")
ggo <- groupGO(gene=gene$ENTREZID, organism="human", ont="BP",level=3, readable=TRUE)
#head(summary(ggo))

# ggo@result[1:5,1:4] # is a S4 object so use '@' instead of '$' 
ggo <- ggo@result[ggo@result$Count > 0,]  # only want results with our proteins!
ggo<-as.data.frame(ggo)           # convert to df
goDone <- ggo[order(-ggo$Count),]   # sort results in descending order i.e. most important first

xtable(goDone[1:10,1:4])


# ------------ Galantamine interacting proteins enrichment
dGene<-gala[,1]
gene <- bitr(dGene,fromType = "SYMBOL",toType="ENTREZID",annoDb="org.Hs.eg.db")
ggo <- groupGO(gene=gene$ENTREZID, organism="human", ont="BP",level=3, readable=TRUE)
#head(summary(ggo))

# ggo@result[1:5,1:4] # is a S4 object so use '@' instead of '$' 
ggo <- ggo@result[ggo@result$Count > 0,]  # only want results with our proteins!
ggo<-as.data.frame(ggo)           # convert to df
goGala <- ggo[order(-ggo$Count),]   # sort results in descending order i.e. most important first

xtable(goGala[1:10,1:4])


# ------------ Rivastigmine interacting proteins enrichment
dGene<-riva[,1]
gene <- bitr(dGene,fromType = "SYMBOL",toType="ENTREZID",annoDb="org.Hs.eg.db")
ggo <- groupGO(gene=gene$ENTREZID, organism="human", ont="BP",level=3, readable=TRUE)
#head(summary(ggo))

# ggo@result[1:5,1:4] # is a S4 object so use '@' instead of '$' 
ggo <- ggo@result[ggo@result$Count > 0,]  # only want results with our proteins!
ggo<-as.data.frame(ggo)           # convert to df
goRiva <- ggo[order(-ggo$Count),]   # sort results in descending order i.e. most important first

xtable(goRiva[1:10,1:4])


# plot a nice combined histogram of the GO enrichments for the three drugs
g1 <- data.frame(len=goGala$Count) # Ensure length is used otherwise ggplot2 wont bloody work.
r1 <- data.frame(len=goRiva$Count)
d1 <- data.frame(len=goDone$Count)

names(g1)[names(g1)=="goGala$Count"] <- "count"  # Replace the stupid names "goXXX$Count" with "count"
names(r1)[names(r1)=="goRiva$Count"] <- "count"
names(d1)[names(d1)=="goDone$Count"] <- "count"

# Give each drug its full name.
g1$drug <- 'Galantamine'
r1$drug <- 'Rivastigmine'
d1$drug <- 'Doneprezil'
  
#and combine into your new data frame drugLen
drugLen <- rbind(d1,r1,g1)
ggplot(drugLen, aes(len, fill = drug)) + geom_bar(pos="dodge" ) + xlab("GO terms with matching proteins") +  ylab("Frequency")


#----------- DOSE package ----------
alzDO="DOID:13938"
class(alzDO) = "DO"
TERM2NAME(alzDO)

#------- mygene package -------------
dGene<-c("APP","PSEN1","PSEN2","APOE");
t<-getGene(dGene,fields=c("symbol"))
t$symbol




