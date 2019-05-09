#Clare Adams
#2019 05 09
#Playing around tihe Decipher in the name of making primers


library(DECIPHER)
?DesignPrimers

setwd("~/Documents/Bioinformatics Scripts & Notes/DECIPHERtest/DECIPHERtest")

# STEP 1
# READ IN FILES AND GET THEM INTO AN OK FORMAT

#Specify the fasta sequences of the target area
fas <- "PauaCOXII.fasta"

# specify a path for where to write the sequence database

#Creating the seq database in the memory
dbConn <- dbConnect(SQLite(), ":memory:")
N <- Seqs2DB(fas, "FASTA", dbConn, "")
#Outputs number of total sequences in the table of sequences

N
#should equal 108 now

#Every sequence here is its own thing
desc <- as.character(seq_len(N))

#Every sequence here is grouped by a distance matrix and then clustered based on the distance matrix / IDs
dna <- readDNAStringSet(fas)
d <- DistanceMatrix(dna, type="dist")
IdClusters(d, method="complete", cutoff=0.05, showPlot=TRUE, type="dendrogram")
#Chuck the clusters into an object or something that can be called
desc<-IdClusters(d, method="complete", cutoff=0.05, showPlot=TRUE, type="dendrogram")
#Can get at lower branches by specifying which branch you want to look at
unique(desc[[2]])

#Decided to go with the other ones because apparently everything is its own thing anyways (though probs not true?)

#ID the sequences with a number
desc <- as.character(seq_len(N))
#Show the unique descriptors (aka the number of each file)
unique(desc)

#Add the unique names for the sequences (in this case, their number) to the database as an ID for each seq
Add2DB(data.frame(identifier=desc), dbConn)

# STEP 2
# Set up some parameters for designing primers

#> # Designing primers for sequencing experiments:
#> TYPE <- "sequence"
#> MIN_SIZE <- 300 # base pairs
#> MAX_SIZE <- 700
#> RESOLUTION <- 5 # k-mer signature
#> LEVELS <- 5 # max number of each k-mer

TYPE <- "sequence"
MIN_SIZE <-175
MAX_SIZE <-350
RESOLUTION <- 5
LEVELS <- 5

#Now, design primers

primers <- DesignSignatures(dbConn,
                            type=TYPE,
                            minProductSize=MIN_SIZE,
                            maxProductSize=MAX_SIZE,
                            resolution=RESOLUTION,
                            levels=LEVELS)
#It does things now

# STEP 3
# Look at things

#This code gives out a forward and reverse primer
primers[which.max(primers$score),]

#PSET to 1 to get the best/first primer
PSET <- 1

#Some code which does things apparently
dna <- SearchDB(dbConn,
                remove="all",
                nameBy="identifier",
                clause="row_names =
                (select min(row_names) from Seqs as S
                where S.identifier = Seqs.identifier)",
                verbose=FALSE)
#code which pulls out the f primer and r primer
f_primer <- DNAStringSet(primers$forward_primer[PSET])
r_primer <- DNAStringSet(primers$reverse_primer[PSET])

#something about
patterns <- c(f_primer,
              reverseComplement(r_primer))
BrowseSeqs(dna,
           patterns=patterns)

#Other stuff
PSET <- which.max(primers$score)
f_primer <- DNAString(primers$forward_primer[PSET])
r_primer <- DNAString(primers$reverse_primer[PSET])
r_primer <- reverseComplement(r_primer)
ids <- dbGetQuery(dbConn, "select distinct identifier from Seqs")
ids <- ids$identifier
if (TYPE=="sequence") {
  signatures <- matrix(0, nrow=4^RESOLUTION, ncol=length(ids))
} else if (TYPE=="melt") {
  signatures <- matrix(0, nrow=length(RESOLUTION), ncol=length(ids))
} else { # TYPE=="length"
  signatures <- matrix(0, nrow=length(RESOLUTION) - 1, ncol=length(ids))
}
colnames(signatures) <- abbreviate(ids, 15)

for (i in seq_along(ids)) {
  dna <- SearchDB(dbConn, identifier=ids[i], remove="all", verbose=FALSE)
  amplicons <- matchLRPatterns(f_primer, r_primer,
                               MAX_SIZE, unlist(dna),
                               max.Lmismatch=2, max.Rmismatch=2,
                               Lfixed="subject", Rfixed="subject")
  amplicons <- as(amplicons, "DNAStringSet")
  if (length(amplicons)==0)
    next
  if (TYPE=="sequence") {
    signature <- oligonucleotideFrequency(amplicons, RESOLUTION)
    signatures[, i] <- colMeans(signature)
  } else if (TYPE=="melt") {
    signature <- MeltDNA(amplicons, "melt curves", RESOLUTION)
    # weight melting curves by their amlicon's width
    signature <- t(signature)*width(amplicons)
    signatures[, i] <- colSums(signature)/sum(width(amplicons))
  } else { # TYPE=="length"
    signature <- .bincode(width(amplicons), RESOLUTION)
    for (j in signature[which(!is.na(signature))])
      signatures[j, i] <- signatures[j, i] + 1/length(signature)
  }
}

if (TYPE=="sequence") {
  d <- dist(t(signatures), "minkowski", p=1) # L1-Norm
  IdClusters(as.matrix(d), showPlot=T, verbose=FALSE)
  mtext(paste(RESOLUTION, "-mer Profile Distance", sep=""),
        side=2, padj=-4)
} else if (TYPE=="melt") {
  matplot(RESOLUTION, signatures, type="l",
          xlab="Temperature (degrees Celsius)", ylab="Average Helicity")
} else { # TYPE=="length"
  if (length(ids) > 20) {
    plot(NA,xlim=c(0.5, length(ids) + 0.5), ylim=range(RESOLUTION),
         xlab="Group Index", ylab="Amplicon Length",
         yaxs="i", xaxs="i")
    axis(1, at=1:length(ids), labels=FALSE, tck=-0.01)
  } else {
    plot(NA,
         xlim=c(0.5, length(ids) + 0.5), ylim=range(RESOLUTION),
         xlab="", ylab="Amplicon Length",
         yaxs="i", xaxs="i", xaxt="n")
    axis(1, at=1:length(ids), labels=abbreviate(ids, 7), las=2)
  }
  xaxs <- RESOLUTION[-1] - diff(RESOLUTION)/2 # average lengths
  for (i in seq_along(ids)) {
    w <- which(signatures[, i] > 0)
    if (length(w) > 0)
      segments(i - 0.45, xaxs[w], i + 0.45, xaxs[w], lwd=2)
  }
}