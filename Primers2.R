#Clare Adams
#2019 05 09
#DECIPHER Package play Round 2
#Using the Design Group-Specific Primers

#It seems like they want to capitalize on the "mismatches near the 3' end" to increase specificity by decreasing polymerase binding/extension

#Install the OligoArrayAux to predict primer binding

library(DECIPHER)
system("/Users/clareadams/Documents/Bioinformatics/DECIPHERtest/DECIPHERtest/bin/hybrid-min -V")

#Set working directory
setwd("~/Documents/Bioinformatics/DECIPHERtest/DECIPHERtest")

#Read in your stuff

fas <- "PauaCOXII.fasta"
fas2 <- "mito.all.fasta"
fas2<-"PauaCRaligned.fasta"
fas2<-"CR-unique-longer.fasta"

#Create a sequence database in memory
dbConn <- dbConnect(SQLite(), ":memory:")
Seqs2DB(fas2, "FASTA", dbConn, "mtDNA")

#Defining groups of related sequences

# get the FASTA record description
desc <- dbGetQuery(dbConn, "select description from Seqs")

# parse the sequence description to obtain the species name
desc <- unlist(lapply(strsplit(desc$description, "mtDNA", fixed=TRUE),
                      function(x) return(x[length(x)])))
#fix some formatting stuff
desc <- gsub("sp. ", "", desc, perl=TRUE)
desc <- gsub("sp_", "", desc, perl=TRUE)

#do some other stuff
desc <- unlist(lapply(strsplit(desc, " ", fixed=TRUE), function(x) return(x[1])))

#see each unique name; probably the first thing
unique(desc)

#[1] "H-rufescens"            "H-discus-hannai"        "H-discus-hannai-H-iris"
#[4] "H-rubra"                "H-iris"


#Need to add the names to the database as identifiers of each sequence
Add2DB(data.frame(identifier=desc), dbConn)

# Step 2 - Designing Primers

#Make tiles/k-mers of ~26-27 nucleotides each with ~10 permutations
tiles <- TileSeqs(dbConn, add2tbl="Tiles", minCoverage=1, processors = 8)

#Alternatively, you can run this; minimum coverage is a little less strict with 12 processors going (for boros)
#tiles <- TileSeqs(dbConn, add2tbl="Tiles", minCoverage=.99, processors = 12)

#Examine the tiles
head(tiles)

#Design all possible primers
primers <- DesignPrimers(tiles, identifier="H-iris",
minCoverage=1, minGroupCoverage=1, processors = 12)

#See primers
primers[1,]

#Designing the best possible primer set
primers <- DesignPrimers(tiles, identifier="H-iris", minCoverage=1,
minGroupCoverage=1, numPrimerSets=5, maxSearchSize=20, processors = 12)

#Look at your primers
head(primers)

