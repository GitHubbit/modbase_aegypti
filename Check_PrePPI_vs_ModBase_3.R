# this script will check to see if ModBase can generate more structural models than what PrePPI has given
# us = PrePPI

rm(list=ls())

install.packages("data.table")
install.packages("tidyr")
install.packages('VennDiagram')
install.packages("stringr")
library(data.table)
library(tidyr)
library(stringr)
library(VennDiagram)

# put checks in to see if were able to map uniprots to any hfpd or other id

setwd('/Users/Kamileh/Code/shapira_lab/PrePPI/aegypti')

modbase.summary <- fread("./aedes_aegypti_2016_214300.summary.txt", header=TRUE)
nrow(modbase.summary) # how many Modbase models do we start out with
modbase.summary <- subset(modbase.summary, modbase.summary$evalue < 0.0001 & modbase.summary$ga341 > 0.7) # filter modbase file for reliable models only (GA341 > 0.7 & evalue < 0.0001)
nrow(modbase.summary) # how many RELIABLE Modbase models do we end up with 


#### output pdb files to textfile for structural neighbor search part....


#################
gene.id <- substr(modbase.summary$database_id, start = 1, stop = 10) # extract protein id from database_id 
uniq.gene.ids <- data.frame(unique(gene.id)) # get the unique gene ids
# length(uniq.gene.ids)
write.table(uniq.gene.ids, "uniq.gene.ids.txt", sep="\n",  # output it to text file
            row.names = FALSE, 
            col.names = FALSE,
            quote = FALSE)

gene.ids.to.uniprots <- fread("./aegypti.modbase.to.uniprot.mapped", header=TRUE) # go to http://www.uniprot.org/uploadlists/ and paste output with paramters GeneName to UniprotKB mapping, get back file and open read into table
colnames(gene.ids.to.uniprots)[2] <- "gene.ids"

n_occur <- data.frame(table(gene.ids.to.uniprots$gene.ids)) # find out which gene ids were duplicated
duplicate.gene.ids <- gene.ids.to.uniprots[gene.ids.to.uniprots$gene.ids %in% n_occur$Var1[n_occur$Freq > 1],] # find out which db ids were duplicated
num.of.duplicate.gene.ids <- data.frame(unique(duplicate.gene.ids$gene.ids)) # gene ids which returned more than 1 uniprot
# find out if there are any duplicate uniprots
n_occur <- data.frame(table(gene.ids.to.uniprots$Entry)) # find out which uniprots were duplicated
duplicate.uniprots <- gene.ids.to.uniprots[gene.ids.to.uniprots$Entry %in% n_occur$Var1[n_occur$Freq > 1],]

# make output file of uniprots and corresponding pdb model file name
# -- paste the gene ids to the modbase summary dataframe, so that we can line up the column of pdb file names and uniprots against it
modbase.summary$pdb.file.paths <- modbase.summary$database_id
modbase.summary$pdb.file.paths <- paste("/ifs/c2b2/whatever/", modbase.summary$pdb.file.paths, sep = "")
modbase.summary$pdb.file.paths <- paste(modbase.summary$pdb.file.paths, ".pdb.xz", sep = "")

## fix uniprots returned so that each uniprot is in separate column (otherwise have comma-separated uniprots returned)
gene.ids.missing.in.modbase <- data.frame(uniq.gene.ids$unique.gene.id.[!(uniq.gene.ids$unique.gene.id. %in% gene.ids.to.uniprots$gene.ids)]) # gene ids submitted to UniProt that returned no uniprots (some are missing bc 2 gene ids will return the same uniprot (true for 8 uniprots))
colnames(gene.ids.missing.in.modbase) <- "gene.ids"
modbase.gene.ids <- data.frame(str_split_fixed(gene.ids.to.uniprots$gene.ids, ",", 2))
modbase.gene.ids1 <- data.frame(gene.ids.to.uniprots$Entry, modbase.gene.ids$X1)
modbase.gene.ids2 <- data.frame(gene.ids.to.uniprots$Entry, modbase.gene.ids$X2)
modbase.gene.ids2 <- modbase.gene.ids2[-(which(modbase.gene.ids2$modbase.gene.ids.X2 == "")), ] # drop the empty cells bc they're the ones in the first dataframe of same name 
colnames(modbase.gene.ids1) <- c("uniprot", "gene.id")
colnames(modbase.gene.ids2) <- c("uniprot", "gene.id")
gene.ids.to.uniprots <- rbind(modbase.gene.ids1, modbase.gene.ids2)
gene.ids.missing.in.modbase <- data.frame(uniq.gene.ids$unique.gene.id.[!(uniq.gene.ids$unique.gene.id. %in% gene.ids.to.uniprots$gene.id)]) #give me all the gene ids submitted to UniProt that returned no Uniprot

#gene.ids.missing.in.modbase <- data.frame(uniq.gene.ids$unique.gene.id.[!(uniq.gene.ids$unique.gene.id. %in% gene.ids.to.uniprots$gene.ids)])
#gene.ids.missing.in.modbase <- subset(gene.ids.missing.in.modbase, !(gene.ids.missing.in.modbase$gene.ids %in% $gene.id)) # list of modbase gene ids excluding 8 that are result of duplicate gene ids matching multiple uniprots

# figure out if data frame has NA
no.pdbs <- fread("./preProc_noPDB.txt", header=FALSE)
pdb.failed <- fread("./preProc_PDBFailed.txt", header=FALSE)
no.pdbs <- data.frame(unique(substr(no.pdbs$V1, start = 1, stop = 6))) #extract hfpd from hfpd+domain ids
pdb.failed <- data.frame(unique(substr(pdb.failed$V1, start = 1, stop = 6))) #extract hfpd from hfpd+domain ids

colnames(no.pdbs)[1] <- "hfpd"
colnames(pdb.failed)[1] <- "hfpd"
map.list <- read.table("./map_list", sep="\t", colClasses='character') #map_list is the file that maps btw hfpds and uniprots
map.list.uniprots <- data.frame(str_split_fixed(map.list$V2, ">", 2))
map.list$uniprots <- map.list.uniprots$X2
map.list$V2 <- NULL
colnames(map.list)[1] <- "hfpd"

# get me the uniprots of proteins with no pdb
uniprots.without.pdbs <- data.frame(map.list$uniprots[map.list$hfpd %in% no.pdbs$hfpd])

# get me the uniprots of proteins with pdbs that failed
uniprots.failed.pdbs <- data.frame(map.list$uniprots[map.list$hfpd %in% pdb.failed$hfpd])

# How many uniprots with no pdb (it's # of TRUE matches)
table(map.list$hfpd %in% no.pdbs$hfpd)

# How many uniprots failed (it's # of TRUE matches)
table(map.list$hfpd %in% pdb.failed$hfpd)

# total number of uniprots that could not be modeled (pdb failed + no pdbs)
nrow(uniprots.failed.pdbs) + nrow(uniprots.without.pdbs)

# dataframe of uniprots that could not be modeled (pdb failed + no pdbs)
colnames(uniprots.failed.pdbs) <- "uniprot"
colnames(uniprots.without.pdbs) <- "uniprot"
uniprots.not.modeled.by.us <- rbind(uniprots.failed.pdbs, uniprots.without.pdbs)

# discard all uniprots that are duplicated in the ModBase summary but absent from our (map.list)

# uniprots of all proteins not modeled by PrePPI that could be modeled by ModBase
gene.ids.to.uniprots.but.not.modeled.by.us <- data.frame(uniprots.not.modeled.by.us$uniprot[uniprots.not.modeled.by.us$uniprot %in% 
                                                                                                    gene.ids.to.uniprots$uniprot])
# uniprots of proteins not modeled by PrePPI (just due to failed pdbs) that could be modeled by ModBase
gene.ids.to.uniprots.but.failed.pdbs <- data.frame(uniprots.failed.pdbs$uniprot[uniprots.failed.pdbs$uniprot %in% 
                                                                                        gene.ids.to.uniprots$uniprot])
# uniprots of proteins not modeled by PrePPI (just due to no pdbs) that could be modeled by ModBase
gene.ids.to.uniprots.but.no.pdbs <- data.frame(uniprots.without.pdbs$uniprot[uniprots.without.pdbs$uniprot %in% 
                                                                                     gene.ids.to.uniprots$uniprot])


############################################################################################################
### successfully modeled uniprots
preppi.modeled.protein.list <- fread("./modeled_protein_list.txt", colClasses = 'character', header=FALSE)
colnames(preppi.modeled.protein.list) <- "hfpd"
preppi.modeled.uniprots <- data.frame(map.list$uniprots[map.list$hfpd %in% preppi.modeled.protein.list$hfpd])    
write.table(preppi.modeled.uniprots, "preppi.uniprots.txt", sep="\n",  # output it to text file
            row.names = FALSE, 
            col.names = FALSE,
            quote = FALSE)

write.table(gene.ids.to.uniprots$Entry, "modbase.uniprots.txt", sep="\n",  # output it to text file
            row.names = FALSE, 
            col.names = FALSE,
            quote = FALSE)

colnames(preppi.modeled.uniprots) <- "uniprots"
preppi.uniprots.overlapping.modbase.uniprots <- data.frame(preppi.modeled.uniprots$uniprots %in% gene.ids.to.uniprots$uniprot)
#preppi.uniprots.overlapping.modbase.uniprots <- data.frame(gene.ids.to.uniprots$uniprot %in% preppi.modeled.uniprots$uniprots)  
colnames(preppi.uniprots.overlapping.modbase.uniprots) <- "uniprots"

# VENN DIAGRAM
grid.newpage()
draw.pairwise.venn(nrow(preppi.modeled.uniprots), 
                   nrow(gene.ids.to.uniprots), 
                   nrow(gene.ids.to.uniprots), 
                   category = c("PrePPI Modeled Uniprots", "ModBase Model Uniprots"), 
                   lty = rep("blank", 2), 
                   fill = c("light blue", "pink"), 
                   alpha = rep(0.5, 1),
                   cat.cex = 1,
                   cat.pos = c(-10, 10))


############################################################################################################

# VENN DIAGRAM
grid.newpage()
draw.pairwise.venn(nrow(uniprots.without.pdbs), 
                   nrow(gene.ids.to.uniprots), 
                   nrow(gene.ids.to.uniprots.but.no.pdbs), 
                   category = c("Not PrePPI Modeled: No PDBs", "ModBase Model-able"), 
                   lty = rep("blank", 2), 
                   fill = c("light blue", "pink"), 
                   alpha = rep(0.5, 1),
                   cat.cex = 1,
                   cat.pos = c(-10, 10))

grid.newpage()
draw.pairwise.venn(nrow(uniprots.failed.pdbs), 
                   nrow(gene.ids.to.uniprots), 
                   nrow(gene.ids.to.uniprots.but.failed.pdbs), 
                   category = c("Not PrePPI Modeled: Failed PDBs", "ModBase Model-able"), 
                   lty = rep("blank", 2), 
                   fill = c("light green", "pink"), 
                   alpha = rep(0.5, 2),
                   cat.cex = 1,
                   cat.pos = c(-10, 10))

grid.newpage()
draw.pairwise.venn(nrow(uniprots.not.modeled.by.us), 
                   nrow(gene.ids.to.uniprots), 
                   nrow(gene.ids.to.uniprots.but.not.modeled.by.us), 
                   main = "Not Modeled vs Model-able by ModBase",
                   category = c("Overall Not PrePPI Modeled", "ModBase Modeled"), 
                   lty = rep("blank", 2), 
                   fill = c("orange", "pink"), 
                   alpha = rep(0.5, 2),
                   cat.cex = 1,
                   cat.pos = c(-10, 10))




