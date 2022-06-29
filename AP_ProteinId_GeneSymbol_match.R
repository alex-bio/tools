## Gene.Symbol to UniProt.Id matching script
## 2022-04-13 AP
## Someone sent you a list of proteins with only Gene.Symbol as an identifier. You recall that Gene.Symbol can often be interchangeable and may refer to specific isoforms, so some information may be lost. What do you do? Run this script to match informationless Gene.Symbol to a UniProt.Id! This script deals with duplicates...basically by not dealing with them (detail below). Anywhere there is EDIT:: written is where your participation is required.


#### Set your working directory ####
## This script sources files from your computer, but you have to tell it where those files are. Within the R window, click Session --> Set Working Directory --> Choose Directory. Choose where your files are located.

#### Loading your data ####
## EDIT:: first change the title of your file in quotes. Do not remove the quotes, and make sure it has .csv as the end. You can just "Save As" your Excel file as .csv. Your headers/column names BEFORE LOADING should be "Gene.Symbol" (without the quotes, capitalization is important). 
Gene.Symbols <- read.csv("hankum.csv",header=T)
## EDIT:: uploading your reference database. In this case, this is UniProt Human Genes (reviewed only) from 2022-04-10. If you are using what AP put in the Dropbox, then you don't have to change anything. Else, change the part below in quotes to the file name. Keep the quotes and the .csv. Your header/column names BEFORE LOADING should be "Protein.Id" for the UniProt ID number (6 alphanumeric symbols) and "Gene.Symbol" for the gene symbol column (capitalization is important). 
ref <- read.csv("uniprot-proteome_UP000005640.csv",header=T)

#### Data Processing; no edits ####
## packages you may need installed
list.of.packages <- c("stringi", "stringr","plyr","dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
## adding all functions of these packages to this session
library(stringi)
library(stringr)
library(plyr)
library(dplyr)
## only pull reviewed UniProt entries
ref <- subset(ref, ref$Status=="reviewed")
## eliminate the semicolons and split up any genes that are written "geneA/geneA1"
ref$Gene.Symbol <- str_replace_all(ref$Gene.Symbol, c("/" = " ", ";" = ""))
## now split each Gene.Symbol entry such that each gene name is its own character string (right now, they are one long character string)
ref$Gene.Symbol <- str_split(ref$Gene.Symbol, " ")
## only need first two columns
ref <- ref[c("Protein.Id","Gene.Symbol")]
## Now change your own Gene.Symbols list into a vector for next steps.
Gene.Symbol.vec <- as.vector(Gene.Symbols$Gene.Symbol)

## making a list for the upcoming for loop
list.precursor <- vector("list",nrow(ref))
## this loop is some spaghetti ass code right here...but it works
for (i in 1:nrow(ref)){
  y <- ref[i, c("Gene.Symbol")]
  list.precursor[[i]] <- intersect(unlist(y), Gene.Symbol.vec)
}
## replace empty/non-matches with NA for filtering out
list.precursor <- lapply(list.precursor, function(x) if(identical(x, character(0))) NA_character_ else x)
## collapse everything
## This step will create a new column "keywordtag" for each matching Gene.Symbol. *If you have the same gene that goes by two Gene.Symbols, and both of these Gene.Symbols are in your Gene.Symbol list (for example, GeneA also goes by GeneB and both GeneA and GeneB are in your list), then list.precursor will contain list elements that are greater than zero. These next steps pull those repeated elements out. If there are values here, then only the first Gene.Symbol matched will show up in your final matched list.
## first dealing with this^ case
same.gene.n.times <- Filter(function(x) length(x) > 1, list.precursor)
same.gene.n.times <- as.data.frame((as.matrix(do.call(rbind.data.frame, same.gene.n.times))))
dummy <- ncol(same.gene.n.times)
colnames(same.gene.n.times) <- c("ThisGeneSameAs","This",rep(paste("AndThis"),dummy-2))
file_name4 <- paste0(format(Sys.time(), "%Y-%m-%d_%H%M%S_"), "synonymous_GeneSymbols.csv")
write.csv(same.gene.n.times, file_name4)
## now pasting matches to the reference/database dataframe, but only first column
keywordtag.precursor <- do.call(rbind.data.frame, list.precursor)
ref$keywordtag <- keywordtag.precursor[1]

## pull out the matches
matches <- subset(ref, !is.na(ref$keywordtag))
## listed lists, so unlist to find duplicates. A way to think about this is that for every UniProt entry, if one of the Gene.Symbol matches, then it is pasted as a pair. All matching pairs are pulled out, so 1 Gene.Symbol may match multiple UniProt.Id. 
matches.v <- unlist(matches$keywordtag)
## how many duplicates? 
sum(matches.v %in% matches.v[duplicated(matches.v)])
## what are the duplicates? This writes the duplicates to a .csv file in your designated directory. This ain't no fancy AI that can magically assigned duplicates to the correct UniProt entry.
matches.table <- data.frame(table(matches$keywordtag))
duplicates <- matches.table[matches.table$Freq >1,]
names(duplicates) <- c("Gene.Symbol", "Appears.n.times")
file_name <- paste0(format(Sys.time(), "%Y-%m-%d_%H%M%S_"), "duplicates.csv")
write.csv(duplicates, file_name)
## The duplicates are carried into the final matched list, but highlighted as separate .csv file. If you want to manually annotate duplicates later, do it after the code runs.
## this "matches.final" dataframe matched the "keywordtag", which is from the original list you fed in to the Uniprot Id found in the Uniprot database. "keywordtag" is renamed to "Gene.Symbol" below.
## this dataframe gets rid of the column of all possible gene names from UniProt (keeps the one you input)
matches.final <- as.data.frame(as.matrix(matches[,c("Protein.Id", "keywordtag" )]))
colnames(matches.final) <- c("Protein.Id","Gene.Symbol")
file_name2 <- paste0(format(Sys.time(), "%Y-%m-%d_%H%M%S_"), "ProteinId_GeneSymbol_matches.csv")
write.csv(matches.final, file_name2)

## also I guess a list to show you any Gene.Symbols that weren't matched to the database.
no.match.final <- matches.final$Gene.Symbol
no.match.original.list <- Gene.Symbols$Gene.Symbol
## matching the final matched Gene.Symbol against what you input in your orignal list; "hankum.csv" is the example
no.match <- setdiff(no.match.original.list, no.match.final)
file_name3 <- paste0(format(Sys.time(), "%Y-%m-%d_%H%M%S_"), "not_matched.csv")
write.csv(no.match, file_name3)

## boom baby
## ta-da!