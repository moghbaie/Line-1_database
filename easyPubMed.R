# Mehrnoosh Oghbaie
# Date: 06/07/2018
# Querying keywords to find pubmed ID

library(easyPubMed)
library(rcrossref)
library(httr) #for the BROWSE() function
library(XML)
library(rentrez)
library(dplyr)
library(magrittr)


#searching for related papers
setwd("C:/Users/moghb/OneDrive/Documents/Line1")

my_query0 <- "LINE1 OR L1-Retrotransposons OR Line-1 OR Line_1 OR LINE1" #"3389"
my_query <- "LINE1 OR L1-Retrotransposons OR Line-1  OR Line_1 OR LINE1 OR long interspersed nuclear element OR long interspersed element" #"3731"

getPubmed_ID <- function(my_query){
  my_entrez_id <- get_pubmed_ids(my_query)
  my_abstracts_txt <- fetch_pubmed_data(my_entrez_id, format= "abstract", retstart=0, retmax=as.integer(as.character(my_entrez_id$Count)))
  pubmed_IDlist <- as.integer(
    gsub("\\PMID: |\\s|\\[Indexed for MEDLINE]","", 
         grep(pattern = "PMID:", my_abstracts_txt, value = TRUE)
    )
  )
  print(complete <- ifelse(length(pubmed_IDlist) == my_entrez_id$Count, "Complete", "missing records"))
  return(pubmed_IDlist)
}


pubmed_IDlist <- getPubmed_ID(my_query)
pubmed_IDlist_homo <- getPubmed_ID(paste(my_query,"AND homo"))
LINE1 <- getPubmed_ID("LINE1 OR Line-1 OR Line_1")
L1_Retrotransposons <- getPubmed_ID("L1-Retrotransposons OR L1_Retrotransposons")
long_interspersed_nuclear_element <- getPubmed_ID("long interspersed nuclear element")
long_interspersed_element <- getPubmed_ID("long interspersed element")


PMCID <- lapply(pubmed_IDlist, function(x) id_converter(x,"pmid")$records$pmcid)
PMCID[PMCID=="NULL"] <- NA
sortpubdate <- do.call(rbind, 
                       lapply(pubmed_IDlist, function(i) entrez_summary(db="pubmed", id=i)$sortpubdate)
)


paper_list <- data.frame("pubMedID"=pubmed_IDlist,
                         "sortpubdate" = as.Date(sortpubdate, "%Y/%m/%d"),
                         "PMCID"= unlist(PMCID),
                         "pubMedID_homo" = rep(NA, length(pubmed_IDlist)),
                         "LINE1"= rep(NA, length(pubmed_IDlist)),
                         "L1_Retrotransposons"= rep(NA, length(pubmed_IDlist)),
                         "long_interspersed_nuclear_element" = rep(NA, length(pubmed_IDlist)),
                         "long_interspersed_element" = rep(NA, length(pubmed_IDlist))
                         )

paper_list$pubMedID_homo[!is.na(match(pubmed_IDlist, pubmed_IDlist_homo))] <- 1
paper_list$LINE1[!is.na(match(pubmed_IDlist, LINE1))] <- 1
paper_list$L1_Retrotransposons[!is.na(match(pubmed_IDlist, L1_Retrotransposons))] <- 1
paper_list$long_interspersed_nuclear_element[!is.na(match(pubmed_IDlist, long_interspersed_nuclear_element))] <- 1
paper_list$long_interspersed_element[!is.na(match(pubmed_IDlist, long_interspersed_element))] <- 1


write.csv(paper_list, file = "Output/papers.csv")

paper <- paper_list%>%
  select("sortpubdate")%>%
  format("%Y")
paper <- data.frame(table(paper))

paper_homo <- paper_list%>%
  filter(pubMedID_homo==1)%>%
  select("sortpubdate")%>%
  format("%Y")
paper_homo <- data.frame(table(paper_homo))

paper$homo <- match(paper$paper,paper_homo$paper_homo)
paper[is.na(paper)]=0
colnames(paper) <- c("Year","Non-homo_sapiens",'homo_sapiens')
rownames(paper)<- paper$paper

jpeg('Image/line1_papers_distribution.jpg')
barplot(as.matrix(t(paper[,-1])), col=c("orange","red"),
        main="Distribution of Line1 papers during times",
        legend.text = TRUE, args.legend = list(x = "top"))

dev.off()


