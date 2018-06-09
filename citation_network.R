# Mehrnoosh Oghbaie
# Date: 06/07/2018
# Build citation network

require(plyr)
require(XML)
require(RCurl)
require(readxl)
require(utils)
require(igraph)

setwd("C:/Users/moghb/OneDrive/Documents/Line1")
# We start with a list PubMed ids, for which we would like to have the citation graph

papers <- read_csv("~/Line1/Output/papers.csv", col_types = cols(X1 = col_skip()))

citegraph <- NULL
for(pmid in papers$pubMedID){
  # Get more info with the pmid id
  tryCatch({
    
    # We build the URL to retrieve the citing pmids
    path <- paste("https://www.ncbi.nlm.nih.gov/pubmed?linkname=pubmed_pubmed_citedin&from_uid=",
                  pmid,
                  "&report=uilist&format=text&dispmax=2000", sep="")
    f <- file(path)
    data <- readLines(f, warn = F)
    citing <-strsplit(xmlToList(xmlParse(data, asText = T))[1], "\n")[[1]]
    close(f)
    
    # Then we create a table containing ll these pmids
    if(length(citing) < 2000){
      for(pm in citing){
        citegraph <- rbind(citegraph, data.frame(pmid=as.numeric(pmid), citing=as.numeric(pm)))
      }
    }
    
    message(paste0(pmid," done: ",length(citing)," citations"))
    
  }, warning = function(w) {
    message(w)
  }, error = function(e) {
    message(e)
  })
}

#write.csv(citegraph, file= "Output/citation_network.csv")
citation_network <- read_csv("~/Line1/Output/citation_network.csv", 
                             col_types = cols(X1 = col_skip()))
colnames(citation_network)<- c("from","to")
citation_network_L1 <- citation_network[citation_network$to %in% papers$pubMedID,]
citation_network_L1$citedate <- format(papers$sortpubdate[match(citation_network_L1$to, papers$pubMedID)],"%Y")
net <- graph_from_data_frame(citation_network_L1, directed = TRUE)

summary(degree(net, mode = c("out")))
summary(degree(net, mode = c("in")))

jpeg('Image/citation_distribution.jpg')
hist(degree(net, mode = c("out")), breaks=20, col="blue",lty="blank",
     main="Histagram of citations of Line-1 papers", xlab="Number of citations")
legend("topright", inset=.05, title="Summary of citation among Line-1 papers",
       legend=c(paste0("Max: ",summary(degree(net, mode = c("out")))["Max."]),
                paste0("Mean: ",round(summary(degree(net, mode = c("out")))["Mean"],3)),
                paste0("Median: ",summary(degree(net, mode = c("out")))["Median"]),
                paste0("Min: ",summary(degree(net, mode = c("out")))["Min."])), 
       horiz=FALSE)
dev.off()

cite_pmid <- data.frame(degree(net, mode = c("out")))
colnames(cite_pmid) <- c("degree_out")
cite_pmid$pubdate <- format(
  papers$sortpubdate[match(rownames(cite_pmid),as.character(papers$pubMedID))],
  "%Y")

colfunc<-colorRampPalette(c("white","red"))
barplot(as.matrix(table(cite_pmid)), col = (colfunc(20)))
