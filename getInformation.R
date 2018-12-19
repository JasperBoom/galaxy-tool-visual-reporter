#!/usr/bin/Rscript
# Copyright @The R Foundation.

# Author: Jasper Boom

# Prequisites:
# - sudo apt-get install r-base
# - sudo chmod -R 777 usr/lib/R
# - sudo chmod -R 777 usr/share/R
# - sudo apt-get install libcurl4-gnutls-dev
# - sudo apt-get install libssl-dev
# - source("http://bioconductor.org/biocLite.R")
# - biocLite("phyloseq")
# - biocLite("optparse")

# Galaxy prequisites:
# - sudo ln -s /path/to/folder/galaxy-tool-visual-reporter/getInformation.R
#              /usr/local/bin/getInformation.R

library("optparse")
library("ggplot2")
library("phyloseq")
library("cluster")
library("methods")

# The getFigures function.
# This function creates the phyloseq plots using the phyloseq OTU object.
# It accumulates all OTUs up to a user specified taxon level.
# It creates a new phyloseq object with relative abundances by dividing all
# numbers by their total. These two objects are then used to create bar plots,
# alpha diversity plots and heatmaps. All of these plots are outputted to a
# pdf file. Additionally, the plots will be outputted to seperate postscript
# files.
getFigures <- function(objPhyseq, flOutput, strTaxon, intPrune, flAbsoluteBar,
                       flRelativeBar, flRichnessPlot, flAbsoluteHeatmap,
                       flRelativeHeatmap){
    objPhyseq = tax_glom(objPhyseq, strTaxon)
    setTaxa = tax_table(objPhyseq)[1:length(taxa_names(objPhyseq)),
                                   as.integer(intPrune)]
    objRelative = transform_sample_counts(objPhyseq, function(x) x / sum(x))
    plotAbsoluteBar <- plot_bar(objPhyseq, fill = strTaxon)
    plotRelativeBar <- plot_bar(objRelative, fill = strTaxon)
    plotRichness <- plot_richness(objPhyseq, measures=c("Observed", "Shannon",
                                                        "Simpson"))
    pdf(flOutput)
    print(plotAbsoluteBar + ggtitle("Absolute abundance bar plot"))
    print(plotRelativeBar + ggtitle("Relative abundance bar plot"))
    print(plotRichness + ggtitle("Alpha diversity plot(s)"))
    heatmap(otu_table(objPhyseq), main = "Absolute abundance heatmap",
            Rowv=NA, labRow = setTaxa)
    heatmap(otu_table(objRelative), main = "Relative abundance heatmap",
            Rowv=NA, labRow = setTaxa)
    postscript(flAbsoluteBar)
    print(plotAbsoluteBar + ggtitle("Absolute abundance bar plot"))
    postscript(flRelativeBar)
    print(plotRelativeBar + ggtitle("Relative abundance bar plot"))
    postscript(flRichnessPlot)
    print(plotRichness + ggtitle("Alpha diversity plot(s)"))
    postscript(flAbsoluteHeatmap)
    heatmap(otu_table(objPhyseq), main = "Absolute abundance heatmap",
            Rowv=NA, labRow = setTaxa)
    postscript(flRelativeHeatmap)
    heatmap(otu_table(objRelative), main = "Relative abundance heatmap",
            Rowv=NA, labRow = setTaxa)
}

# The getTaxaTableLca function.
# This function creates a matrix from every taxonomy column in a LCA file.
# These 7 matrices are combined into one matrix with 7 columns. The OTU
# names and taxonomy level names are added to this matrix. And the phyloseq
# function creates a taxonomy table from that matrix. The new phyloseq taxon
# object is returned.
getTaxaTableLca <- function(dfOtuTable, dfOtuTableWithRowNames,
                            intColCount, dfOtuTableTrimmedBlast){
    mtxKingdom = as.matrix(dfOtuTable[intColCount+5])
    mtxPhylum = as.matrix(dfOtuTable[intColCount+6])
    mtxClass = as.matrix(dfOtuTable[intColCount+7])
    mtxOrder = as.matrix(dfOtuTable[intColCount+8])
    mtxFamily = as.matrix(dfOtuTable[intColCount+9])
    mtxGenus = as.matrix(dfOtuTable[intColCount+10])
    mtxSpecies = as.matrix(dfOtuTable[intColCount+11])
    intNumberOfRows = nrow(dfOtuTableWithRowNames)
    mtxTaxa = matrix(c(mtxKingdom, mtxPhylum, mtxClass,
                       mtxOrder, mtxFamily, mtxGenus,
                       mtxSpecies), nrow=intNumberOfRows, ncol=7)
    vcrTaxa = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus",
                "Species")
    colnames(mtxTaxa) = vcrTaxa
    rownames(mtxTaxa) = dfOtuTableTrimmedBlast[,1]
    dfTax = tax_table(mtxTaxa)
    return(dfTax)
}

# The getTaxaTableOtu function.
# This function loops through the taxonomy column and seperates the taxonomic
# names based on the slash(/) character per row. Each name is put into a vector
# based on their position. These positions correlate to their taxonomic level.
# A new matrix is created with a column for every taxonomic level. This new
# matrix is then used to create a phyloseq taxonomic table. This table is
# returned.
getTaxaTableOtu <- function(strLastColumn, dfOtuTableWithRowNames,
                         dfOtuTableTrimmedBlast){
    vcrIdentifications = c()
    for (strLine in strLastColumn){
        strTemp = strsplit(strLine, " / ")
        for (intPosition in c(1, 2, 3, 4, 5, 6 , 7)){
            vcrIdentifications = c(vcrIdentifications,
                                   strTemp[[1]][intPosition])
        }
    }
    intNumberOfRows = nrow(dfOtuTableWithRowNames)
    mtxTaxa = matrix(vcrIdentifications, ncol = 7, nrow = intNumberOfRows,
                     byrow = TRUE)
    vcrTaxa = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus",
                "Species")
    colnames(mtxTaxa) = vcrTaxa
    rownames(mtxTaxa) = dfOtuTableTrimmedBlast[,1]
    dfTax = tax_table(mtxTaxa)
    return(dfTax)
}

# The getPhyloseqTable function.
# This function sets a new working directory and loads in the input file.
# The headers in the input file are selected and the sample file names are
# isolated. Based on the number of samples, the taxonomy column is isolated.
# The BLAST identifications are split off and the OTU identifiers are set as row
# indexes. This altered OTU table is transformed into a matrix and then into a
# Phyloseq OTU table. The BLAST identifications are send to the getTaxaTableOtu
# function. The result from that function and the Phyloseq OTU table are
# combined into a Phyloseq object. This object is send to the getFigures
# function.
getPhyloseqTable <- function(strWorking, flInput, flOutput, strFormat, strTaxon,
                             intPrune, flAbsoluteBar, flRelativeBar,
							 flRichnessPlot, flAbsoluteHeatmap,
							 flRelativeHeatmap){
    setwd(strWorking)
    dfOtuTable = read.csv(flInput, header = TRUE, sep="\t",
                          check.names = FALSE)
    strColNames = colnames(dfOtuTable)
    intColCount = 0
    for (strName in strColNames){
        if (strName != "" & strName != "OccurrenceStatus" &
            substring(strName, 1, 1) != "#"){
            intColCount = intColCount + 1
        }
    }
    if (strFormat == "otu_old"){
        strLastColumn = dfOtuTable[,intColCount+12]
    } else if (strFormat == "otu_new"){
        strLastColumn = dfOtuTable[,intColCount+11]
    }
    dfOtuTableTrimmedBlast = dfOtuTable[,c(0:intColCount+1)]
    dfOtuTableWithRowNames = data.frame(dfOtuTableTrimmedBlast[,-1],
                                        row.names=dfOtuTableTrimmedBlast[,1],
                                        check.names = FALSE)
    mtxOtuTable = data.matrix(dfOtuTableWithRowNames)
    dfOtuPhyloseq = otu_table(mtxOtuTable, taxa_are_rows = TRUE)
    if (strFormat == "otu_old" | strFormat == "otu_new"){
        dfTaxonomy = getTaxaTableOtu(strLastColumn, dfOtuTableWithRowNames,
                                     dfOtuTableTrimmedBlast)
    } else if (strFormat == "lca"){
        dfTaxonomy = getTaxaTableLca(dfOtuTable, dfOtuTableWithRowNames,
                                     intColCount, dfOtuTableTrimmedBlast)
    }
    objPhyseq = phyloseq(dfOtuPhyloseq, dfTaxonomy)
    getFigures(objPhyseq, flOutput, strTaxon, intPrune, flAbsoluteBar,
               flRelativeBar, flRichnessPlot, flAbsoluteHeatmap,
			   flRelativeHeatmap)
}

# The argvs function.
parseArgvs <- function(){
    parser <- OptionParser()
    parser <- add_option(parser, c("-v", "--version"), action = "store_true",
                         default = FALSE,
                         help="Show program's version number and exit")
    parser <- add_option(parser, c("-w", "--fisWorking"), default = "NONE",
                         help="The location of the working directory")
    parser <- add_option(parser, c("-i", "--fisInput"), default = "NONE",
                         help="The location of the input file(s)")
    parser <- add_option(parser, c("-o", "--fosOutput"), default = "NONE",
                         help="The location of the output file(s)")
    parser <- add_option(parser, c("-f", "--disFormat"), default = "NONE",
                         help="The format of the input file(s) [otu_old/otu_new/lca]")
    parser <- add_option(parser, c("-t", "--disTaxon"), default = "NONE",
                         help="The taxon level for graphic output [Kingdom/Phylum/Class/Order/Family/Genus/Species]")
    parser <- add_option(parser, c("-p", "--disPrune"), default = "NONE",
                         help="The taxon level number for graphic output [1/2/3/4/5/6/7]")
    parser <- add_option(parser, c("-a", "--fosAbsoluteBar"), default = "NONE",
                         help="The location of the absolute bar output file(s)")
    parser <- add_option(parser, c("-b", "--fosRelativeBar"), default = "NONE",
                         help="The location of the relative bar output file(s)")
    parser <- add_option(parser, c("-r", "--fosRichness"), default = "NONE",
                         help="The location of the richness plot output file(s)")
    parser <- add_option(parser, c("-e", "--fosAbsoluteMap"), default = "NONE",
                         help="The location of the absolute heatmap output file(s)")
    parser <- add_option(parser, c("-m", "--fosRelativeMap"), default = "NONE",
                         help="The location of the relative heatmap output file(s)")
    return(parser)
}

# The main function.
main <- function(){
    argvs = parse_args(parseArgvs())
    if ( argvs$version ){
        print("getInformation.R [0.1.0]")
        quit()
    }
    if ( argvs$fisInput != "NONE" & argvs$fosOutput != "NONE" & 
         argvs$fisWorking != "NONE" & argvs$disFormat != "NONE" &
         argvs$disTaxon != "NONE" & argvs$disPrune != "NONE" &
         argvs$fosAbsoluteBar != "NONE" & argvs$fosRelativeBar != "NONE" &
         argvs$fosRichness != "NONE" & argvs$fosAbsoluteMap != "NONE" &
         argvs$fosRelativeMap != "NONE"){
        getPhyloseqTable(argvs$fisWorking, argvs$fisInput, argvs$fosOutput,
                         argvs$disFormat, argvs$disTaxon, argvs$disPrune,
                         argvs$fosAbsoluteBar, argvs$fosRelativeBar,
                         argvs$fosRichness, argvs$fosAbsoluteMap,
                         argvs$fosRelativeMap)
    }
}

main()

# Additional information:
# =======================
#
# File names can not start with a "#".
# All columns in a OTU table should have a header starting with "#".
