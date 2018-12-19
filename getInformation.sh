#!/bin/bash
# Copyright @2007-2019 Free Software Foundation

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

# The getTaxonNumber function.
# This function sets a integer based on the chosen taxon level.
getTaxonNumber() {
    if [ "${disTaxon}" == "Kingdom" ]
    then
        disPrune="1"
    fi
    if [ "${disTaxon}" == "Phylum" ]
    then
        disPrune="2"
    fi
    if [ "${disTaxon}" == "Class" ]
    then
        disPrune="3"
    fi
    if [ "${disTaxon}" == "Order" ]
    then
        disPrune="4"
    fi
    if [ "${disTaxon}" == "Family" ]
    then
        disPrune="5"
    fi
    if [ "${disTaxon}" == "Genus" ]
    then
        disPrune="6"
    fi
    if [ "${disTaxon}" == "Species" ]
    then
        disPrune="7"
    fi
}

# The getFormatFlow function.
# This function creates a temporary storage directory in the output directory.
# It then calls the getInformation.R script with the correct input values.
# After the script is finished, output is send to the expected location and
# the temporary storage directory is deleted.
getFormatFlow() {
    getTaxonNumber
    strDirectory=${fosOutput::-4}
    mkdir -p "${strDirectory}_temp"
    getInformation.R -i ${fisInput} -o ${strDirectory}_temp/flNewOutput.pdf \
                     -w ${strDirectory}_temp -f ${disFormat} \
                     -t ${disTaxon} -p ${disPrune} -a ${fosAbsoluteBar} \
                     -b ${fosRelativeBar} -r ${fosRichness} \
                     -e ${fosAbsoluteMap} -m ${fosRelativeMap} 2>&1 >/dev/null
    cat ${strDirectory}_temp/flNewOutput.pdf > ${fosOutput}
    rm -rf ${strDirectory}_temp
}

# The main function.
main() {
    getFormatFlow
}

# The getopts function.
while getopts ":i:o:f:t:a:b:r:e:m:vh" opt; do
    case ${opt} in
        i)
            fisInput=${OPTARG}
            ;;
        o)
            fosOutput=${OPTARG}
            ;;
        f)
            disFormat=${OPTARG}
            ;;
        t)
            disTaxon=${OPTARG}
            ;;
        a)
            fosAbsoluteBar=${OPTARG}
            ;;
        b)
            fosRelativeBar=${OPTARG}
            ;;
        r)
            fosRichness=${OPTARG}
            ;;
        e)
            fosAbsoluteMap=${OPTARG}
            ;;
        m)
            fosRelativeMap=${OPTARG}
            ;;
        v)
            echo ""
            echo "getInformation.sh [0.1.0]"
            echo ""

            exit
            ;;
        h)
            echo ""
            echo "Usage: getInformation.sh [-h] [-v][-i INPUT] [-o OUTPUT]"
            echo "                         [-f FORMAT] [-t TAXON]"
            echo "                         [-a ABSOLUTEBAR] [-b RELATIVEBAR]"
            echo "                         [-r RICHNESS] [-e ABSOLUTEMAP]"
            echo "                         [-m RELATIVEMAP]"
            echo ""
            echo "Optional arguments:"
            echo " -h                    Show this help page and exit"
            echo " -v                    Show the software's version number"
            echo "                       and exit"
            echo " -i                    The location of the input file(s)"
            echo " -o                    The location of the output file(s)"
            echo " -f                    The format of the input"
            echo "                       file(s) [otu_old/otu_new/lca]"
            echo " -t                    The taxon level for graphic output"
            echo "                       [Kingdom/Phylum/Class/Order/Family/"
            echo "                       Genus/Species]"
            echo " -a                    The location of the absolute bar"
            echo "                       output file(s)"
            echo " -b                    The location of the relative bar"
            echo "                       output file(s)"
            echo " -r                    The location of the richness plot"
            echo "                       output file(s)"
            echo " -e                    The location of the absolute heatmap"
            echo "                       output file(s)"
            echo " -m                    The location of the relative heatmap"
            echo "                       output file(s)"
            echo ""
            echo "The Statistical Analysis tool will utilize the Phyloseq R"
            echo "package to create multiple plots based on a OTU table."
            echo ""
            echo "These plots are:"
            echo " - A bar plot: a plot to quickly and easily show informative"
            echo "   summary graphics of the differences in taxa abundance"
            echo "   between samples in an experiment."
            echo " - A richness plot: a number of alpha-diversity metrics"
            echo "   (Observed, Shannon and Simpson)"
            echo " - A heatmap plot: a heat map is a false color image with a"
            echo "   dendrogram added to the top."
            echo ""
            echo "All of these plots will be generated using both the absolute"
            echo "numbers in a OTU table as well as relative numbers for a OTU"
            echo "table."
            echo ""
            echo "Sample names can not start with a '#'."
            echo "All columns in a OTU table should have a header starting"
            echo "with '#'."
            echo ""
            echo "Source(s):"
            echo " - McMurdie PJ, Holmes S, Phyloseq: An R Package for"
            echo "   Reproducible Interactive Analysis and Graphics of"
            echo "   Microbiome Census Data."
            echo "   PLOS One. 2013; 8(4)."
            echo "   doi: 10.1371/journal.pone.0061217"
            echo "   https://joey711.github.io/phyloseq/"
            echo ""

            exit
            ;;
        \?)
            echo ""
            echo "You've entered an invalid option: -${OPTARG}."
            echo "Please use the -h option for correct formatting information."
            echo ""

            exit
            ;;
        :)
            echo ""
            echo "You've entered an invalid option: -${OPTARG}."
            echo "Please use the -h option for correct formatting information."
            echo ""

            exit
            ;;
    esac
done

main

# Additional information:
# =======================
#
# Sample names can not start with a "#".
# All columns in a OTU table should have a header starting with "#".
