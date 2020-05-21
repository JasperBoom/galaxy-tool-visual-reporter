# galaxy-tool-visual-reporter

__!!Repo can be deleted, copy is also kept [here](https://github.com/JasperBoom/galaxy-tools-naturalis-internship)!!.__

Use Phyloseq to a statistical analysis resulting in multiple plots.  
The Statistical Analysis tool will utilize the Phyloseq R package to create multiple plots based on a OTU table.

Sample names can not start with a "#".  
All columns in a OTU table should have a header starting with "#".

# Getting started

### Prerequisites
Download and install the following software according to the following steps.
```
sudo apt-get install r-base
sudo apt-get install libcurl4-gnutls-dev
sudo apt-get install libssl-dev
```
These packages should be installed in R.
```
source("http://bioconductor.org/biocLite.R")
biocLite("phyloseq")
biocLite("optparse")
```

### Installing
Download and install the tool according to the following steps.
```
sudo mkdir -m 755 /home/Tools
cd /home/Tools
sudo git clone https://github.com/JasperBoom/galaxy-tool-visual-reporter
sudo chmod -R 755 galaxy-tool-visual-reporter
```
The following file in the galaxy-tool-visual-reporter folder should be made avaible from any location.
```
sudo ln -s /home/Tools/galaxy-tool-visual-reporter/getInformation.R /usr/local/bin/getInformation.R
```
Continue with the tool installation
```
sudo mkdir -m 755 /home/galaxy/tools/directoryname
sudo cp /home/Tools/galaxy-tool-visual-reporter/getInformation.sh /home/galaxy/tools/directoryname/getInformation.sh
sudo cp /home/Tools/galaxy-tool-visual-reporter/getInformation.xml /home/galaxy/tools/directoryname/getInformation.xml
```
Edit the following file in order to make galaxy display the tool.
```
/home/galaxy/config/tool_conf.xml
```
```
<tool file="airdentification/getInformation.xml"/>
```

## Source(s)
* __McMurdie PJ, Holmes S__,  
  Phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data.  
  PLOS One. 2013; 8(4). __doi: 10.1371/journal.pone.0061217__  
  [Phyloseq](https://joey711.github.io/phyloseq/)
* __Giardine B, Riemer C, Hardison RC, Burhans R, Elnitski L, Shah P__,  
  Galaxy: A platform for interactive large-scale genome analysis.  
  Genome Research. 2005; 15(10) 1451-1455. __doi: 10.1101/gr.4086505__  
  [Galaxy](https://www.galaxyproject.org/)

```
Copyright (C) 2018 Jasper Boom

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License version 3 as
published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
```
