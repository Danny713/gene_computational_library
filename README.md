# gene_computational_library
This is a project developed as my Junior Independent Work

This library requires R (https://www.r-project.org/). 

Install the following R packages:
  install.packages("data.table")
  install.packages("rJava")
  install.packages("pryr")
  source("http://bioconductor.org/biocLite.R")
  biocLite("WGCNA")

To use the command-line client provided:
Compilation:
    javac -cp Import/*:. CorrelationClient.java

Execution:
    ./run -cp Import/*:. CorrelationClient
    
Some of the logic are extracted from a project called GenEx Project.
