# Helpful Code for BCB 722 Topics in Population Genetics: Assignment 1

### Getting started ####
rm(list=ls())  # clear the working environment
setwd("/path/to/your/directory")  # point to the directory containing your data 

### Assignment Hints ####
# Most of the provided code must be modified and expanded to complete the assignment
# Make use of R internet resources and the help files for functions which can be viewed by entering "?howdoesthisfuncwork()" into the console, or by searching in RStudio's "help" pane.


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#


# 2.1) Read in a VCF file as a data frame ####
## This function will do the job:
read.vcf <- function(file, special.char="##", ...) {
  search.term=paste0(special.char, ".*")
  all.lines=readLines(file)
  clean.lines=gsub(search.term, "",  all.lines)
  clean.lines=gsub("#CHROM", "CHROM", clean.lines)
  read.table(..., text=paste(clean.lines, collapse="\n"))
}

data <- read.vcf(file="MyVCFfile.vcf", header=TRUE, stringsAsFactors=FALSE)  # call the function we just defined

## VCF files are formatted such that each row is a different variant or SNP, columns 1-9 contain info for each SNP, and columns 10 and up are each sampled genotypes for each SNP

dim(data)  # check the dimensions of the data.frame you just created
str(data)  # check out an overview of the data.frame you just created 

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#


# 2.2) Divide data into windows, create a results data frame, and count SNPs


## *USE THIS for the snm.vcf practice analysis
# This will create a simple results table with a single row, which represents one big genomic window with all the data

results <- data.frame(Start=min(data$POS), End=max(data$POS), SNPs=0)  # make a table (data.frame) to collect our results


## *UNCOMMENT AND USE THESE for your analysis of the 1000 genomes data

windows <- seq(min(data$POS), max(data$POS), by=5000)  # get window start positions *USE THIS ONE for your analysis of the 1000 genomes data
results <- data.frame(Start=windows, End=(windows+5000), SNPs=rep(0, length(windows)))  # Get Window End positions by adding 5000 to every Start position
## start with 0s in the SNP column, which will be replaced with the SNP counts you generate next

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#


## Counting SNPs (segregating sites) with a for loop:

for (i in 1:nrow(results)) { # Looping through each row of the results data.frame (windows), NOT the original input data
  d=subset(data, (data$POS>=results$Start[i] & data$POS<=results$End[i]))  # Get the subset of data within the window
  results$SNPs[i] <- nrow(d) # Count the number of SNPs in the window, add to results table
} #End for

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#


# 3) Calculating Summary Statistics: Do for EACH WINDOW you defined earlier - Note that for the snm data you only have one window, so you'll only have one of each result.

## Watterson's theta: 
### You'll need the number of chromosomes in the data sample (2N)

num.samples <- ncol(data)-9
results$Nchr <- 2*(num.samples)

#Enter your code to calculate Watterson's Theta here




#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#


## Pi

results$Pi=rep(0, nrow(results))  # an empty column of zeros to collect the the incoming Pi values
for (i in 1:nrow(results)) { # loop for every WINDOW
  d=subset(data, (data$POS>=results$Start[i] & data$POS<=results$End[i])) # get the subset of data in the window
  
  # Here, you should devise a way to get allele frequencies (counts) for each row in the subset of data (d). This vector of counts can be variable "j"
  # Helpful functions: strsplit(), as.numeric(), sum()
  # Calculate Pi using the "j" vector and the results$Nchr values.
  # Collect the Pi values in the results$Pi column
  
  #Enter your code to do the above here
  
} #End for

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#


## Tajima's D


### Use the variance function below, along with the values of Waterson's Theta, and Pi you calculated above to find Tajima's D for each window
# This variance function takes your number of chromosomes and number of segregating sites (SNPs) as input. 
variance.d <- function(n,S) {
  a1=sum(1/(seq(from=1, to=(n-1), by=1)))
  a2=sum(1/((seq(from=1, to=(n-1), by=1))**2))
  b1=(n+1)/(3*(n-1))
  b2=(2*((n**2)+n+3))/((9*n)*(n-1))
  c1=b1 - (1/a1)
  c2=b2-((n+2)/(a1*n)) + (a2/(a1**2))
  e1=c1/a1
  e2=c2/((a1**2)+a2)
  var=(e1*S) + (e2*S*(S-1))
  return(var)
} #End function


#Empty vector for Tajima's D values
results$TajimasD = rep(0, nrow(results))

#Code for calculating Tajima's D and inserting into results$TajimasD here



#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#


# 4) Plot Tajima's D vs. Genomic position ####
# It is ok if your Tajima's D values are all negative, SNPs in the 1000 genomes project seem to be rare
# However, they should be ~ around -3 to 1
plot(results$Start, results$TajimasD, pch=20,xlab="Position", ylab="Tajima's D", type="l", ylim=c(-3,1))

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#


## Add rectangles to the plot to annotate some genes

rect(left,bottom,right,top, col="carolina.blue")  # modify to use coordinates of your genes of interest



#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

### Acknowledgement: ####
# Elements of this assignment were inspired by Dr. Elizabeth Cooper at the University of North Carolina at Charlotte

