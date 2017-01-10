# Some of the methods are from the GenEx Project - written by Rachel Xu '17

###source("queue.R")

library(pryr)
library(rJava)
.jinit() # this starts the JVM

.jaddClassPath("/Users/home_folder/Desktop/Princeton/Junior/Independent\ Work/source_file_original")
.jclassPath()
#dyn.load("../src/pvalAndTriFunctions.so")

testinggggg <-function() {
	return("ok")
}

loadAggDT <- function(filePath) {
	print("starting new impl")
	aggDT <- read.csv(
		file=filePath, 
		header=TRUE, 
		sep=",")
	aggDT <- data.table(aggDT)

	return(aggDT)
}


corrData <- function(
    aggDT,
    P = 0.05,
    tau = 0.9
    )
{
	edgeDF <- data.frame(index1=integer(), index2=integer(), pVal=double()) # create an empty dataframe
    corrEnv = tryCatch(
    { #correlation
	    corrEnv <- new.env()
	    corrEnv$index <- as.character(aggDT$Gene)

	    x <- t(aggDT[, Gene := NULL])
	    rownames(x) <- NULL


	    NCOL <- ncol(x)
	    DFR <- (nrow(x) - 2)
	    LEN <- (NCOL * (NCOL - 1))/2 

	    ## preallocate vectors and matrices of dimension into new environment
	   
	    corrEnv$pVec <- vector(mode = 'numeric', length =  LEN)
	    corrEnv$corrVec <- vector(mode = 'numeric', length = LEN)

	    #Data structures used for BFS:
	    visited <- logical(length = NCOL) #boolean array (visited) of size NCOL
	    queue <- new.queue()
	    visitCounter = 0
	    lastCheckedIndex = 1 # R begins index at 1

	   	enqueue(queue, 1)
	    while (visitCounter < NCOL) {
	    	isPriority = FALSE
	    	if ((priorityGene = getPriorityGene(corrEnv$index)) > 0 && !visited[priorityGene]) {
	    		#empty queue and insert into queue priorityGene
	    		empty(queue)
	    		enqueue(queue, priorityGene)
	    		isPriority = TRUE
	    	} else if (is.empty(queue)) {
	    		#enqueue a new gene to start
	    		for (i in lastCheckedIndex:NCOL) {
	    			if (!visited[i]) {
	    				lastCheckedIndex = i + 1
	    				enqueue(queue, i)
	    				break
	    			}
	    		}
	    	}

			newGene <- dequeue(queue)
			if (!visited[newGene]) {
				visited[newGene] = TRUE
				visitCounter <- visitCounter + 1

				for (i in 1:NCOL) {
					if (!visited[i]) {
						COL <- spearman(x[, newGene], x[, i])
						pVal <- getp2(COL, DFR)

						vecIndex <- getIndex(newGene, i, NCOL)
						# add COL to corrEnv$corrVec
						corrEnv$corrVec[vecIndex] <- COL[, 1]
						# add pVAL to corrEnv$pVec
						corrEnv$pVec[vecIndex] <- pVal

						if (pVal <= P && abs(COL[, 1]) >= tau) {
				        	# add to edge data frame
				        	edgeDF[nrow(edgeDF) + 1,] <- c(newGene, i, COL[, 1])
				        	#append to queue
				        	enqueue(queue, i)
				        }
					}
				}
				# Update user if gene is a priority gene
				if (isPriority && visitCounter < NCOL) {
					sendUpdatedEdgeList(edgeDF, corrEnv$index, FALSE)
				}
			}
	    }
	    sendUpdatedEdgeList(edgeDF, corrEnv$index, TRUE)

	    return(corrEnv)

	}, error = function(e) {
		print(e)
		#return error
		return(-1)
	}, finally = {
		#clean up
		gc()
	}) #END tryCatch
}

getIndex <- function(first, second, N) {
	if (first < second) {
		#first is col, second is row
		return((N*(first-1)) - ((first*(first-1))/2) + (second - first))
	} else {
		#second is col, first is row
		return((N*(second-1)) - ((second*(second-1))/2) + (first - second))
	}
}

# output any matrix vectors + index to .csv file
# works for *SYMMETRICAL* matrices
outMatrix <- function(vec, index, filepath) {
	result = tryCatch(
    { #write out correlation data
    	N <- length(index)
    	#check if vector is correct length
    	if ((N * (N-1)/2) != length(vec)) {
    		stop('Vector length incorrect')
    	}
	    #create and fill strict upper triangle of matrix
	    mat <- matrix(, nrow = length(index), ncol = length(index))
	    mat[lower.tri(mat, diag=FALSE)] <- vec

	    #add column names
	    colnames(mat) <- index
	    #write out
	    options(digits = 3)	    
	    write.table(mat, file= filepath, quote = FALSE, sep = "\t", eol = "\n", na = " ", row.names = FALSE, fileEncoding = "")

	    #return success
	    return(1)

	    }, error = function(e) {
	    	print(e)
		#return -1 if error
		return(-1)
		}, finally = {
		#close connections
		closeAllConnections()
	}) #END tryCatch

	return(result)
}

sendUpdatedEdgeList <- function(edgeDF, name, isLast) {
	# o is an object array that can be passed to java
	    # array of columns (index1, index2, pVal)
	corrJavaObj <- .jnew("Correlate")
	o=.jarray(lapply(edgeDF, .jarray)) #lapply does it by column
	.jcall(corrJavaObj,"V","handleUpdatedEdgeList",o,name,isLast)
}

getPriorityGene <- function(array) {
	corrJavaObj <- .jnew("Correlate")
	priorityGene <- .jcall(corrJavaObj,"S","getPriorityGene")
	if (is.null(priorityGene)) return(-1)

	index = match(tolower(priorityGene), tolower(array))
	if (is.na(index)) return(-1)
	return(index)
}

#spearman correlation
spearman <- function(x, y) {
    cor(x, y, method = 'spearman')
}

getp2 <- function(cor, df) {
	pVec <- vector(mode = 'numeric', length = 1)
	.Call("pval", cor, df, pVec)
}






