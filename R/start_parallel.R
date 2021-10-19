#
#============================ Start parallel computing for variable selection
# Borrowed from 'clustvarsel', author: Luca Scrucca


start_parallel <- function(parallel = TRUE, ...)
    # Start parallel computing
{
    # set default parallel functionality depending on system OS:
    parallelType <- if(.Platform$OS.type == "windows") 
        "snow" else "multicore"
    
    # get the current number of cores available
    numCores <- parallel::detectCores()
    
    # set parameters for parallelization
    if ( is.logical(parallel) ) { 
        NULL 
    } else if ( is.numeric(parallel) ) {
        numCores <- as.integer(parallel)
        parallel <- TRUE  
    } else parallel <- FALSE
    
    attr(parallel, "type") <- parallelType
    attr(parallel, "cores") <- numCores
    
    # start "parallel backend" if needed
    if ( parallel )
    { 
        if ( parallelType == "snow" )
        { 
            # snow functionality on Unix-like systems & Windows
            cl <- parallel::makeCluster(numCores, type = "PSOCK")
            attr(parallel, "cluster") <- cl
            doParallel::registerDoParallel(cl, cores = numCores)
            #
        } else if(parallelType == "multicore") {
            # multicore functionality on Unix-like systems
            cl <- parallel::makeCluster(numCores, type = "FORK")
            doParallel::registerDoParallel(cl, cores = numCores) 
            attr(parallel, "cluster") <- cl   
        }
    }
    
    return(parallel)
}
