
#' Run MCMC iterations
#' 
#' @param track Dataframe of data, with columns "x", "y", "time", and "ID"
#' @param nbStates Number of states
#' @param nbIter Number of iterations
#' @param inits List of initial parameters 
#' (beta, sigma, Q, state)
#' @param priors List of parameters of prior distributions, with components:
#' \itemize{
#'   \item{"mean":} Vector of means for normal priors on movement parameters, of length 2*nbStates
#'   \item{"sd":} Vector of standard deviations for normal priors on movement parameters, of length
#' 2*nbStates
#'   \item{"shape":} Vector of shapes of gamma priors for the transition rates
#'   \item{"rate":} Vector of rates of gamma priors for the transition rates
#'   \item{"con":} Vector of concentrations of Dirichlet priors for transition probabilities
#' }
#' @param props List of parameters of proposal distributions, with components:
#' \itemize{
#'   \item{"betaSD":} Scalar standard deviation for normal proposal distribution of beta
#'   \item{"sigmaSD":} Scalar standard deviation for normal proposal distribution of sigma
#'   \item{"updateLim":} Vector of two values: min and max length of updated state sequence
#'   \item{"updateProbs":} Probability for each element of updateLim[1]:updateLim[2] (if NULL,
#'   all values are equiprobable)
#' }
#' @param tunes List of tuning parameters, with components:
#' \itemize{
#'   \item{"thinStates":} Thinning factor for the posterior state sequences (needed because
#' of memory limitations)
#' }
#' @param kalmanpars List of parameters of the Kalman filter, with components:
#' \itemize{
#'   \item{"Hmat":} Matrix of observation error variance (four columns, and one row 
#' for each row of data)
#'   \item{"a0":} List of initial state estimate vectors, named for each individual
#'   \item{"P0":} List of initial estimate covariance matrices, named for each individual
#' }
#' @param updateState Logical. If FALSE, the state process is not updated
#' (for exploratory analysis only)
#' 
#' @importFrom MASS mvrnorm
#' @importFrom stats dnorm runif rnorm rexp rgamma
#' @export
#' 
#' @useDynLib MScrawl
runMCMC <- function(track, nbStates, nbIter, inits, priors, props, tunes, kalmanpars,
                    updateState=TRUE)
{
    
    ######################
    ## Unpack arguments ##
    ######################
    # initial parameters
    beta <- inits$beta
    sigma <- inits$sigma
    param <- c(beta,sigma)
    Q <- rep(list(inits$Q), length( unique( track$ID ) ) )
    names(Q) <- unique( track$ID ) 
    state0 <- inits$state
    if(is.null(beta) | is.null(sigma) | is.null(Q) | is.null(state0))
        stop("'inits' should have components beta, sigma, Q, and state")
    
    # Kalman filter parameters
    Hmat <- kalmanpars$Hmat
    a0 <- kalmanpars$a0
    P0 <- kalmanpars$P0
    if(is.null(Hmat) | is.null(a0) | is.null(P0))
        stop("'kalmanpars' should have components Hmat, a0, and P0")
    
    # unpack prior parameters
    priorMean <- priors$mean
    priorSD <- priors$sd
    priorShape <- priors$shape
    priorRate <- priors$rate
    priorCon <- priors$con
    if(is.null(priorMean) | is.null(priorSD) | is.null(priorShape) | 
       is.null(priorRate) | is.null(priorCon)) {
        stop("'priors' should have components mean, sd, shape, rate, and con")
    }
    
    # unpack proposal parameters
    betaPropSD <- props$betaSD
    sigmaPropSD <- props$sigmaSD
    updateLim <- props$updateLim
    updateProbs <- props$updateProbs
    if(is.null(betaPropSD) | is.null(sigmaPropSD) | is.null(updateLim))
        stop("'props' should have components betaSD, sigmaSD, and updateLim")
    
    if(is.null(updateProbs))
        updateProbs <- rep(1, length(updateLim[1]:updateLim[2]))/length(updateLim[1]:updateLim[2])
    
    if(length(updateLim[1]:updateLim[2])!=length(updateProbs))
        stop("'updateProbs' has the wrong length")
    
    # unpack tuning parameters
    thinStates <- tunes$thinStates
    
    ####################
    ## Initialisation ##
    ####################
    # Prepare data structures
    if( !( "ID" %in% colnames( track ) ) ) {
      warning( "'track' should have column 'ID', assuming data is one individual"  )
      track$ID <- 1
    }
    
    ids <- as.character(unique(track$ID))
    obs <- list()
    switch <- list()
    data <- list()
    HmatAll <- list()
    
    for( i in ids ) {
      nbObs <- nrow(track[ which( track$ID == i ), ])
      obs[[ i ]] <- matrix( c( track[ which( track$ID == i ), "x" ],
                               track[ which( track$ID == i ), "y" ], 
                               track[ which( track$ID == i ), "time" ], 
                               state0[ which( track$ID == i )  ] ), 
                            ncol = 4 )
      colnames( obs[[ i ]] ) <- c( "x", "y", "time", "state" )
      indSwitch <- which( obs[[ i ]][ -1, "state" ] != obs[[ i ]][ -nbObs, "state" ] ) + 1
      switch[[ i ]] <- matrix( c( obs[[ i ]][ indSwitch, "time" ] - 0.001, rle( obs[[ i ]][ , "state" ] )$values[ -1 ] ), ncol = 2 )
      colnames( switch[[ i ]] ) <- c( "time", "state" )
      data[[ i ]] <- rbind( obs[[i]] ,cbind( NA, NA, switch[[ i ]] ) )
      data[[ i ]] <- data[[ i ]][ order( data[[ i ]][ , "time" ] ), ]
      
      # initialise Hmat (rows of 0s for transitions)
      HmatAll[[ i ]] <- matrix( 0, nrow( data[[ i ]] ), 4 )
      HmatAll[[ i ]][ which( !is.na( data[[ i ]][ , "x" ] ) ), ] <- Hmat[ which( track$ID == i ), ]
    }
    
    # initial likelihood
    oldllk <- lapply( ids, function( id ) { kalman_rcpp( data = data[[ id ]], param = param, Hmat = HmatAll[[ id ]], a0 = a0[[ id ]], P0 = P0[[ id ]] ) } )
    oldllk <- do.call( 'sum', oldllk )
    # initial log-prior
    oldlogprior <- sum( dnorm( log( param ), priorMean, priorSD, log = TRUE ) )
    
    ###############################
    ## Loop over MCMC iterations ##
    ###############################
    allparam <- matrix(NA,nrow=nbIter,ncol=2*nbStates)
    allrates <- array( NA, dim = c( nbIter, nbStates*(nbStates-1), length( ids ) ) )
    allstates <- matrix(NA,nrow=nbIter/thinStates,ncol=nbObs) # uses a lot of memory!
    accSwitch <- rep(0,nbIter)
    accParam <- rep(0,nbIter)
    allLen <- rep(NA,nbIter)
    
    t0 <- Sys.time()
    for(iter in 1:nbIter) {
        if(iter%%100==0)
            cat("\rIteration ", iter, "/", nbIter, "... ", round(Sys.time()-t0,2),
                " -- accSwitch = ", round(sum(accSwitch)/iter*100), "%",
                " -- accParam = ", round(sum(accParam)/iter*100), "%", sep="")
        
        ######################################
        ## 1. Update discrete state process ##
        ######################################
        if(updateState) {
          
            # pick an individual at random for which we will update the state sequence
            id <- sample( ids, 1 ) 
            upState <- updateState( obs = obs[[ id ]], switch = switch[[ id ]], updateLim = updateLim, 
                                   updateProbs = updateProbs, Q = Q[[ id ]])
            newData <- data
            newData[[ id ]] <- upState$newData
            newSwitch <- switch
            newSwitch[[ id ]] <- upState$newSwitch
            allLen[iter] <- upState$len
            
            # update Hmat (rows of 0s for transitions)
            newHmatAll <- HmatAll
            newHmatAll[[ id ]] <- matrix( 0, nrow( newData[[ id ]] ), 4 )
            newHmatAll[[ id ]][ which( !is.na( newData[[ id ]] [ , "x" ] ) ), ] <- Hmat[ which( track$ID == id ), ] 
            
            # Calculate acceptance ratio
            newllk <- lapply( ids, function( id ) { kalman_rcpp( data = newData[[ id ]], param = param, Hmat = newHmatAll[[ id ]], a0 = a0[[ id ]], P0 = P0[[ id ]] ) } )
            newllk <- do.call( 'sum', newllk )
            logHR <- newllk - oldllk
            
            if(log(runif(1))<logHR) {
                # Accept new state sequence
                accSwitch[iter] <- 1
                switch <- newSwitch
                data <- newData
                obs <- lapply(data,function(data) {data[!is.na(data[,"x"]),]})
                oldllk <- newllk
                HmatAll <- newHmatAll
            }
        }
        
        ###################################
        ## 2. Update movement parameters ##
        ###################################
        # On working scale
        betaprimeW <- rnorm(nbStates,log(beta),betaPropSD)
        sigmaprimeW <- rnorm(nbStates,log(sigma),sigmaPropSD)
        newlogprior <- sum(dnorm(c(betaprimeW,sigmaprimeW),priorMean,priorSD,log=TRUE))
        
        # On natural scale
        betaprime <- exp(betaprimeW)
        sigmaprime <- exp(sigmaprimeW)
        
        # Calculate acceptance ratio
        newllk <- lapply( ids, function( id ) { kalman_rcpp( data = data[[ id ]], param = c( betaprime, sigmaprime ), Hmat = HmatAll[[ id ]], a0 = a0[[ id ]], P0 = P0[[ id ]] ) } )
        newllk <- do.call( 'sum', newllk )
        logHR <- newllk + newlogprior - oldllk - oldlogprior

        if(log(runif(1))<logHR) {
            # Accept new parameter values
            accParam[iter] <- 1
            beta <- betaprime
            sigma <- sigmaprime
            param <- c(beta,sigma)
            oldllk <- newllk
            oldlogprior <- newlogprior
        }
        
        ###############################
        ## 3. Update switching rates ##
        ###############################
        Q <- lapply( ids, function( id ) { updateQ( nbStates = nbStates, data = data[[ id ]], switch = switch[[ id ]], 
                     priorShape = priorShape, priorRate = priorRate,
                     priorCon = priorCon ) } )
        
        #########################
        ## Save posterior draw ##
        #########################
        allparam[iter,] <- param
        allrates[ iter, , ] <- matrix( unlist( lapply( Q, function( q ){ q[ !diag( nbStates ) ] } ) ), ncol = length( ids ), nrow = nbStates*(nbStates-1) ) 
        if(iter%%thinStates==0)
            allstates[iter/thinStates,] <- obs[,"state"]
    }
    cat("\n")
    
    return(list(allparam = allparam,
                allrates = allrates,
                allstates = allstates,
                accSwitch = accSwitch,
                allLen = allLen))
}
