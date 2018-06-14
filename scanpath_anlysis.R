##################################################################################
# NAME:         .R
# AUTHOUR:      Alan Davies
# DATE:         31/10/2017
# INSTITUTION:  Interaction Analysis and Modelling Lab (IAM), University of Manchester
# DESCRIPTION:  
#               
##################################################################################

#---------------------------------------------------------------------------------
# FUNCTION:     loadPackages(package.args)
# INPUT:        vector
# OUTPUT:       void
# DESCRIPTION:  Loads required packages.
#                
#---------------------------------------------------------------------------------
loadPackages <- function(package.args)
{ 
    for(i in package.args)
    {
        if(!is.element(i, .packages(all.available = TRUE)))
        {
            cat("\nPackage <", i, "> not found, attempting to add it...")
            install.packages(i)
        }
        library(i, character.only = TRUE)
    }
}

#---------------------------------------------------------------------------------
# FUNCTION:     initialize()
# INPUT:        void
# OUTPUT:       void
# DESCRIPTION:  Set up function for adding packages and other source data
#               
#---------------------------------------------------------------------------------
initialize <- function()
{
    # load packages
    package.args <- c("ggplot2") #, "dplyr", "grDevices", "lattice", "plotrix", "plyr", 
                      #"tidyr", "lme4", "stargazer", "lmeresampler") 
    
    loadPackages(package.args)
    
    # set wd
    setwd(path.expand("~"))
}

#---------------------------------------------------------------------------------
# FUNCTION:     main()
# INPUT:        void
# OUTPUT:       void
# DESCRIPTION:  Main function. 
#               Makes all subsequent function calls.     
#---------------------------------------------------------------------------------
main <- function()
{
    initialize()
    
    # load additional source files with stimuli metadata
    common_source_data <- paste0(getwd(), "/experiment_history.R")
    source(common_source_data) 
    
    
}

# run main
main()
