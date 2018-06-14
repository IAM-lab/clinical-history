##################################################################################
# NAME:         DBSCAN_history.R
# AUTHOUR:      Alan Davies
# DATE:         02/12/2017
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
    package.args <- c("dbscan", "png")
    #"ggplot2", "binom", "tidyr", "reshape2", "plyr", "nlme", "lme4", "stargazer", "lmeresampler", "PMCMR") 
    #, "dplyr", "grDevices", "lattice", "plotrix", 
    
    loadPackages(package.args)
    
    # set wd
    setwd(path.expand("~"))
}

#---------------------------------------------------------------------------------
# FUNCTION:     openDataFile()
# INPUT:        String
# OUTPUT:       data.frame
# DESCRIPTION:  Returns a data frame of the requested file
#                
#---------------------------------------------------------------------------------
openDataFile <- function(file_path)
{
    data_file <- read.csv(paste0(file_path), header = TRUE, na.strings = c(" ", "NA", "-"))
    return(data_file)
}

#---------------------------------------------------------------------------------
# FUNCTION:     swapStimuliNames()
# INPUT:        void
# OUTPUT:       vector
# DESCRIPTION:  Swap stimuli label for more meaninful ECG name
#               using a substitution vector 
#---------------------------------------------------------------------------------
swapStimuliNames <- function()
{
    translation_vector <- c("AF.png" = "Atrial fibrillation", 
                            "anteriorSTEMI.png" = "Anterior STEMI",
                            "anterolateralSTEMI.png" = "Anterolateral STEMI",
                            "highlateralSTEMI.png" = "High lateral STEMI",
                            "inferiorSTEMIandAF.png" = "Inferior STEMI and AF",
                            "inferolateralSTEMI.png" = "Inferolateral STEMI",
                            "lateralSTEMI.png" = "Lateral STEMI",
                            "LBBB.png" = "LBBB",
                            "RBBB.png" = "RBBB")
    
    return(translation_vector)
}

#---------------------------------------------------------------------------------
# FUNCTION:     runDBSCAN()
# INPUT:        int, data.frame, String, list
# OUTPUT:       void
# DESCRIPTION:  Runs and displays results of DBSCAN algorithm, overlaying clusters
#               ontop of the stimuli images (paintings)
#---------------------------------------------------------------------------------
runDBSCAN <- function(condition, data, stimuli_name, args)
{
    trans_vec <- swapStimuliNames()
    title_str <- paste0("Condition ", condition, ": ", trans_vec[stimuli_name])
    
    # load background image 
    bkg_img <- readPNG(paste0(getwd(), "/Final-PHD-analysis/stimuli/", stimuli_name))
    
    # run DBSCAN and output results
    matrix_data <- as.matrix(na.omit(data))
    db <- dbscan(matrix_data, eps = args$eps, minPts = args$minPts)
    cat("\n", title_str, "\n")
    print(db)
    
    # output plots and add bkg image to them
    #print(pairs(matrix_data, col = db$cluster + 1L))
    #print(plot(matrix_data, col = db$cluster + 1L, main = title_str))  #res$cluster, main = title_str))
    #limits <- par()
    #rasterImage(bkg_img, limits$usr[1], limits$usr[3], limits$usr[2], limits$usr[4])
    #print(grid())
    
    # print raw fixation plot
    #print(plot(matrix_data, col = "blue", main = title_str))
    #rasterImage(bkg_img, limits$usr[1], limits$usr[3], limits$usr[2], limits$usr[4])
    #print(grid())
}

#---------------------------------------------------------------------------------
# FUNCTION:     runOPTICS()
# INPUT:        int, data.frame, String, list
# OUTPUT:       void
# DESCRIPTION:  Runs and displays results of OPTICS algorithm
#           
#---------------------------------------------------------------------------------
runOPTICS <- function(condition, data, stimuli_name, args)
{
    trans_vec <- swapStimuliNames()
    title_str <- paste0("Condition ", condition, ": ", trans_vec[stimuli_name])
    
    matrix_data <- as.matrix(na.omit(data))
    opt <- optics(matrix_data, eps = args$eps, minPts = args$minPts) #, xi = args$xi)
    
    cat("\n", title_str, "\n")
    print(opt)
    print(plot(opt))
}

#---------------------------------------------------------------------------------
# FUNCTION:     generatekNNDistPlots()
# INPUT:        data.frame, list
# OUTPUT:       void
# DESCRIPTION:  Produces kNN plot to determine optimal eps value by visualising
#               "knee" in plot curve
#---------------------------------------------------------------------------------
generatekNNDistPlots <- function(data, args)
{
    # plot kNN to determine optimal eps value
    matrix_data <- as.matrix(na.omit(data))
    print(kNNdistplot(matrix_data, k = args$minPts))
    #print(abline(h = 30, lty = 2))
}

#---------------------------------------------------------------------------------
# FUNCTION:     extractUniqueColData()
# INPUT:        data.frame, String
# OUTPUT:       vector
# DESCRIPTION:  Return a list of unique items from a given
#               data set by a given column name
#---------------------------------------------------------------------------------
extractUniqueColData <- function(data, col_label)
{
    return(as.vector(unlist(data[!duplicated(data[[col_label]]), 
                                 col_label], use.names = FALSE)))
}

#---------------------------------------------------------------------------------
# FUNCTION:     main()
# INPUT:        void
# OUTPUT:       void
# DESCRIPTION:  Main function. 
#               Makes all subsequent function calls
#---------------------------------------------------------------------------------
main <- function()
{
    initialize()
    #data <- loadDataFiles()
    
    data <- list()
    data[["condition1"]] <- openDataFile(paste0(getwd(), "/Final-PHD-analysis/data/cond1.csv"))
    data[["condition2"]] <- openDataFile(paste0(getwd(), "/Final-PHD-analysis/data/cond2.csv"))
    
    # rename the fixation columns
    colnames(data[["condition1"]])[12] <- "FixationPointX"
    colnames(data[["condition1"]])[13] <- "FixationPointY"
    colnames(data[["condition2"]])[12] <- "FixationPointX"
    colnames(data[["condition2"]])[13] <- "FixationPointY"
    
    stimuli <- extractUniqueColData(data[["condition1"]], "MediaName")
    cond_str <- "condition"
    stimulus_data <- list()
    algorithm_args <- list()
    k <- 1
    
    # set the properties of the algorithms 
    algorithm_args[["eps"]] <- 30               # 15
    algorithm_args[["minPts"]] <- 4             # 4
    algorithm_args[["xi"]] <- 0.05              # 0.05
    
    # loop over conditions 
    for(i in 1:length(data))
    {
        fixation_points_df <- NULL
        fixation_points <- NULL
       
        # get condition and subset
        cond_data <- data[[paste0(cond_str, i)]]
        
        # get unique participants and name list elements the same
        participants <- extractUniqueColData(cond_data, "ParticipantName")
        
        # loop over stimuli (paintings)
        for(j in 9:9) #length(stimuli))
        {
            # get the current stimuli
            current_stimuli <- cond_data[cond_data$MediaName == stimuli[j], ]
            current_stimuli <- current_stimuli[!duplicated(current_stimuli), ]
            
            stimulus_data[[k]] <- current_stimuli
            k <- k + 1
            
            for(l in 1:length(participants))
            {
                # aggregate fixation data accross all stimuli
                participant_FD <- current_stimuli[current_stimuli$ParticipantName == participants[l], ]
                if(length(participant_FD$GazeEventDuration) > 0)
                {
                    if(!is.na(participant_FD$FixationPointX) && !is.na(participant_FD$FixationPointY))
                    {
                        # get all x and y fixation points for each participant
                        x <- participant_FD$FixationPointX
                        y <- participant_FD$FixationPointY
                        
                        fixation_points <- cbind(x, y)
                    }
                }
                # put them together in a single data frame
                fixation_points_df <- rbind(fixation_points_df, fixation_points) 
            }

            # find optimal eps values
            generatekNNDistPlots(fixation_points_df, algorithm_args)
            
            # run DBSCAN and OPTICS algorithms
            runDBSCAN(i, fixation_points_df, stimuli[j], algorithm_args)
            runOPTICS(i, fixation_points_df, stimuli[j], algorithm_args)
        }
    }
}

# run main
main()