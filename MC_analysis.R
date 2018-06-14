##################################################################################
# NAME:         MC_analysis_1.R
# AUTHOUR:      Alan Davies
# DATE:         31/10/2017
# INSTITUTION:  Interaction Analysis and Modelling Lab (IAM), University of Manchester
# DESCRIPTION:  Markov chain analysis of bi-gram transition for the second ECG 
#               experiement. Analsis using top-down AOIs generated in Tobii Studio.
#               Uses the Jensen-Shannon distance and the Hellinger distance
#               metrics that can be chosen (input at runtime).
#               Analysis for HPC vs no HPC for those seeing either first.
##################################################################################

#---------------------------------------------------------------------------------
# FUNCTION:     openDataFiles()
# INPUT:        String, vector
# OUTPUT:       list
# DESCRIPTION:  Returns a list of data files accessable by label
#                
#---------------------------------------------------------------------------------
openDataFiles <- function(file_path, conditions)
{
    data_files <- list()
 
    for(i in 1:length(conditions))
    {
        # load the file
        data_file <- read.csv(paste0(file_path, conditions[i]), 
                              header = TRUE, na.strings = c(" ", "NA", "-"))
        
        # store it in list by label
        data_files[[paste0("condition_", i)]] <- data_file
    }

    return(data_files)
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
# FUNCTION:     getPresentationSequence()
# INPUT:        data.frame, String
# OUTPUT:       vector
# DESCRIPTION:  Returns a vector of participant ID's for a specified presentation
#               sequence
#---------------------------------------------------------------------------------
getPresentationSequence <- function(data, presentation_seq)
{
    data <- data[!duplicated(data$Participant), ]
    return(as.vector(data[data$PresentationSequence == presentation_seq, "Participant"]))
}

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
            #cat("\nPackage <", i, "> not found, attempting to add it...\n\n")
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
    # required packages list
    package.args <- c("ggplot2", "crayon", "dplyr", "plyr", "tidyr") #, "dplyr", "grDevices", "lattice", "plotrix", "plyr", 
                      #"tidyr", "lme4", "stargazer", "lmeresampler") 
    
    # load packages and set working dir
    loadPackages(package.args)
    setwd(path.expand("~"))
}

#---------------------------------------------------------------------------------
# FUNCTION:     convertToMatrix()
# INPUT:        data.frame, BOOL
# OUTPUT:       matrix
# DESCRIPTION:  Returns matrix representing probabilty of AOI transitions
#               
#---------------------------------------------------------------------------------
convertToMatrix <- function(data, diagonal = TRUE)
{
    tmp0 <- data
    tmp <- tmp0 %>% group_by(participant) %>% mutate(to = lead(AOI))
    tmp2 <- tmp[complete.cases(tmp), ]
    with(tmp2, table(AOI, to))
    out_mat <- as.matrix(with(tmp2, table(AOI, to)))
    
    # can remove diagonal
    if(!diagonal) diag(out_mat) <- 0
    
    # produce probability matrix
    prob_mat <- out_mat / rowSums(out_mat)
    return(prob_mat)
}

#---------------------------------------------------------------------------------
# FUNCTION:     convertToPrior()
# INPUT:        matrix
# OUTPUT:       matrix
# DESCRIPTION:  Bayesian prior (Dirichlet) method
#               
#---------------------------------------------------------------------------------
convertToPrior <- function(data)
{
    m <- 0
    matrix_len <- nrow(data) 
    for(i in 1:matrix_len)
    {
        m <- sum(data[i, ] * 10)
        for(j in 1:matrix_len)
            data[i, j] <- ((data[i, j] * 10) + 1) / (m + matrix_len)
    }
    return(data)
}

#---------------------------------------------------------------------------------
# FUNCTION:     calculateJSDistance()
# INPUT:        matrix, matrix
# OUTPUT:       double
# DESCRIPTION:  Computes correct Jensen-Shannon distance between two matrices
#               
#---------------------------------------------------------------------------------
calculateJSDistance <- function(m1, m2)
{
    results <- list()
    JS_distance <- 0
    summed_JS_distance <- 0
    matrix_length <- nrow(m1)
    coeff1 <- 0
    coeff2 <- 0
    
    for(i in 1:matrix_length)
    {
        for(j in 1:matrix_length)
        {
            # calculate KLD per row
            coeff1 <- coeff1 + (m1[i, j] * log(m1[i, j] / (0.5 * (m1[i, j] + m2[i, j]))))
            coeff2 <- coeff2 + (m2[i, j] * log(m2[i, j] / (0.5 * (m1[i, j] + m2[i, j])))) 
        }
        JS_distance <- sqrt(0.5 * (coeff1 + coeff2))
        summed_JS_distance <- summed_JS_distance + JS_distance
        JS_distance <- 0
        coeff1 <- 0
        coeff2 <- 0
    }
    average_distance <- summed_JS_distance / matrix_length
    return(average_distance)
}

#---------------------------------------------------------------------------------
# FUNCTION:     calculateHellingerDistance()
# INPUT:        matrix, matrix
# OUTPUT:       Summed Hellinger Distance
# DESCRIPTION:  Defined as: 1/sqrt(2) * sqrt(sum(square(sqrt(pi - squrt(qi))))
#---------------------------------------------------------------------------------
calculateHellingerDistance <- function(m1, m2)
{
    length_of_matrix <- nrow(m1)
    HPQ <- 0
    HJ <- 0
    
    for(i in 1:length_of_matrix)
    {
        for(j in 1:length_of_matrix)
        {
            HJ <- HJ + ((sqrt(m1[i, j]) - sqrt(m2[i, j])) ^ 2)
        }
        HPQ <- HPQ + ((1 / sqrt(2)) * sqrt(HJ))
        HJ <- 0
    }
    return(HPQ / length_of_matrix) 
}

#---------------------------------------------------------------------------------
# FUNCTION:     calculatePvalue()
# INPUT:        vector, vector
# OUTPUT:       void
# DESCRIPTION:  Calculate p-value (% of values > correct/incorrect value)
#               
#---------------------------------------------------------------------------------
calculatePvalue <- function(correct_and_incorrect, shuffled_distances)
{
    gtr <- length(shuffled_distances[shuffled_distances > correct_and_incorrect])
    pvalue <- gtr / length(shuffled_distances) 
    return(pvalue)
}

#---------------------------------------------------------------------------------
# FUNCTION:     generateDensityPlot()
# INPUT:        double, vector, list, list
# OUTPUT:       void
# DESCRIPTION:  Output density plot (distribution plot) for the results of the
#               permutation test. Purple line indicates primary group distance
#---------------------------------------------------------------------------------
generateDensityPlot <- function(distance_result, shuffled_distances, stimuli_data, settings)
{
    metric <- ifelse(settings$distance_metric == 1, "Jensen-Shannon Distance", "Hellinger Distance")
    density_plot <- density(shuffled_distances)
    plot(density_plot, type = "n", main = stimuli_data$label, xlab = metric, panel.first = grid())
    polygon(density_plot, col = "lightgray", border = "grey")
    rug(shuffled_distances, col = ifelse(shuffled_distances == distance_result, 'blue', 'red'))
    print(abline(v = distance_result, col = "purple"))
    print(density_plot)
}

#---------------------------------------------------------------------------------
# FUNCTION:     generateReport()
# INPUT:        list, list, list
# OUTPUT:       void
# DESCRIPTION:  Output results of distance and p-value 
#               
#---------------------------------------------------------------------------------
generateReport <- function(settings, stimuli_data, report_args)
{
    cat("\n\nData:", stimuli_data$label, "\n")
    
    if(settings$distance_metric == "1")
    {
        df <- data.frame(Jsd = report_args$distance, p.value = report_args$pvalue, Permutations = settings$num_perms)
    }
    else
    {
        df <- data.frame(Hd = report_args$distance, p.value = report_args$pvalue, Permutations = settings$num_perms)
    }
    print(df) 
    cat("\nGroup 1 (n):", report_args$groupsize1)
    cat("\nGroup 2 (n):", report_args$groupsize2, "\n\n")  
}

#---------------------------------------------------------------------------------
# FUNCTION:     generateSettingsReport()
# INPUT:        list
# OUTPUT:       void
# DESCRIPTION:  Output relevant settings used for the study 
#               
#---------------------------------------------------------------------------------
generateSettingsReport <- function(settings)
{
    cat("\n\nStudy settings:\n\n")
    cat("Number of permutations:", settings$num_perms)
    cat("\nDistance:", ifelse(settings$distance_metric == "1", "Jensen_Shannon", "Hellinger"))
    if(settings$distance_metric == "1") cat("\nPrior:", toupper(settings$prior))
    cat("\nDiagonal transitions removed:", toupper(settings$diagonal), "\n\n")
}

#---------------------------------------------------------------------------------
# FUNCTION:     promptForSettings()
# INPUT:        void
# OUTPUT:       list
# DESCRIPTION:  Prompt the user to enter settings for the various test
#               parameters
#---------------------------------------------------------------------------------
promptForSettings <- function()
{
    settings <- list()
    cat(cyan(rep("-", 40)))
    cat(cyan("\nPlease enter test perameters (HPC = History of Presenting Complaint)\n"))
    cat(cyan(rep("-", 40)), "\n")
    settings[["num_perms"]] <- readline("Enter number of permutations: ")  
    settings[["diagonal"]] <- readline("Remove diagonal transitions (y/n): ")  
    settings[["distance_metric"]] <- readline("Enter 1 for Jensen-Shannon or 2 for Hellinger distance: ")  
    settings[["echo_locations"]] <- readline("Echo locations [show AOI sequences] (y/n): ")
    settings[["trans_plots"]] <- readline("Output transition matricies (y/n): ")
    settings[["density_plots"]] <- readline("Output distribution plots (y/n): ")
    if(settings[["distance_metric"]] == "1")
        settings[["prior"]] <- readline("Apply Bayesian prior [dirichlet method] (y/n): ")  
    
    return(settings)
}

#---------------------------------------------------------------------------------
# FUNCTION:     getUniqueData()
# INPUT:        data.frame, String
# OUTPUT:       vector
# DESCRIPTION:  Returns a vector of unique values from a specified data frame
#               i.e. all participants or stimuli in a data set
#---------------------------------------------------------------------------------
getUniqueData <- function(data, parameter)
{
    return(as.vector(unlist(data[!duplicated(data[[parameter]]), parameter], use.names = FALSE)))
}

#---------------------------------------------------------------------------------
# FUNCTION:     getAOIHits()
# INPUT:        data.frame, data.frame, list
# OUTPUT:       data.frame
# DESCRIPTION:  Get the locations for each fixation within the stimulus 
#               
#---------------------------------------------------------------------------------
getAOIHits <- function(data, AOI_locations, condition, settings) 
{
    AOI <- NULL
    AOI_locations <- AOI_locations[AOI_locations$condition == condition, ]
    
    echo <- ifelse(settings$echo_locations == "y", TRUE, FALSE) 
    
    stimuli <- getUniqueData(AOI_locations, "ECG")
    from_col <- which(colnames(data) == "ParticipantName")
    to_col <- which(colnames(data) == "SaccadicAmplitude")
    
    # extract data from 1 column to another & rename cols
    main_data <- data[c(from_col:to_col)]
    colnames(main_data)[4] <- "MediaPosX"
    colnames(main_data)[5] <- "MediaPosY"
    colnames(main_data)[12] <- "PosX"
    colnames(main_data)[13] <- "PosY"
    
    if(echo) cat(blue("\n Computing locations:\n\n"))
    
    for(i in 1:length(stimuli))
    {
        # get the data for the AOI and stimuli
        if(echo) cat(blue("\n\nStimuli: ", stimuli[i]), "\n\n")
        AOI_stimuli <- AOI_locations[AOI_locations$ECG == stimuli[i], ]
        current_stimuli <- main_data[main_data$MediaName == stimuli[i], ]
        
        for(j in 1:nrow(current_stimuli))
        {
            # detect fixations within AOIs
            current_row <- current_stimuli[j, ]
            AOI <- c(AOI, detectHit(AOI_stimuli, current_row, echo))
        }
    }

    # append the AOI column to data and remove NA fixations from data
    main_data <- cbind(main_data, AOI)
    main_data <- main_data[complete.cases(main_data[, "AOI"]), ]
    
    return(main_data)
}

#---------------------------------------------------------------------------------
# FUNCTION:     detectHit()
# INPUT:        data.frame, data.frame
# OUTPUT:       String
# DESCRIPTION:  Return the name of the AOI that the fixation is detected 
#               within
#---------------------------------------------------------------------------------
detectHit <- function(AOI_data, current_row, echo)
{
    AOI <- NA
    
    # important need to add offset to fixations to compute correct position
    fixation_x <- (current_row$PosX + current_row$MediaPosX)
    fixation_y <- (current_row$PosY + current_row$MediaPosY)

    # if any are NA return NA
    if((is.null(fixation_x) || is.null(fixation_y)) || (is.na(fixation_x) || is.na(fixation_y)))
        return(NA)
    
    # if fixations are outside the stimuli boundaries then return NA
    if((fixation_x < current_row$MediaPosX) || (fixation_x > (current_row$MediaPosX + current_row$MediaWidth)) ||
       (fixation_y < current_row$MediaPosY) || (fixation_y > (current_row$MediaPosY + current_row$MediaHeight)))
    {
        return(NA)
    }
    
    for(i in 1:nrow(AOI_data))
    {
        AOI_row <- AOI_data[i, ]
      
        # add the offsets
        x_offset <- current_row$MediaPosX + AOI_row$left
        y_offset <- current_row$MediaPosY + AOI_row$top
        
        # detect hit inside AOI
        if((fixation_x > x_offset) && (fixation_x < (x_offset + AOI_row$width)) &&
           (fixation_y > y_offset) && (fixation_y < (y_offset + AOI_row$height))) 
        {
            AOI <- toString(AOI_row$AOI)
            if(echo) cat(AOI)
            break
        }
    }
    return(AOI)
}

#---------------------------------------------------------------------------------
# FUNCTION:     randomPermutations()
# INPUT:        data.frame, list, int, int
# OUTPUT:       vector
# DESCRIPTION:  Generate 2 random sub groups the same size as the original groups
#               and compare distances for the selected amout of permutations.
#---------------------------------------------------------------------------------
randomPermutations <- function(data, settings, group1_size, stimuli)
{
    progress <- 0
    distance <- NULL
    perms <- as.numeric(settings$num_perms)
    
    # create progress bar
    pb_title <- paste0("Computing Distances ", stimuli$label)
    progress_bar <- winProgressBar(title = pb_title, min = 0, max = perms, width = 300)
    getWinProgressBar(progress_bar)
    
    # loop over the permutations  
    for(i in 1:perms) 
    {
        participants <- as.data.frame(data[!duplicated(data$participant), "participant"])
        colnames(participants) <- "participant"
        
        # get a random subset of participants and store in group 1 whatever is left put in second group
        group1_participants <- participants[sample(unique(nrow(participants)), group1_size), ]
        group2_participants <- as.vector(participants[!(participants$participant %in% group1_participants), ])   
        
        # get the selected groups data
        group1 <- data[which(data$participant %in% group1_participants), ]
        group2 <- data[which(data$participant %in% group2_participants), ]
        
        # build Markov chains
        mat_diag <- ifelse(settings$diagonal == "y", FALSE, TRUE)
        chain_1 <- convertToMatrix(group1, mat_diag)
        chain_2 <- convertToMatrix(group2, mat_diag)
        
        # remove any NaN
        chain_1[is.na(chain_1)] <- 0
        chain_2[is.na(chain_2)] <- 0
        
        # get distance 
        if(settings$distance_metric == "1")
        {
            if(settings$prior == "y")
            {
                # convert to prior
                chain_1 <- convertToPrior(chain_1)
                chain_2 <- convertToPrior(chain_2)                       
            }
            
            # Jensen-Shannon
            distance_results <- calculateJSDistance(chain_1, chain_2)
        }
        else
        {
            # Hellinger
            distance_results <- calculateHellingerDistance(chain_1, chain_2)
        }
        
        # store computed distance in vector and return
        distance <- c(distance, distance_results)
        
        # display progress bar
        setWinProgressBar(progress_bar, i, title = paste(round(i / perms * 100, 0), "% processed [Computing Distance]"))
    }
    
    # close the progress bar & return result
    close(progress_bar)
    return(distance)
}

#---------------------------------------------------------------------------------
# FUNCTION:     transitionPlots()
# INPUT:        data.frame, data.frame, list
# OUTPUT:       void
# DESCRIPTION:  Output transition plots for the 2 selected groups 
#               normalised by max value 
#---------------------------------------------------------------------------------
transitionPlots <- function(c1, c2, stimuli_data)
{
    lead_names <- c("I", "II", "III", "aVR", "aVL", "aVF", "V1", "V2", "V3", "V4", "V5", "V6", "RS1", "RS2", "RS3")
    c1 <- c1[c1$MediaName == stimuli_data$file_name, c("participant", "AOI")]
    c2 <- c2[c2$MediaName == stimuli_data$file_name, c("participant", "AOI")]
    
    # convert to frequencies 
    c1 <- statetable.msm(c1$AOI, subject = 1)
    c2 <- statetable.msm(c2$AOI, subject = 1)
    
    # normalise by max value
    max_val <- max(max(c1), max(c2))
    
    for(i in 1:2)
    {
        plot_title <- stimuli_data$label
        if(i == 1){
            data <- c1
            plot_title <- paste0(plot_title, " \nwith HPC")
        } else {
            data <- c2
            plot_title <- paste0(plot_title, " \nwithout HPC")
        }
       
        # output transition matrix
        tp <- levelplot(data, col.regions = colorpanel(max_val, "white", "grey10"), 
                        at = unique(c(seq(from = 0, to = max_val, by = 5))),
                        main = plot_title, ylab = "AOI (destination)", xlab="AOI (Origin)",
                        scales = list(x = list(labels = lead_names, cex = .5), y = list(labels = lead_names, cex = .5)))
        
        print(tp)
    }
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
    report_args <- list()
    primary_conditions <- c("cond1.csv", "cond2.csv")
    
    initialize()
    
    if(interactive()){
        settings <- promptForSettings()
    }
    
    # load additional source files with stimuli metadata
    common_source_data <- paste0(getwd(), "/Final-PHD-analysis/experiment_history.R")
    source(common_source_data) 
    
    # open the data files
    ECG_data_files <- openDataFiles(paste0(getwd(), "/Final-PHD-analysis/data/"), primary_conditions)
    
    # open the answers and the AOI locations
    given_answers <- openDataFile(paste0(getwd(), "/Final-PHD-analysis/data/gs/given_answers.csv"))
    AOI_locations <- openDataFile(paste0(getwd(), "/Final-PHD-analysis/data/AOI.csv"))
    
    # get the people that saw HPC first and those that say no HPC first
    saw_HPC <- getPresentationSequence(given_answers, "1-2")
    saw_noHPC <- getPresentationSequence(given_answers, "2-1")
    
    # get the data for the groups
    group_HPC <- ECG_data_files[[1]][which(ECG_data_files[[1]]$ParticipantName %in% saw_HPC), ]
    group_noHPC <- ECG_data_files[[2]][which(ECG_data_files[[2]]$ParticipantName %in% saw_noHPC), ]
    
    group_HPC <- getAOIHits(group_HPC, AOI_locations, 1, settings)
    group_1 <- group_HPC[ ,c("ParticipantName", "MediaName", "AOI")]
    colnames(group_1)[1] <- "participant"
    group_noHPC <- getAOIHits(group_noHPC, AOI_locations, 2, settings)
    group_2 <- group_noHPC[ ,c("ParticipantName", "MediaName", "AOI")]
    colnames(group_2)[1] <- "participant"
    
    for(i in 1:length(getStimuliList()))
    {
        # get data for the current stimuli
        current_stimuli_metadata <- getExperimentSetupData(getStimuliList()[i])
        
        # extract stimulus data
        chain_1 <- group_1[group_1$MediaName == current_stimuli_metadata$file_name, c("participant", "AOI")]
        chain_2 <- group_2[group_2$MediaName == current_stimuli_metadata$file_name, c("participant", "AOI")]
        combined_data <- rbind(chain_1, chain_2)
        
        mat_diag <- ifelse(settings$diagonal == "y", FALSE, TRUE) 
     
        # convert to matrix form
        chain_1 <- convertToMatrix(chain_1, mat_diag)
        chain_2 <- convertToMatrix(chain_2, mat_diag)
        
        # produce transition plots
        if(settings$trans_plots == "y")
        {
            transitionPlots(group_1, group_2, current_stimuli_metadata)
        }
        
        # remove any NA's from the chains
        chain_1[is.na(chain_1)] <- 0
        chain_2[is.na(chain_2)] <- 0
        
        if(settings$distance_metric == "1") 
        {
            if(settings$prior == "y")
            {
                # if using Jsd and a prior then convert to prior
                chain_1 <- convertToPrior(chain_1)
                chain_2 <- convertToPrior(chain_2)    
            }
            
            # calculate distance (Jsd)
            distance_result <- calculateJSDistance(chain_1, chain_2)
        }
        if(settings$distance_metric == "2")
        {
            # calculate distance (Hd)
            distance_result <- calculateHellingerDistance(chain_1, chain_2)
        }
        
        # run distance on random groups for n permutations
        shuffled_distances <- randomPermutations(combined_data, settings, length(saw_HPC), current_stimuli_metadata)
        
        if(settings$density_plots == "y")
        {
            # output a distribution plot
            generateDensityPlot(distance_result, shuffled_distances, current_stimuli_metadata, settings)
        }
        
        # build report arguments
        report_args$distance <- distance_result
        report_args$pvalue <- calculatePvalue(distance_result, shuffled_distances) 
        report_args$groupsize1 <- length(saw_HPC)
        report_args$groupsize2 <- length(saw_noHPC)
        
        # generate a textual report
        generateReport(settings, current_stimuli_metadata, report_args)
    }  
    
    # make a final report on the settings used for the study
    generateSettingsReport(settings)
}

# run main
main()
