##################################################################################
# NAME:         rhythm_stip_analysis.R
# AUTHOUR:      Alan Davies
# DATE:         04/12/2017
# INSTITUTION:  Interaction Analysis and Modelling Lab (IAM), University of Manchester
# DESCRIPTION:  Analysis of the rhythm strip portion of the experiment
#               
#               
#               
#               
#               
##################################################################################

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
    package.args <- c("ggplot2", "crayon", "dplyr", "tidyr", "grDevices", "lattice", "msm", "stats")
    
    # load packages and set working dir
    loadPackages(package.args)
    setwd(path.expand("~"))
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
getAOIHits <- function(data, AOI_locations, echo = FALSE) 
{
    AOI <- NULL
    
    stimuli <- getUniqueData(data, "MediaName")
    from_col <- which(colnames(data) == "ParticipantName")
    to_col <- which(colnames(data) == "SaccadicAmplitude")
    
    # extract data from 1 column to another & rename cols
    main_data <- data[c(from_col:to_col)]
    colnames(main_data)[4] <- "MediaPosX"
    colnames(main_data)[5] <- "MediaPosY"
    colnames(main_data)[16] <- "PosX"
    colnames(main_data)[17] <- "PosY"
    
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
            if(echo) cat(AOI,"")
            break
        }
    }
    return(AOI)
}

#---------------------------------------------------------------------------------
# FUNCTION:     processWaveformQuestion()
# INPUT:        data.frame
# OUTPUT:       void
# DESCRIPTION:  Produce a histogram for participant responses to question about
#               which part of the waveform they paid most attention too
#---------------------------------------------------------------------------------
processWaveformQuestion <- function(data)
{
    question_responses <- NULL
    
    # rename column
    colnames(data)[2] <- "question"
    
    participants <- extractUniqueColData(data, "ParticipantName")
    for(i in 1:length(participants))
    {
        # add participants response to vector
        current_participant <- data[data$ParticipantName == participants[i], ]
        answer <- extractUniqueColData(current_participant, "question")
        question_responses <- c(question_responses, answer)
    }
    
    # output histogram
    print(barplot(table(question_responses), xlab = "Waveform component", ylab = "Count"))
}

#---------------------------------------------------------------------------------
# FUNCTION:     calculateProportionOfTime()
# INPUT:        data.frame, BOOL
# OUTPUT:       data.frame
# DESCRIPTION:  Calculate the proportion of time each participants has spent
#               looking at either P wave or rest of waveform and overall for
#               the entire rhythm strip
#---------------------------------------------------------------------------------
calculateProportionOfTime <- function(data, show_stats = TRUE)
{
    prop_data <- NULL
    prop_p_rows <- NULL
    prop_qrs_rows <- NULL
    participants <- getUniqueData(data, "ParticipantName")
    
    for(i in 1:length(participants))
    {
        # get the current participant
        current_participant <- data[data$ParticipantName == participants[i], ]
        
        for(j in 1:nrow(current_participant))
        {
            # get the current row
            current_row <- current_participant[j, ]
            
            # find all P or QRS
            if(toString(current_row$AOI) %in% c("p1", "p2", "p3", "p4", "p5"))
                prop_p_rows <- rbind(prop_p_rows, current_row)
            
            if(toString(current_row$AOI) %in% c("qrs1", "qrs2", "qrs3", "qrs4", "qrs5"))
                prop_qrs_rows <- rbind(prop_qrs_rows, current_row)
            
        }
      
        # calculate aggregate values per participant
        duration_p <-  sum(prop_p_rows$GazeEventDuration)
        duration_qrs <- sum(prop_qrs_rows$GazeEventDuration)
        total_time <- duration_p + duration_qrs
        prop_p <- (duration_p / total_time) * 100
        prop_qrs <- (duration_qrs / total_time) * 100
        
        # build data frame
        prop_data <- rbind(prop_data, list(participant = participants[i],
                                           stimuli = toString(current_participant[1, "MediaName"]),
                                           duration_p = duration_p,
                                           duration_qrs = duration_qrs,
                                           prop_p = prop_p,
                                           prop_qrs = prop_qrs))
        
    }
    
    # convert back to data frame
    prop_data <- as.data.frame(prop_data)
    prop_data$participant <- unlist(prop_data$participant)
    prop_data$stimuli <- unlist(prop_data$stimuli)
    prop_data$duration_p <- unlist(prop_data$duration_p)
    prop_data$duration_qrs <- unlist(prop_data$duration_qrs)
    prop_data$prop_p <- unlist(prop_data$prop_p)
    prop_data$prop_qrs <- unlist(prop_data$prop_qrs)
    
    # output descriptive stats (M, SD and Mdn)
    if(show_stats)
    {
        cat("\n Fixation duration (p): M = ", mean(prop_data$duration_p), " SD = ", sd(prop_data$duration_p))
        cat("\n Fixation duration (qrs): M = ", mean(prop_data$duration_qrs), " SD = ", sd(prop_data$duration_qrs))
        cat("\n Proportion (p): M = ", mean(prop_data$prop_p), " SD = ", sd(prop_data$prop_p), " Mdn = ", median(prop_data$prop_p))
        cat("\n Proportion (qrs): M = ", mean(prop_data$prop_qrs), " SD = ", sd(prop_data$prop_qrs), "Mdn = ", median(prop_data$prop_qrs))
    }
    return(prop_data)
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
    #common_source_data <- paste0(getwd(), "/Final-PHD-analysis/experiment_history.R")
    #source(common_source_data) 
    
    # open the data files
    eye_tracking_data <- openDataFile(paste0(getwd(), "/Final-PHD-analysis/data/cond3.csv"))
    AOI_locations <- openDataFile(paste0(getwd(), "/Final-PHD-analysis/data/AOI_RS.csv"))
    
    # show histogram of responses to waveform question
    #processWaveformQuestion(eye_tracking_data)
    
    STE <- eye_tracking_data[eye_tracking_data$MediaName == "STE_scaled.PNG", ]
    NSR <- eye_tracking_data[eye_tracking_data$MediaName == "NSR_scaled.PNG", ]
    
    # get hits
    hit_data_STE <- getAOIHits(STE, AOI_locations)
    hit_data_NSR <- getAOIHits(NSR, AOI_locations)
    
    # remove this participant as wasn't in both groups
    hit_data_NSR <- subset(hit_data_NSR, ParticipantName != "P27F")
    
    # calculate proportion of total time spent looking 
    prop_STE <- calculateProportionOfTime(hit_data_STE, TRUE)
    prop_NSR <- calculateProportionOfTime(hit_data_NSR, TRUE)
    
    # switch of scientific notation for p-values
    options(scipen = 999)
    
    # do paired t tests to look at proportion of time in each strip
    print(t.test(prop_STE$prop_p, prop_NSR$prop_p, paired = TRUE))
    print(t.test(prop_STE$prop_qrs, prop_NSR$prop_qrs, paired = TRUE))
    
}

# run main
main()
