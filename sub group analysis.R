##################################################################################
# NAME:         sub group analysis.R
# AUTHOUR:      Alan Davies
# DATE:         25/04/2018
# INSTITUTION:  Interaction Analysis and Modelling Lab (IAM), University of Manchester
# DESCRIPTION:  
#               
#               
#               
#               
##################################################################################
library(tidyverse)
library(crayon)

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
# FUNCTION:     getAccuracyData()
# INPUT:        data.frame, int
# OUTPUT:       data.frame
# DESCRIPTION:  Get the accuracy data for the selected condition
#               1 = history, 2 = no-history
#---------------------------------------------------------------------------------
getAccuracyData <- function(accuracy, selected_condition)
{
  participants <- getUniqueData(accuracy, "Participant")
  participants <- participants[-1]
  
  condition <- accuracy %>% 
    filter(Condition == selected_condition & Participant != "P01F") %>% 
    select(c(1, 3:ncol(accuracy)))
  
  condition_long <- condition %>% select(2:ncol(condition)) %>% gather(`Answer 1`:`Answer 9`, accuracy)
  colnames(condition_long)[1] <- "stimuli"
  
  condition_long$stimuli <- recode(condition_long$stimuli, `Answer 1` = "anterolateral STEMI",
                                   `Answer 2` = "LBBB",
                                   `Answer 3` = "lateral STEMI",
                                   `Answer 4` = "AF",
                                   `Answer 5` = "RBBB",
                                   `Answer 6` = "inferior STEMI and AF",
                                   `Answer 7` = "anterior STEMI",
                                   `Answer 8` = "high lateral STEMI",
                                   `Answer 9` = "inferolateral STEMI")
  
  condition_long <- cbind(rep(participants, 9), condition_long)
  colnames(condition_long)[1] <- "participant"
  return(condition_long)
}

#---------------------------------------------------------------------------------
# FUNCTION:     recodeStimuli()
# INPUT:        data.frame
# OUTPUT:       vector
# DESCRIPTION:  Return recoded stimuli vector 
#               
#---------------------------------------------------------------------------------
recodeStimuli <- function(data_vector)
{
    return(recode(data_vector, "anterolateralSTEMI.png" = "anterolateral STEMI",
                               "LBBB.png" = "LBBB",
                               "lateralSTEMI.png" = "lateral STEMI",
                               "AF.png" = "AF",
                               "RBBB.png" = "RBBB",
                               "inferiorSTEMIandAF.png" = "inferior STEMI and AF",
                               "anteriorSTEMI.png" = "anterior STEMI",
                               "highlateralSTEMI.png" = "high lateral STEMI",
                               "inferolateralSTEMI.png" = "inferolateral STEMI"))
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
    settings$num_perms <- readline("Enter number of permutations: ")  
    settings$AOI_types <- readline("(1) Top-down (from file) AOIs, (2) bottom-up (grid) AOIs: ")
    settings$echo_locations <- readline("Echo locations [show AOI sequences] (y/n): ")
    settings$density_plots <- readline("Output distribution plots (y/n): ")
    return(settings)
}

#---------------------------------------------------------------------------------
# FUNCTION:     convertToMatrix()
# INPUT:        data.frame, BOOL
# OUTPUT:       matrix
# DESCRIPTION:  Returns matrix representing probabilty of AOI transitions
#               
#---------------------------------------------------------------------------------
convertToMatrix <- function(data)
{
    tmp0 <- data
    tmp <- tmp0 %>% group_by(ParticipantName) %>% mutate(to = lead(AOI))
    tmp2 <- tmp[complete.cases(tmp), ]
    with(tmp2, table(AOI, to))
    out_mat <- as.matrix(with(tmp2, table(AOI, to)))
    
    # produce probability matrix
    prob_mat <- out_mat / rowSums(out_mat)
    return(prob_mat)
}

#---------------------------------------------------------------------------------
# FUNCTION:     getAccuracyGroup()
# INPUT:        data.frame, data.frame, int
# OUTPUT:       data.frame
# DESCRIPTION:  Get accuracy group and combine with ET data
#               
#---------------------------------------------------------------------------------
getAccuracyGroup <- function(answer_data, eye_data, accuracy_group)
{
    filtered_answers <- NULL
    current_stim <- NULL
    final_data <- NULL
 
    study_stimuli <- getUniqueData(answer_data, "stimuli")
    
    # recode the stimuli names
    eye_data$MediaName <- recodeStimuli(eye_data$MediaName)
    
    for(i in 1:length(study_stimuli))
    {
        # for each stimuli get all the people based on accuracy grouping and combine
        filtered_answers <- answer_data %>% filter(accuracy == accuracy_group & stimuli == study_stimuli[i]) %>% select(c(1, 2, 3))
        participants <- filtered_answers$participant 
        current_stim <- eye_data %>% filter(MediaName == study_stimuli[i] & ParticipantName %in% participants) %>% select(c(1, 3:7, 11:13)) 
        final_data <- rbind(final_data, current_stim)
    }
    
    # rename cols
    colnames(final_data)[3] <- "MediaPosX"
    colnames(final_data)[4] <- "MediaPosY"
    colnames(final_data)[8] <- "PosX"
    colnames(final_data)[9] <- "PosY"

    return(final_data)
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

    echo <- ifelse(settings$echo_locations == "y", TRUE, FALSE) 
    if(echo) cat(blue("\n Computing locations:\n\n"))
    AOI_locations <- AOI_locations[AOI_locations$condition == condition, ]
    stimuli <- getUniqueData(AOI_locations, "ECG")
    
    for(i in 1:length(stimuli))
    {
        # get the data for the AOI and stimuli
        if(echo) cat(blue("\n\nStimuli: ", stimuli[i]), "\n\n")
        AOI_stimuli <- AOI_locations[AOI_locations$ECG == stimuli[i], ]
        current_stimuli <- data[data$MediaName == stimuli[i], ]
        
        for(j in 1:nrow(current_stimuli))
        {
            # detect fixations within AOIs
            current_row <- current_stimuli[j, ]
            AOI <- c(AOI, detectHit(AOI_stimuli, current_row, echo))
        }
    }
    
    # append the AOI column to data and remove NA fixations from data
    data <- cbind(data, AOI)
    data <- data[complete.cases(data[, "AOI"]), ]
    
    return(data)
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
# FUNCTION:     generateDensityPlot()
# INPUT:        double, vector, list, list
# OUTPUT:       void
# DESCRIPTION:  Output density plot (distribution plot) for the results of the
#               permutation test. Purple line indicates primary group distance
#---------------------------------------------------------------------------------
generateDensityPlot <- function(distance_result, shuffled_distances, stimuli_data, settings)
{
    density_plot <- density(shuffled_distances)
    plot(density_plot, type = "n", main = stimuli_data$label, xlab = "Hellinger Distance", panel.first = grid())
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
    df <- data.frame(Hd = report_args$distance, p.value = report_args$pvalue, Permutations = settings$num_perms)
    print(df) 
    cat("\nGroup 1 (n):", report_args$groupsize1)
    cat("\nGroup 2 (n):", report_args$groupsize2, "\n\n")  
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
        participants <- as.data.frame(data[!duplicated(data$ParticipantName), "ParticipantName"])
        colnames(participants) <- "participant"
        
        # get a random subset of participants and store in group 1 whatever is left put in second group
        group1_participants <- participants[sample(unique(nrow(participants)), group1_size), ]
        group2_participants <- as.vector(participants[!(participants$participant %in% group1_participants), ])   
        
        # get the selected groups data
        group1 <- data[which(data$participant %in% group1_participants), ]
        group2 <- data[which(data$participant %in% group2_participants), ]
        
        # build Markov chains
        chain_1 <- convertToMatrix(group1)
        chain_2 <- convertToMatrix(group2)
        
        # remove any NaN
        chain_1[is.na(chain_1)] <- 0
        chain_2[is.na(chain_2)] <- 0
        
        distance_results <- calculateHellingerDistance(chain_1, chain_2)
        
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
# FUNCTION:     runPermutations()
# INPUT:        data.frame, data.frame, list
# OUTPUT:       void
# DESCRIPTION:  Run permutation tests 
#                    
#---------------------------------------------------------------------------------
runPermutations <- function(group1, group2, settings)
{
    report_args <- list()
    
    for(i in 1:length(getStimuliList()))
    {
        current_stimuli_metadata <- getExperimentSetupData(getStimuliList()[i])
        chain_1 <- group1 %>% filter(MediaName == current_stimuli_metadata$ref_name) %>% select(c("ParticipantName", "AOI"))
        chain_2 <- group2 %>% filter(MediaName == current_stimuli_metadata$ref_name) %>% select(c("ParticipantName", "AOI"))
        combined_data <- rbind(chain_1, chain_2)
        
        # extract num people in both groups
        num_people_in_group <- length(getUniqueData(chain_1, "ParticipantName"))
        remaining_people <- length(getUniqueData(chain_2, "ParticipantName"))
        
        # convert to matrix form
        chain_1 <- convertToMatrix(chain_1)
        chain_2 <- convertToMatrix(chain_2)  
        
        # remove any NA's from the chains
        chain_1[is.na(chain_1)] <- 0
        chain_2[is.na(chain_2)] <- 0
        
        # calculate distance and perumutations
        distance_result <- calculateHellingerDistance(chain_1, chain_2)
        shuffled_distances <- randomPermutations(combined_data, settings, num_people_in_group, current_stimuli_metadata)
        
        if(settings$density_plots == "y")
            generateDensityPlot(distance_result, shuffled_distances, current_stimuli_metadata, settings)
        
        # build report arguments
        report_args$distance <- distance_result
        report_args$pvalue <- calculatePvalue(distance_result, shuffled_distances) 
        report_args$groupsize1 <- num_people_in_group
        report_args$groupsize2 <- remaining_people
        
        # generate a textual report
        generateReport(settings, current_stimuli_metadata, report_args)
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
  # eps values (radius) from DBSCAN
  gridsizes <- c(30, 30, 30, 40, 30, 30, 40, 30, 30) * 2 
  
  if(interactive()){
      settings <- promptForSettings()
  }
  
  # load additional source and data files with stimuli metadata
  common_source_data <- paste0(getwd(), "/Final-PHD-analysis/experiment_history.R")
  data_files <- paste0(getwd(), "/Final-PHD-analysis/load_data_files.R")
  source(common_source_data) 
  source(data_files) 
  
  # recode the stimuli names in the AOI locations
  AOI_locations$ECG <- recodeStimuli(AOI_locations$ECG)
  #stimuli_list <- getUniqueData(AOI_locations, "ECG")
  
  # get the accuracy data (cond 1 or 2)
  history_data <- getAccuracyData(answers, 1)
  no_history_data <- getAccuracyData(answers, 2)
  
  # generate accuracy by condition sub groups
  history_correct <- getAccuracyGroup(history_data, condition_1, 1)
  history_incorrect <- getAccuracyGroup(history_data, condition_1, 0)
  no_history_correct <- getAccuracyGroup(no_history_data, condition_2, 1)
  no_history_incorrect <- getAccuracyGroup(no_history_data, condition_2, 0)
  
  # generate AOIs
  if(settings$AOI_types == "1"){
    history_correct <- getAOIHits(history_correct, AOI_locations, 1, settings)
    history_incorrect <- getAOIHits(history_incorrect, AOI_locations, 1, settings)
    no_history_correct <- getAOIHits(no_history_correct, AOI_locations, 2, settings)
    no_history_incorrect <- getAOIHits(no_history_incorrect, AOI_locations, 2, settings)
  } else {
      
  }
  
  # run permutation tests
  runPermutations(history_correct, history_incorrect, settings)
  runPermutations(no_history_correct, no_history_incorrect, settings)
 
  #utils::View(AOI_locations)
  #'utils::View(history_correct)
  #utils::View(history_incorrect)
  #utils::View(no_history_correct)
  #utils::View(no_history_incorrect)
  
  # TODO: Look for use of getAOIHits() function in MC_analysis.R
}

# run main
main()