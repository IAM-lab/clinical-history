##################################################################################
# NAME:         standard_analysis.R
# AUTHOUR:      Alan Davies
# DATE:         31/10/2017
# INSTITUTION:  Interaction Analysis and Modelling Lab (IAM), University of Manchester
# DESCRIPTION:  Standard statistical analysis techniques for analysis of 
#               second ECG experiment.
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
    package.args <- c("ggplot2", "binom", "tidyr", "reshape2", "plyr", "nlme", "lme4",
                      "stargazer", "lmeresampler", "PMCMR", "data.table", "plotrix") 
                    
    loadPackages(package.args)
    
    # set wd
    setwd(path.expand("~"))
}

#---------------------------------------------------------------------------------
# FUNCTION:     loadECGAnswers()
# INPUT:        BOOL
# OUTPUT:       data.frame
# DESCRIPTION:  Returns a data frame for the answers to the ECGs given by participants
#               as either coded or raw responses
#---------------------------------------------------------------------------------
loadECGAnswers <- function(coded = TRUE)
{
    to_get <- ifelse(coded, "answers_coded", "anaswers")
    data_path <- paste0(getwd(), "/Final-PHD-analysis/data/gs/", to_get, ".csv")
    data <- read.csv(data_path, header = TRUE, na.strings = c(" ", "NA", "-"))
    return(data)
}

#---------------------------------------------------------------------------------
# FUNCTION:     loadGivenAnswers()
# INPUT:        void
# OUTPUT:       data.frame
# DESCRIPTION:  Returns a data frame for the answers to the ECGs given by participants
#       
#---------------------------------------------------------------------------------
loadGivenAnswers <- function()
{
    data_path <- paste0(getwd(), "/Final-PHD-analysis/data/gs/given_answers.csv")
    data <- read.csv(data_path, header = TRUE, na.strings = c(" ", "NA", "-"))
    return(data)
}

#---------------------------------------------------------------------------------
# FUNCTION:     loadEyeTrackingData()
# INPUT:        void
# OUTPUT:       list
# DESCRIPTION:  Returns a list of data frames for the 3 different files containing
#               the eye-tracking data for the 3 conditions (history, no-history and
#               rhythm-strip)
#---------------------------------------------------------------------------------
loadEyeTrackingData <- function()
{
    data_files <- list()
    data_path <- paste0(getwd(), "/Final-PHD-analysis/data/")
    data_files[["cond_1"]] <- read.csv(paste0(data_path, "cond1.csv"), header = TRUE, na.strings = c(" ", "NA", "-"))
    data_files[["cond_2"]] <- read.csv(paste0(data_path, "cond2.csv"), header = TRUE, na.strings = c(" ", "NA", "-"))
    data_files[["cond_3"]] <- read.csv(paste0(data_path, "cond3.csv"), header = TRUE, na.strings = c(" ", "NA", "-"))
    return(data_files)
}

#---------------------------------------------------------------------------------
# FUNCTION:     loadAOILocations()
# INPUT:        void
# OUTPUT:       data.frame
# DESCRIPTION:  Loads the AOI location and dimension data
#               
#               
#---------------------------------------------------------------------------------
loadAOILocations <- function()
{
    data_path <- paste0(getwd(), "/Final-PHD-analysis/data/AOI.csv")
    data <- read.csv(data_path, header = TRUE, na.strings = c(" ", "NA", "-"))
    return(data) 
}

#---------------------------------------------------------------------------------
# FUNCTION:     loadSurveyData()
# INPUT:        void
# OUTPUT:       data.frame
# DESCRIPTION:  Return the formatted survey results as a data frame 
#               
#---------------------------------------------------------------------------------
loadSurveyData <- function()
{
    data <- read.csv(paste0(getwd(), "/Final-PHD-analysis/data/survey results.csv"), 
                     header = TRUE, na.strings = c(" ", "NA"))
    
    # rename the columns
    colnames(data) <- c("timestamp", "difficulty", "confidence", "student", 
                        "role", "roleduration", "interpretationfreq", "ratedexperience",
                        "hourstraining", "trainingformat", "system", "systemtype",
                        "systemorigin", "systemchanged", "howchanged", "whychanged",
                        "leads", "calibration", "autouseful", "autoaccurate", 
                        "difficultywithout", "difficultywith", "pid", "gender", "age")
    
    return(data)
}

#---------------------------------------------------------------------------------
# FUNCTION:     accuracyPerStimulus()
# INPUT:        data.frame
# OUTPUT:       void
# DESCRIPTION:  Displays the number of correct and incorrect per condition
#               per ECG     
#---------------------------------------------------------------------------------
accuracyPerStimulus <- function(data)
{
    cond_1 <- data[data$Condition == 1, ]
    cond_2 <- data[data$Condition == 2, ]
    
    for(i in 1:2)
    {
        if(i == 1){
            condition <- cond_1
        } else {
            condition <- cond_2
        }
        condition <- condition[ ,c(3:ncol(condition))]
        cat("\n\nCondition: ", i, " total = ", nrow(condition), "\n\n")
        for(j in 1:ncol(condition))
        {
            cat("\n", getAnswerNames()[j], "\t\t Correct = ", 
                sum(as.vector(condition[ ,j])), ", Incorrect = ", 
                nrow(condition) - sum(as.vector(condition[ ,j])))
        }
    }
}

#---------------------------------------------------------------------------------
# FUNCTION:     accuracyPerParticipant()
# INPUT:        data.frame
# OUTPUT:       void
# DESCRIPTION:  Create boxplot for history and no-histroy groups showing the 
#               overall % accuracy accross all stimuli. Check for differences
#               withb MWU test
#---------------------------------------------------------------------------------
accuracyPerParticipant <- function(data)
{
    accuracy <- NULL
    cond_1 <- data[data$Condition == 1, ]
    cond_2 <- data[data$Condition == 2, ]
    
    participants <-  as.vector(data[!duplicated(data$Participant), "Participant"])
    participants <- participants[-1]
    
    for(i in 1:2)
    {
        if(i == 1){
            condition <- cond_1
        } else {
            condition <- cond_2
        }
        for(j in 1:length(participants))
        {
            # build data frame
            current_participant <- condition[condition$Participant == participants[j], c(3:ncol(condition))]
            accuracy <- rbind(accuracy, list(Participant = participants[j],
                                             pc_accuracy = (rowSums(current_participant) / 9) * 100,
                                             condition = i))
        }
    }
    #unlist data
    accuracy <- as.data.frame(accuracy)
    accuracy$Participant <- unlist(accuracy$Participant)
    accuracy$pc_accuracy <- unlist(accuracy$pc_accuracy) 
    accuracy$condition <- unlist(accuracy$condition)
    
    # box plot for accuracy with lablled outliers 
    accuracy_plot <- ggplot(accuracy, aes(x = factor(condition), y = pc_accuracy, fill = factor(condition))) +
                     geom_boxplot(alpha = 0.7) + 
                     geom_text(aes(label = ifelse(pc_accuracy < 20, as.character(Participant), '')), 
                               position = position_jitter(width = .06)) +
                     theme(text = element_text(size = 15), axis.text = element_text(colour = "black", size = 15),
                     legend.position = "bottom") +
                     labs(x = "Condition", y = "Interpretation accuracy (%)") +
                     scale_x_discrete(labels = c("History", "No-history"))
    
    print(accuracy_plot)

    print(boxplot(accuracy$pc_accuracy, ylab = "Accuracy (%)", xlab = "Participants"))
    cat("\n\nMean accuracy (all participants) = ", mean(accuracy$pc_accuracy), "SD = ", sd(accuracy$pc_accuracy), "\n")
    
    # Compare the overall % accuracy of both groups
    hx_accuracy <- as.vector(accuracy[accuracy$condition == 1, "pc_accuracy"])
    nohx_accuracy <- as.vector(accuracy[accuracy$condition == 2, "pc_accuracy"])
    print(wilcox.test(hx_accuracy, nohx_accuracy, paired = TRUE))
    
    cat("\n\nHx accuracy (mean) = ", mean(hx_accuracy), " SD = ", sd(hx_accuracy), "\n")
    cat("nno Hx accuracy (mean) = ", mean(nohx_accuracy), " SD = ", sd(nohx_accuracy), "\n")
    
    given_answers <- loadGivenAnswers()
    saw_HPC <- getPresentationSequence(given_answers, "1-2")
    saw_noHPC <- getPresentationSequence(given_answers, "2-1")
    
    HPC <- accuracy[which(accuracy$Participant %in% saw_HPC), ]
    NHPC <- accuracy[which(accuracy$Participant %in% saw_noHPC), ]
    
    cat("\n\n Saw Hx first accuracy (mean) = ", mean(HPC$pc_accuracy), " SD = ", sd(HPC$pc_accuracy), "\n")
    cat("n Saw no Hx first accuracy (mean) = ", mean(NHPC$pc_accuracy), " SD = ", sd(NHPC$pc_accuracy), "\n")
    print(wilcox.test(HPC$pc_accuracy, NHPC$pc_accuracy))
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
# FUNCTION:     accuracyPerStimuli()
# INPUT:        data.frame
# OUTPUT:       void
# DESCRIPTION:  Produce barplot per condition for the proportion of people making
#               correct interpretations per ECG
#---------------------------------------------------------------------------------
accuracyPerStimuli <- function(data)
{
    accuracy <- NULL
    cond_1 <- data[data$Condition == 1, ]
    cond_2 <- data[data$Condition == 2, ]
    
    # remove P1
    data <- data[-1, ]
    
    for(i in 1:2)
    {
        if(i == 1){
            condition <- cond_1
        } else {
            condition <- cond_2
        }
        condition <- condition[ ,c(3:ncol(condition))]
        for(j in 1:ncol(condition))
        {
           
            # build data frame
            current_stimuli <- as.vector(condition[ ,j])
            sum_stim <- sum(current_stimuli)
            num_stim <- length(current_stimuli)

            # calculate proportional CI using Wilson method            
            CI <- binom.confint(sum_stim, 29, conf.level = 0.95, methods = "wilson") 
            accuracy <- rbind(accuracy, list(stimuli =  getAnswerNames()[j],
                                             pc_correct = (sum_stim / num_stim) * 100,
                                             CI_wilson_lower = CI$lower * 100,
                                             CI_wilson_upper = CI$upper * 100,
                                             cond = i))
        }
    }
    #unlist data
    accuracy <- as.data.frame(accuracy)
    accuracy$stimuli <- unlist(accuracy$stimuli)
    accuracy$pc_correct <- unlist(accuracy$pc_correct) 
    accuracy$CI_wilson_lower <- unlist(accuracy$CI_wilson_lower)
    accuracy$CI_wilson_upper <- unlist(accuracy$CI_wilson_upper)
    accuracy$cond <- unlist(accuracy$cond)
  
    # produce grouped bar plot for condition and accuracy
    accuracy_plot_1 <- ggplot(accuracy, aes(x = factor(stimuli), y = pc_correct, fill = factor(cond))) +
                             geom_bar(stat = "identity", position = position_dodge(0.9)) +
                             theme(text = element_text(size = 15), 
                                   axis.text.x = element_text(angle = 70, hjust = 1, colour = "black"), 
                                   axis.text.y = element_text(colour = "black")) +
                             geom_errorbar(aes(ymin = CI_wilson_lower, ymax = CI_wilson_upper), 
                                   position = position_dodge(0.9), width = 0.25) +
                             labs(x = "Stimuli (ECG)", y = "Overall interpretation accuracy (%)") +
                             scale_y_continuous(limits = c(0, 100)) 
    
    print(accuracy_plot_1)
}

#---------------------------------------------------------------------------------
# FUNCTION:     surveyResultsBoxplots()
# INPUT:        data.frame
# OUTPUT:       void
# DESCRIPTION:  Make a horizontal stacked barplot, showing the % of scores for
#               a given rating (1-5 Likert scale) for survey questions
#---------------------------------------------------------------------------------
surveyResultsBoxplots <- function(data)
{
    question_data <- NULL
    
    # get all the likert scale questions
    data <- data[ ,c("difficulty", "confidence", "ratedexperience", "autouseful",
                     "autoaccurate", "difficultywithout", "difficultywith")]
    
    # store the questions 
    questions <- c("How would you rate the difficulty of the on-screen\n ECG interpretation task you just took part in?",
                   "How confident do you typically feel when interpreting ECGs?",
                   "How would you rate your experience level with ECG interpretation?",
                   "How useful do you find automated ECG interpretation outputs?",
                   "How accurate do you believe automated ECG interpretation to be?",
                   "How would you rate the difficulty of the interpretation\n task WITHOUT the history of the patient's presenting complaint?",
                   "How would you rate the difficulty of the interpretation task\n WITH the history of the patient's presenting complaint?")
    
    # count number of occurances of each rating for the questions
    for(i in 1:ncol(data))
    {
        current_col <- as.vector(data[ ,i])
        question_data <- rbind(question_data, list(question = questions[i],
                                                   "1" = (length(which(current_col == 1)) / length(current_col)) * 100,
                                                   "2" = (length(which(current_col == 2)) / length(current_col)) * 100,
                                                   "3" = (length(which(current_col == 3)) / length(current_col)) * 100,
                                                   "4" = (length(which(current_col == 4)) / length(current_col)) * 100,
                                                   "5" = (length(which(current_col == 5)) / length(current_col)) * 100))
    }
    
    # unlist data
    question_data <- as.data.frame(question_data)
    question_data$question <- unlist(question_data$question)
    question_data[["1"]] <- unlist(question_data[["1"]])
    question_data[["2"]] <- unlist(question_data[["2"]])
    question_data[["3"]] <- unlist(question_data[["3"]])
    question_data[["4"]] <- unlist(question_data[["4"]])
    question_data[["5"]] <- unlist(question_data[["5"]])
    
    # convert data format into summary of values for each question
    question_data <- as.data.frame(melt(question_data, id = c("question")))
    
    # plot data
    stacked_plot <- ggplot(question_data, aes(x = factor(question), y = value, fill = factor(variable))) + 
        geom_bar(stat = "identity", width = .5) +
        theme(aspect.ratio = 0.9, axis.text = element_text(size = 12, colour = "black"), panel.background = element_blank()) +
        scale_fill_brewer() +
        geom_text(aes(label = ifelse(value > 0, paste0(round(value),"%"), "")), 
                  position = position_stack(vjust = .5)) +
        coord_flip() +
        labs(x = "", y = "") +
        theme(legend.position = "bottom") +
        guides(fill = guide_legend(title = ""))
    
    print(stacked_plot)
}


#---------------------------------------------------------------------------------
# FUNCTION:     getStimuli()
# INPUT:        data.frame
# OUTPUT:       vector
# DESCRIPTION:  Returns a vector of unique stimuli names
#               
#---------------------------------------------------------------------------------
getStimuli <- function(data)
{
    data <- data[!duplicated(data$MediaName), ]
    data <- data[data$MediaName != "", ]
    return(as.vector(data$MediaName))
}

#---------------------------------------------------------------------------------
# FUNCTION:     extractStimulus()
# INPUT:        data.frame, String
# OUTPUT:       data.frame
# DESCRIPTION:  Retunrs data for a specific stimulus 
#               
#---------------------------------------------------------------------------------
extractStimulus <- function(data, stimuli_name)
{
    return(data[data$MediaName == stimuli_name, ])
}

#---------------------------------------------------------------------------------
# FUNCTION:     filterStimuliTransitionFixation()
# INPUT:        data.frame
# OUTPUT:       data.frame
# DESCRIPTION:  For each participant remove the first fixation for all subsequent
#               stimului except the first one and recombine. This is done to remove
#               the first fixation from a new stimulus as it is still related to where
#               the participant was looking on the previous stimulus
#---------------------------------------------------------------------------------
filterStimuliTransitionFixation <- function(data, stimuli)
{
    filtered_data <- data[FALSE, ]
    for(i in 1:length(stimuli))
    {
        current_stimuli <- extractStimulus(data, stimuli[i])
        if(i > 1) current_stimuli <- current_stimuli[-1, ]
        filtered_data <- rbind(filtered_data, current_stimuli)
    }
    return(filtered_data)
}

#---------------------------------------------------------------------------------
# FUNCTION:     getUniqueFixationData()
# INPUT:        data.frame
# OUTPUT:       data.frame
# DESCRIPTION:  Returns filtered fixation data (removes any repeated data)
#              
#---------------------------------------------------------------------------------
getUniqueFixationData <- function(data)
{
    filtered_data <- data[FALSE, ]
    stimuli <- getStimuli(data)
    participants <-  as.vector(unlist(data[!duplicated(data$ParticipantName), "ParticipantName"], use.names = FALSE))
    
    # get all the unique fixation data  
    for(i in 1:length(participants))
    {
        participant_data <- data[data$ParticipantName == participants[i], ]
        participant_data <- participant_data[!duplicated(participant_data$FixationIndex), ]
        participant_data <- filterStimuliTransitionFixation(participant_data, stimuli) 
        filtered_data <- rbind(filtered_data, participant_data)
    }
    return(filtered_data)
}

#---------------------------------------------------------------------------------
# FUNCTION:     filterSequenceData()
# INPUT:        void
# OUTPUT:       data.frame
# DESCRIPTION:  Returns filtered sequence data combined with accuracy
#              
#---------------------------------------------------------------------------------
filterSequenceData <- function(data)
{
    AOI_vector <- NULL
    processed_data <- NULL
    final_data <- NULL
    
    # get fixation data
    data <- data[data$GazeEventType == "Fixation", ]
    data[is.na(data)] <- 0
   
    for(i in 1:nrow(data))
    {
        # get hit data by row
        current_row <- data[i, c(15:ncol(data))]
        hit <- which(current_row %in% 1)  
        if(length(hit) > 0) 
        {
            AOI <- substr(colnames(current_row)[hit], 5, 5)
            final_data <- rbind(final_data, list(ParticipantName = as.vector(data[i, "ParticipantName"]),
                                                 MediaName = as.vector(data[i, "MediaName"]),
                                                 FixationIndex = as.vector(data[i, "FixationIndex"]),
                                                 GazeEventType = as.vector(data[i, "GazeEventType"]),
                                                 GazeEventDuration = as.vector(data[i, "GazeEventDuration"]),
                                                 AOI = AOI))
        }
    }
    
    # unlist data
    final_data <- as.data.frame(final_data)

    return(final_data)
}

#---------------------------------------------------------------------------------
# FUNCTION:     buildScanPathSequencesPerStimuli()
# INPUT:        data.frame
# OUTPUT:       void
# DESCRIPTION:  Generate a data frame that contains the scanpaths, both
#               compressed and uncompressed, with lengths and accuracy
#---------------------------------------------------------------------------------
buildScanPathSequencesPerStimuli <- function(data)
{
    merged_data <- NULL
    
    for(i in 1:length(data))
    {
        # get the sequence data for each ECG
        filterd_seq <- filterSequenceData(data[[i]])
        participants <- as.vector(filterd_seq[!duplicated(filterd_seq$ParticipantName), "ParticipantName"])
        
        for(j in 1:length(participants))
        {
            # extract data per-participant
            current_participant <- filterd_seq[filterd_seq$ParticipantName == participants[j], ]
            current_stimuli <- toString(filterd_seq[1, "MediaName"])
            accuracy <- current_participant[1, "accuracy"] 
            scanpath <- paste(current_participant$AOI, collapse = "")
            collapsed <- gsub('([[:alpha:]])\\1+', '\\1', scanpath)
            
            # build a new data frame row by row
            merged_data <- rbind(merged_data, list(participant = participants[j],
                                                   stimuli = current_stimuli,
                                                   scanpath = scanpath,
                                                   scanpath_len = nchar(scanpath),
                                                   collapsed = collapsed,
                                                   collapsed_len = nchar(collapsed),
                                                   accuracy = accuracy))
        }
    }
    return(merged_data)
}

#---------------------------------------------------------------------------------
# FUNCTION:     combineAccuracyData()
# INPUT:        data.frame, data.frame
# OUTPUT:       data.frame
# DESCRIPTION:  Add the accuracy data to the eye-tracking data frame
#               
#---------------------------------------------------------------------------------
combineAccuracyData <- function(eye_tracking_data, answer_data)
{
    combined_data <- NULL
    
    # get the unique stimuli and participants lists
    stimuli <- extractUniqueColData(eye_tracking_data, "MediaName")
    participants <- extractUniqueColData(eye_tracking_data, "ParticipantName")
    
    # Rename the answer columns
    colnames(answer_data) <-  c("Participant", "Condition", "anterolateralSTEMI.png", "LBBB.png",
                                "lateralSTEMI.png", "AF.png", "RBBB.png", "inferiorSTEMIandAF.png",
                                "anteriorSTEMI.png", "highlateralSTEMI.png", "inferolateralSTEMI.png")
    
    for(k in 1:2)
    {
        # get condition at a time
        condition_data <- answer_data[answer_data$Condition == k, ]
        
        for(i in 1:length(participants))
        {
            # get the eye-tracking data a
            current_et_participant <- eye_tracking_data[eye_tracking_data$ParticipantName == participants[i], ]
            current_ans_participant <- condition_data[condition_data$Participant == participants[i], ]
        
            for(j in 1:length(stimuli))
            {
                # get the accuracy
                accuracy <- current_ans_participant[1, stimuli[j]]
                
                if(!is.na(accuracy))
                {
                    # if not NA then add the accuracy and condition data to the eye-tracking data
                    current_et_participant$accuracy <- rep(accuracy, nrow(current_et_participant))
                    current_et_participant$condition <- rep(k, nrow(current_et_participant))
                    combined_data <- rbind(combined_data, current_et_participant)
                }
            }
        }
    }

    return(combined_data)
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
# FUNCTION:     generateModels()
# INPUT:        data.frame
# OUTPUT:       void
# DESCRIPTION:  
#               
#---------------------------------------------------------------------------------
generateModels <- function(data)
{
    # load the given answers with the presentation sequence
    given_answers <- loadGivenAnswers()
    given_answers <- given_answers[-1, ]
    
    # get the people that saw each condition first
    history_participants <- as.vector(given_answers[given_answers$PresentationSequence == "1-2", "Participant"])
    nohistory_participants <- as.vector(given_answers[given_answers$PresentationSequence == "2-1", "Participant"])
    history_participants <- history_participants[!duplicated(history_participants)]
    nohistory_participants <- nohistory_participants[!duplicated(nohistory_participants)]
    clinical_hx <- data[which(data$ParticipantName %in% history_participants), ]
    noclinical_hx <- data[which(data$ParticipantName %in% nohistory_participants), ]
    
    # combibne the data removing the repeted mesure
    combined_data <- rbind(clinical_hx, noclinical_hx)
    combined_data <- combined_data[, c("ParticipantName", "MediaName", "GazeEventDuration", "AOI", "accuracy", "condition")]    
    combined_data <- as.data.frame(combined_data)
    combined_data$ParticipantName <- unlist(combined_data$ParticipantName)
    combined_data$MediaName <- unlist(combined_data$MediaName)
    combined_data$GazeEventDuration <- unlist(combined_data$GazeEventDuration)
    combined_data$AOI <- factor(unlist(combined_data$AOI))
    
    mod_1 <- lmer(accuracy ~ GazeEventDuration + (1|condition/MediaName), data = combined_data, REML = FALSE)
    mod_2 <- lmer(accuracy ~ GazeEventDuration + (1|condition), data = combined_data, REML = FALSE)
    mod_3 <- lmer(accuracy ~ GazeEventDuration + (1|condition) + (1|AOI), data = combined_data, REML = FALSE)
    
    print(summary(mod_3))
    print(stargazer(mod_3, type = "text", digits = 3, star.cutoffs = c(0.05, 0.01, 0.001), digit.separator = ""))
}

#---------------------------------------------------------------------------------
# FUNCTION:     resultsSurvey()
# INPUT:        data.frame
# OUTPUT:       void
# DESCRIPTION:  Output graphs and stats for the the remaining survey data
#               
#---------------------------------------------------------------------------------
resultsSurvey <- function(data)
{
    cat("\n Role duration: M=", mean(data$roleduration), " SD=", sd(data$roleduration), 
        " Mdn=", median(data$roleduration), " min=", min(data$roleduration), " max=", max(data$roleduration))
    
    int_freq <- data$interpretationfreq
    hours_freq <- data$hourstraining
    system_freq <- data$system
    changed_freq <- data$systemchanged
    calib_freq <- data$calibration
    
    p1 <- ggplot(data.frame(int_freq), aes(x = int_freq)) + 
        geom_bar() +
        theme(text = element_text(size = 15), 
              axis.text.x = element_text(angle = 70, hjust = 1, colour = "black"), 
              axis.text.y = element_text(colour = "black")) +      
        labs(x = "Frequency of interpretation", y = "Count")
    
    p2 <- ggplot(data.frame(hours_freq), aes(x = hours_freq)) + 
        geom_bar() +
        theme(text = element_text(size = 15), 
              axis.text.x = element_text(angle = 70, hjust = 1, colour = "black"), 
              axis.text.y = element_text(colour = "black")) +      
        labs(x = "Hours training", y = "Count")
    
    p3 <- ggplot(data.frame(system_freq), aes(x = system_freq)) + 
        geom_bar() +
        theme(text = element_text(size = 15), 
              axis.text.x = element_text(angle = 70, hjust = 1, colour = "black"), 
              axis.text.y = element_text(colour = "black")) +      
        labs(x = "Use a system", y = "Count")
    
    p4 <- ggplot(data.frame(changed_freq), aes(x = changed_freq)) + 
        geom_bar() +
        theme(text = element_text(size = 15), 
              axis.text.x = element_text(angle = 70, hjust = 1, colour = "black"), 
              axis.text.y = element_text(colour = "black")) +      
        labs(x = "Has the system changed", y = "Count")
    
    p5 <- ggplot(data.frame(calib_freq), aes(x = calib_freq)) + 
        geom_bar() +
        theme(text = element_text(size = 15), 
              axis.text.x = element_text(angle = 70, hjust = 1, colour = "black"), 
              axis.text.y = element_text(colour = "black")) +      
        labs(x = "Frequency of checking calibration", y = "Count")
    
    print(p1)
    print(p2)
    print(p3)
    print(p4)
    print(p5)
    
    #print(sum(int_freq == "Monthly"))
    #print(sum(int_freq == "Weekly"))
    
    lead_data <- NULL
    lead_names <- c("I", "II", "III", "aVR", "aVL", "aVF", "V1", "V2", "V3", "V4", "V5", "V6", "RS")
    lead_count <- c(10, 25, 14, 8, 9, 11, 23, 19, 17, 15, 19, 16, 19)
    
    for(i in 1:length(lead_names))
    {
        lead_data <- c(lead_data, rep(toString(lead_names[i]), lead_count[i]))
    }
    
    p6 <- ggplot(data.frame(lead_data), aes(x = lead_data)) + 
        geom_bar() +
        theme(text = element_text(size = 15), 
              axis.text.x = element_text(angle = 70, hjust = 1, colour = "black"), 
              axis.text.y = element_text(colour = "black")) +      
        labs(x = "Leads paid most attention", y = "Count") 
    
    print(p6)
    
    print(data.frame(table(hours_freq)))
    print(data.frame(table(system_freq))) 
    print(data.frame(table(changed_freq)))
    print(data.frame(table(calib_freq)))
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
# INPUT:        data.frame, data.frame
# OUTPUT:       data.frame
# DESCRIPTION:  Get the locations for each fixation within the stimulus 
#               
#---------------------------------------------------------------------------------
getAOIHits <- function(data, AOI_locations, condition) 
{
    AOI <- NULL
    AOI_locations <- AOI_locations[AOI_locations$condition == condition, ]
    
    stimuli <- getUniqueData(AOI_locations, "ECG")
    from_col <- which(colnames(data) == "ParticipantName")
    to_col <- which(colnames(data) == "SaccadicAmplitude")
    
    # extract data from 1 column to another & rename cols
    main_data <- data[c(from_col:to_col)]
    colnames(main_data)[4] <- "MediaPosX"
    colnames(main_data)[5] <- "MediaPosY"
    colnames(main_data)[12] <- "PosX"
    colnames(main_data)[13] <- "PosY"
    
    for(i in 1:length(stimuli))
    {
        # get the data for the AOI and stimuli
        AOI_stimuli <- AOI_locations[AOI_locations$ECG == stimuli[i], ]
        current_stimuli <- main_data[main_data$MediaName == stimuli[i], ]
        
        for(j in 1:nrow(current_stimuli))
        {
            # detect fixations within AOIs
            current_row <- current_stimuli[j, ]
            AOI <- c(AOI, detectHit(AOI_stimuli, current_row))
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
detectHit <- function(AOI_data, current_row)
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
            break
        }
    }
    return(AOI)
}

#---------------------------------------------------------------------------------
# FUNCTION:     fixationDurationPerStimuli()
# INPUT:        data.frame
# OUTPUT:       void
# DESCRIPTION:  Produce plots for condition groups for average FD, TFF and FC
#               
#---------------------------------------------------------------------------------
fixationDurationPerStimuli <- function(data, AOI_locations)
{
    AOI_data <- NULL
    AOIs <- LETTERS[1:15]
    history <- data[[1]]
    nohistory <- data[[2]]
    
    AOI_history <- getAOIHits(history, AOI_locations, 1)
    AOI_nohistory <- getAOIHits(nohistory, AOI_locations, 2)
    AOI_history$group <- rep(1, nrow(AOI_history))
    AOI_nohistory$group <- rep(2, nrow(AOI_nohistory))
    merged_data <- rbind(AOI_history, AOI_nohistory)
    
    stimuli <- extractUniqueColData(AOI_history, "MediaName")
    
    utils::View(merged_data)
    
    for(i in 1:2)
    {
        if(i == 1)
            current_data <- AOI_history
        else 
            current_data <- AOI_nohistory
        
        for(j in 1:length(stimuli))
        {
            current_stimuli <- current_data[current_data$MediaName == stimuli[j], ] 
            AOI_data <- rbind(AOI_data, list(stimuli = stimuli[j],
                                             FC = nrow(current_stimuli),
                                             avg_FD = mean(current_stimuli$GazeEventDuration),
                                             sd = sd(current_stimuli$GazeEventDuration),
                                             se = std.error(current_stimuli$GazeEventDuration),
                                             group = i))
        }
    }
    AOI_data <- as.data.frame(AOI_data)
    AOI_data$stimuli <- unlist(AOI_data$stimuli)
    AOI_data$FC <- unlist(AOI_data$FC)
    AOI_data$avg_FD <- unlist(AOI_data$avg_FD)
    AOI_data$sd <- unlist(AOI_data$sd)
    AOI_data$se <- unlist(AOI_data$se)
    AOI_data$group <- unlist(AOI_data$group)
    
    stimuli_plot <- ggplot(AOI_data, aes(x = factor(stimuli), y = avg_FD, fill = factor(group))) +
        geom_bar(stat = "identity", position = position_dodge(0.9)) +
        theme(text = element_text(size = 15), 
              axis.text = element_text(colour = "black", size = 15),
              axis.text.x = element_text(angle = 70, hjust = 1, colour = "black"), 
              axis.text.y = element_text(colour = "black")) +
        geom_errorbar(aes(ymin = avg_FD - se, ymax = avg_FD + se), position = position_dodge(0.9), width = 0.25) +
        labs(x = "ECG", y = "Mean fixation duration (ms)") +
        scale_x_discrete(labels = swapStimuliNames())

    print(stimuli_plot)
    
    pairwiseTestsPerStimuli(merged_data)
}

#---------------------------------------------------------------------------------
# FUNCTION:     pairwiseTestsPerStimuli()
# INPUT:        data.frame
# OUTPUT:       void
# DESCRIPTION:  Output pairwise tests per ECG using Bonferroni adjustment
#               
#---------------------------------------------------------------------------------
pairwiseTestsPerStimuli <- function(data)
{
    df <- NULL
    stimuli_data <- list()
    stimuli <- extractUniqueColData(data, "MediaName")
    participants <- extractUniqueColData(data, "ParticipantName")
    
    for(k in 1:2)
    {
        condition <- data[data$group == k, ]
        
        for(i in 1:length(participants))
        {
            current_participant <- condition[condition$ParticipantName == participants[i], ]
    
            for(j in 1:length(stimuli))
            {
                current_stimuli <- current_participant[current_participant$MediaName == stimuli[j], ]
                stimuli_data <- c(stimuli_data, sum(current_stimuli$GazeEventDuration))
            }
            if(!is.na(current_participant$group[1]))
            {
                df <- rbind(df, list(participant = participants[i],
                                     anterolateralSTEMI.png = stimuli_data[1],
                                     LBBB.png = stimuli_data[2],
                                     lateralSTEMI.png = stimuli_data[3],
                                     AF.png = stimuli_data[4],
                                     RBBB.png = stimuli_data[5],
                                     inferiorSTEMIandAF.png = stimuli_data[6],
                                     anteriorSTEMI.png = stimuli_data[7],
                                     highlateralSTEMI.png = stimuli_data[8],
                                     inferolateralSTEMI.png = stimuli_data[9],
                                     group = current_participant$group[1]))
            }
            stimuli_data <- NULL
        }
    }
    
    df <- as.data.frame(df)
    df$participant <- unlist(df$participant)
    df$anterolateralSTEMI.png = unlist(df$anterolateralSTEMI.png)
    df$LBBB.png = unlist(df$LBBB.png)
    df$lateralSTEMI.png = unlist(df$lateralSTEMI.png)
    df$AF.png = unlist(df$AF.png)
    df$RBBB.png = unlist(df$RBBB.png)
    df$inferiorSTEMIandAF.png = unlist(df$inferiorSTEMIandAF.png)
    df$anteriorSTEMI.png = unlist(df$anteriorSTEMI.png)
    df$highlateralSTEMI.png = unlist(df$highlateralSTEMI.png)
    df$inferolateralSTEMI.png = unlist(df$inferolateralSTEMI.png)
    df$group = unlist(df$group)
    
    # Bonferroni adjustment (alpha / stimuli)
    cat("\n Pairwise tests (Bonferroni adjustment = ", (0.05 / 9),")\n\n")
    
    group_1 <- df[df$group == 1, ]
    group_2 <- df[df$group == 2, ]
    
    group_1 <- group_1[ ,c(2:ncol(group_1))]
    group_2 <- group_2[ ,c(2:ncol(group_2))]
    
    for(i in 1:length(stimuli))
    {
        cat("\nStimuli = ", stimuli[i], "\n")
        history <- group_1[ ,i]
        nohistory <- group_2[ ,i]
        print(wilcox.test(history, nohistory, paired = TRUE))
    }
}

#---------------------------------------------------------------------------------
# FUNCTION:     fixationDurationPerLead()
# INPUT:        data.frame
# OUTPUT:       void
# DESCRIPTION:  
#               
#---------------------------------------------------------------------------------
fixationDurationPerLead <- function(data, AOI_locations)
{
    AOI_data <- NULL
    AOIs <- LETTERS[1:15]
    
    history <- data[[1]]
    nohistory <- data[[2]]
    
    AOI_history <- getAOIHits(history, AOI_locations, 1)
    AOI_nohistory <- getAOIHits(nohistory, AOI_locations, 2)
    
    stimuli <- extractUniqueColData(AOI_history, "MediaName")
    participants <- extractUniqueColData(AOI_history, "ParticipantName")
    
    for(i in 1:2)
    {
        if(i == 1)
            current_data <- AOI_history
        else 
            current_data <- AOI_nohistory
     
        for(k in 1:length(stimuli))
        {
            current_stimuli <- current_data[current_data$MediaName == stimuli[k], ] 
            
            #for(m in 1:length(participants))
            #{
            #    current_participant <- current_stimuli[current_stimuli$ParticipantName == participants[m], ]
                
                for(j in 1:length(AOIs))
                {
                    current_AOI <- current_stimuli[current_stimuli$AOI == AOIs[j], ]
                    AOI_data <- rbind(AOI_data, list(stimuli = stimuli[k],
                                                     AOI = AOIs[j],
                                                     FC = nrow(current_AOI),
                                                     avg_FD = mean(current_AOI$GazeEventDuration),
                                                     sd = sd(current_AOI$GazeEventDuration),
                                                     total_FD = sum(current_AOI$GazeEventDuration),
                                                     se = std.error(current_AOI$GazeEventDuration),
                                                     group = i))
                }
            #}
        }
    }
    AOI_data <- as.data.frame(AOI_data)
    AOI_data$stimuli <- unlist(AOI_data$stimuli)
    #AOI_data$participant <- unlist(AOI_data$participant)
    AOI_data$AOI <- unlist(AOI_data$AOI)
    AOI_data$FC <- unlist(AOI_data$FC)
    AOI_data$avg_FD <- unlist(AOI_data$avg_FD)
    AOI_data$sd <- unlist(AOI_data$sd)
    AOI_data$se <- unlist(AOI_data$se)
    AOI_data$group <- unlist(AOI_data$group)
    AOI_data$total_FD <- unlist(AOI_data$total_FD)
   
    for(i in 1:length(stimuli))
    {
        current_stimuli <- AOI_data[AOI_data$stimuli == stimuli[i], ] 
        
        stimuli_plot <- ggplot(current_stimuli, aes(x = factor(AOI), y = avg_FD, fill = factor(group))) +
                               geom_bar(stat = "identity", position = position_dodge(0.9)) +
                               theme(text = element_text(size = 15), 
                                     axis.text = element_text(colour = "black", size = 15),
                                     axis.text.x = element_text(angle = 70, hjust = 1, colour = "black"), 
                                     axis.text.y = element_text(colour = "black")) +
                                geom_errorbar(aes(ymin = avg_FD - se, ymax = avg_FD + se), position = position_dodge(0.9), width = 0.25) +
                                labs(title = stimuli[i], x = "ECG lead", y = "Mean fixation duration (ms)") +
                                scale_x_discrete(labels = swapLeadNames())
        
        print(stimuli_plot)
        
    }
    
    #utils::View(AOI_data)
    #aoi_plot <- ggplot(AOI_data, aes(x = facror))   
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
    common_source_data <- paste0(getwd(), "/Final-PHD-analysis/experiment_history.R")
    source(common_source_data) 
    
    # load data files
    eye_tracking <- loadEyeTrackingData()
    answers <- loadECGAnswers()
    survey <- loadSurveyData()
    AOI_locations <- loadAOILocations() 
    
    # fixation duration
    fixationDurationPerStimuli(eye_tracking, AOI_locations)
    #fixationDurationPerLead(eye_tracking, AOI_locations)
    
    #resultsSurvey(survey)
    #print(pairs(survey))
  
    
    # Examine accuracy, overall and by ECG
    #accuracyPerStimulus(answers)
    #accuracyPerParticipant(answers)
    #accuracyPerStimuli(answers)
    
    # Survey analysis
    surveyResultsBoxplots(survey)
    
    stop()
    
    # get scanpath sequence data for both conditions
    sequence_data_1 <-  filterSequenceData(eye_tracking[["cond_1"]])
    sequence_data_2 <-  filterSequenceData(eye_tracking[["cond_2"]])
    
    combined_data <- combineAccuracyData(rbind(sequence_data_1, sequence_data_2), answers)
    

}

# run main
main()
