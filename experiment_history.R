##################################################################################
# NAME:         experiment_history.R
# AUTHOUR:      Alan Davies
# DATE:         24/08/2017 
# INSTITUTION:  Interaction Analysis and Modelling Lab (IAM), University of Manchester
# DESCRIPTION:  Metadata for experiement (patient history).
#               
##################################################################################

#---------------------------------------------------------------------------------
# FUNCTION:     getStimuliList()
# INPUT:        void
# OUTPUT:       void
# DESCRIPTION:  Returns list of all the stimuli in the experiment
#               
#---------------------------------------------------------------------------------
getStimuliList <- function()
{
    study_stimuli <- c("anterolateralSTEMI",
                       "LBBB", 
                       "lateralSTEMI", 
                       "AF", 
                       "RBBB", 
                       "inferiorSTEMIandAF", 
                       "anteriorSTEMI",
                       "highlateralSTEMI", 
                       "inferolateralSTEMI") 
    
    return(study_stimuli)
}

#---------------------------------------------------------------------------------
# FUNCTION:     getExperimentSetupData()
# INPUT:        String
# OUTPUT:       list
# DESCRIPTION:  Returns metradata about the experimental stimuli
#               
#---------------------------------------------------------------------------------
getExperimentSetupData <- function(stimuli_name)
{
    stimuli <- list()
    stimuli[["anterolateralSTEMI"]] <- list(ref_name = "anterolateral STEMI", file_name = "anterolateralSTEMI.png", label = "Anterolateral STEMI", leads = 15)
    stimuli[["LBBB"]] <- list(ref_name = "LBBB", file_name = "LBBB.png", label = "Left Bundle Branch Block (LBBB)", leads = 15)
    stimuli[["lateralSTEMI"]] <- list(ref_name = "lateral STEMI", file_name = "lateralSTEMI.png", label = "Lateral STEMI", leads = 15)
    stimuli[["AF"]] <- list(ref_name = "AF", file_name = "AF.png", label = "Atrial fibrillation", leads = 15)
    stimuli[["RBBB"]] <- list(ref_name = "RBBB", file_name = "RBBB.png", label = "Right Bundle Branch Block (RBBB)", leads = 15)
    stimuli[["inferiorSTEMIandAF"]] <- list(ref_name = "inferior STEMI and AF", file_name = "inferiorSTEMIandAF.png", label = "Inferior STEMI and AF", leads = 15)
    stimuli[["anteriorSTEMI"]] <- list(ref_name = "anterior STEMI", file_name = "anteriorSTEMI.png", label = "Anterior STEMI", leads = 13)
    stimuli[["highlateralSTEMI"]] <- list(ref_name = "high lateral STEMI", file_name = "highlateralSTEMI.png", label = "High lateral STEMI", leads = 15)
    stimuli[["inferolateralSTEMI"]] <- list(ref_name = "inferolateral STEMI", file_name = "inferolateralSTEMI.png", label = "Inferolateral STEMI", leads = 15)

    return(stimuli[[stimuli_name]])   
}

#---------------------------------------------------------------------------------
# FUNCTION:     getAnswerNames()
# INPUT:        void
# OUTPUT:       vector
# DESCRIPTION:  Swap answer label for more meaninful ECG name
#               using a substitution vector 
#---------------------------------------------------------------------------------
getAnswerNames <- function()
{
    translation_vector <- c("LBBB", 
                            "Lateral STEMI",
                            "Atrial fibrillation",
                            "RBBB",
                            "Inferior STEMI and AF",
                            "Anterior STEMI",
                            "High lateral STEMI",
                            "Inferolateral STEMI",
                            "Anterolateral STEMI")
    
    return(translation_vector)
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
    translation_vector <- c("LBBB.png" = "LBBB", 
                            "lateralSTEMI.png" = "Lateral STEMI",
                            "AF.png" = "AF",
                            "RBBB.png" = "RBBB",
                            "inferiorSTEMIandAF.png" = "Inferior STEMI and AF",
                            "anteriorSTEMI.png" = "Anterior STEMI",
                            "highlateralSTEMI.png" = "High lateral STEMI",
                            "inferolateralSTEMI.png" = "Inferolateral STEMI",
                            "anterolateralSTEMI.png" = "Anterolateral STEMI")
    
    return(translation_vector)
}

#---------------------------------------------------------------------------------
# FUNCTION:     getLeadNames()
# INPUT:        String
# OUTPUT:       Stirng
# DESCRIPTION:  Turn the AOI letter into the equivalent lead name
#               
#---------------------------------------------------------------------------------
getLeadNames <- function(lead_id)
{
    switch(lead_id,
        "A" = "I", "B" = "II", "C" = "III", "D" = "aVR", "E" = "aVL",
        "F" = "aVF", "G" = "V1", "H" = "V2", "I" = "V3", "J" = "V4",
        "K" = "V5", "L" = "V6", "M" = "RS1", "N" = "RS2", "O" = "RS3"
    )
}

#---------------------------------------------------------------------------------
# FUNCTION:     swapLeadNames()
# INPUT:        void
# OUTPUT:       Stirng
# DESCRIPTION:  Turn the AOI letter into the equivalent lead name
#               
#---------------------------------------------------------------------------------
swapLeadNames <- function()
{
    translation_vector <- c("A" = "I", "B" = "II", "C" = "III", "D" = "aVR", "E" = "aVL",
                            "F" = "aVF", "G" = "V1", "H" = "V2", "I" = "V3", "J" = "V4",
                            "K" = "V5", "L" = "V6", "M" = "RS1", "N" = "RS2", "O" = "RS3")
    
    return(translation_vector)
}

#---------------------------------------------------------------------------------
# FUNCTION:     displayExperimentComparison()
# INPUT:        int
# OUTPUT:       Stirng
# DESCRIPTION:  Display which things are being compared to which other things
#               for multiple group analysis (C = correct / I = incorrect)
#
#               Saw HPC first?
#               T     |      F
#              -------|------
#               C     |     C
#               I     |     I
#               C     |     I
#               -------------
#               C/I   |
#                     |   C/I
#---------------------------------------------------------------------------------
displayExperimentComparison <- function(study_comparisons)
{
    if(study_comparisons == 1){
        cat("\n\n Saw HPC first: TRUE [correct] vs FALSE [correct]\n\n")
    } 
    else if(study_comparisons == 2){
        cat("\n\n Saw HPC first: TRUE [incorrect] vs FALSE [incorrect]\n\n")
    }
    else if(study_comparisons == 3){
        cat("\n\n Saw HPC first: TRUE [correct] vs FALSE [incorrect]\n\n")
    }
    else if(study_comparisons == 4){
        cat("\n\n Saw HPC first: TRUE [correct] vs [incorrect]\n\n")
    }
    else if(study_comparisons == 5){
        cat("\n\n Saw HPC first: FALSE [correct] vs [incorrect]\n\n")
    }
    else if(study_comparisons == 6){
        cat("\n\n All history [correct] vs [incorrect]\n\n")
    }
}

#---------------------------------------------------------------------------------
# FUNCTION:     getSelectedExperimentalGroups()
# INPUT:        data.frame, data.frame, int
# OUTPUT:       list
# DESCRIPTION:  Returns the correct 2 groups data in a list based on the
#               required experimental comparision
#---------------------------------------------------------------------------------
getSelectedExperimentalGroups <- function(study_comparisons, ...)
{
    dots <- list(...)
    group_data <- list()
    
    if(length(dots) == 2){
        if(study_comparisons == 1){
            group_data[[1]] <- dots[[1]][dots[[1]]$Accuracy == 1, ]
            group_data[[2]] <- dots[[2]][dots[[2]]$Accuracy == 1, ]
        } 
        else if(study_comparisons == 2){
            group_data[[1]] <- dots[[1]][dots[[1]]$Accuracy == 0, ]
            group_data[[2]] <- dots[[2]][dots[[2]]$Accuracy == 0, ]
        }
        else if(study_comparisons == 3){
            group_data[[1]] <- dots[[1]][dots[[1]]$Accuracy == 1, ]
            group_data[[2]] <- dots[[2]][dots[[2]]$Accuracy == 0, ]
        }
        else if(study_comparisons == 4){
            tmp <- dots[[1]]
            group_data[[1]] <- dots[[1]][dots[[1]]$Accuracy == 1, ]
            group_data[[2]] <- tmp[tmp$Accuracy == 0, ]
        }
        else if(study_comparisons == 5){
            tmp <- dots[[2]]
            group_data[[1]] <- dots[[2]][dots[[2]]$Accuracy == 1, ]
            group_data[[2]] <- tmp[tmp$Accuracy == 0, ]
        }
    } else if(length(dots) == 1) {
        if(study_comparisons == 6){
            groups <- dots[[1]]
            group_data[[1]] <- groups[groups$Accuracy == 1, ]
            group_data[[2]] <- groups[groups$Accuracy == 0, ]
        }
    }
    return(group_data)
}


