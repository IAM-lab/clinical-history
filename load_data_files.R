
  # condition 1 (history)
  condition_1 <- read_csv(file = paste0(getwd(), "/Final-PHD-analysis/data/cond1.csv"), col_types = cols(
    .default = col_integer(),
    ParticipantName = col_character(),
    `[Q01]Value` = col_character(),
    MediaName = col_character(),
    SaccadeIndex = col_character(),
    GazeEventType = col_character(),
    SaccadicAmplitude = col_double()
  ))
  
  # condition 2 (no-history)
  condition_2 <- read_csv(file = paste0(getwd(), "/Final-PHD-analysis/data/cond2.csv"), col_types = cols(
    .default = col_integer(),
    ParticipantName = col_character(),
    `[Q01]Value` = col_character(),
    MediaName = col_character(),
    SaccadeIndex = col_character(),
    GazeEventType = col_character(),
    SaccadicAmplitude = col_double()
  ))
  
  # answers
  answers <- read_csv(file = paste0(getwd(), "/Final-PHD-analysis/data/gs/answers_coded.csv"), col_types = cols(
    Participant = col_character(),
    Condition = col_integer(),
    `Answer 1` = col_integer(),
    `Answer 2` = col_integer(),
    `Answer 3` = col_integer(),
    `Answer 4` = col_integer(),
    `Answer 5` = col_integer(),
    `Answer 6` = col_integer(),
    `Answer 7` = col_integer(),
    `Answer 8` = col_integer(),
    `Answer 9` = col_integer()
  ))
  
  # AOI positions and dimensions per stimuli
  AOI_locations <- read_csv(file = paste0(getwd(), "/Final-PHD-analysis/data/AOI.csv"), col_types = cols(
    ECG = col_character(),
    condition = col_integer(),
    AOI = col_character(),
    left = col_integer(),
    top = col_integer(),
    width = col_integer(),
    height = col_integer()
  ))


