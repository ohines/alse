require(dplyr)
# Based on Wallace et al code, with two changes to match (n=1732)
# Modifications made by David Whitney:
# 1. DW commented out a block that imputed Enzyme/Amiodarone..Cordarone
# 2. DW commented out filter on Warfarin Dose < 6
# The above changes yield a match on descriptive statistics reported in e.g.
# Schulz and Moodie (2021), Wallace et al (2018), Chen et al (2016)

# For additional files see data maintained at https://www.pharmgkb.org/downloads
# For data licence see https://www.pharmgkb.org/page/dataUsagePolicy
iwpc_url <- "https://api.pharmgkb.org/v1/download/submission/553247439"

data_directory <- "Data/"
raw_data_path <- paste0(data_directory, "iwpc_raw.xls")
output_path <- paste0(data_directory, "iwpc_clean.parquet")

download.file(iwpc_url, raw_data_path)

col_types <- rep("guess", times = 68)
col_types[37] <- "text"
col_types[58] <- "text"
col_types[68] <- "skip"

dose <- readxl::read_xls(
    raw_data_path,
    sheet = "Subject Data",
    na = "NA",
    col_types = col_types,
    .name_repair = "universal"
)

dose_filtered <- filter(
    dose,
    Subject.Reached.Stable.Dose.of.Warfarin == 1,
    !is.na(INR.on.Reported.Therapeutic.Dose.of.Warfarin),
    Race..OMB. != "Unknown",
    !is.na(Age),
    !is.na(Weight..kg.),
    !is.na(Height..cm.),
    !is.na(Gender),
    !is.na(CYP2C9.consensus)
)

#### genetics vs clinical ####
##### imputation based on page 13 of supplement of NEJM 2009 #########
#### IMPUTE rs9923231

dose_tmp3 <- dose_filtered
index9 <- which(is.na(dose_tmp3$VKORC1..1639.consensus))

dose_tmp3$VKORC1..1639.consensus[intersect(which(dose_tmp3$VKORC1.1173.consensus == "C/C"), index9)] <- "G/G"
dose_tmp3$VKORC1..1639.consensus[intersect(which(dose_tmp3$VKORC1.1173.consensus == "T/T"), index9)] <- "A/A"
dose_tmp3$VKORC1..1639.consensus[intersect(which(dose_tmp3$VKORC1.1173.consensus == "C/T"), index9)] <- "A/G"

index9 <- intersect(which(is.na(dose_tmp3$VKORC1..1639.consensus)), which(dose_tmp3$Race..OMB. == "White"))
dose_tmp3$VKORC1..1639.consensus[intersect(which(dose_tmp3$VKORC1.1542.consensus == "G/G"), index9)] <- "G/G"
dose_tmp3$VKORC1..1639.consensus[intersect(which(dose_tmp3$VKORC1.1542.consensus == "C/C"), index9)] <- "A/A"
dose_tmp3$VKORC1..1639.consensus[intersect(which(dose_tmp3$VKORC1.1542.consensus == "C/G"), index9)] <- "A/G"

dose_tmp3$VKORC1..1639.consensus[intersect(which(dose_tmp3$VKORC1.2255.consensus == "C/C"), index9)] <- "G/G"
dose_tmp3$VKORC1..1639.consensus[intersect(which(dose_tmp3$VKORC1.2255.consensus == "T/T"), index9)] <- "A/A"
dose_tmp3$VKORC1..1639.consensus[intersect(which(dose_tmp3$VKORC1.2255.consensus == "C/T"), index9)] <- "A/G"

dose_tmp4 <- dose_tmp3 %>%
    filter(!is.na(VKORC1..1639.consensus)) %>%
    mutate(
        Enzyme = case_when(
            ## one type is taken then enzyme is one, none of taken is zero, otherwise missing.
            (Rifampin.or.Rifampicin == 1) | (Carbamazepine..Tegretol. == 1) | (Phenytoin..Dilantin. == 1) ~ 1,
            (Rifampin.or.Rifampicin == 0) & (Carbamazepine..Tegretol. == 0) & (Phenytoin..Dilantin. == 0) ~ 0,
            .default = NA
        )
    )

### DW commented out -- it resulted in a different sample from comparison papers
## Updated version: We believe one would not use two type of drug at the same time.
# index10 = intersect(which(is.na(dose_tmp4$Enzyme)),which(dose_tmp4$Amiodarone..Cordarone. == 1))
# index11 = intersect(which(is.na(dose_tmp4$Amiodarone..Cordarone.)),which(dose_tmp4$Enzyme == 1))
# dose_tmp4$Enzyme[index10] = 0
# dose_tmp4$Amiodarone..Cordarone.[index11] = 0
### end DW comment out

df <- filter(
    dose_tmp4,
    !is.na(Amiodarone..Cordarone.),
    !is.na(Enzyme),
    Therapeutic.Dose.of.Warfarin <= 95
    # Therapeutic.Dose.of.Warfarin >= 6
)
# DW removed filter for Dose >= 6 as 5.81 is reported in ODTR comparison papers

### Final dataset construction
traindata <- dplyr::tibble(
    Age = as.numeric(as.factor(df$Age)) - 1,
    Weight = df$Weight..kg.,
    Height = df$Height..cm.,
    Enzyme = df$Enzyme,
    Amiodarone = df$Amiodarone..Cordarone.,
    Gender = as.numeric(as.factor(df$Gender)) - 1,
    Black = as.numeric(df$Race..OMB. == "Black or African American"),
    Asian = as.numeric(df$Race..OMB. == "Asian"),
    VKORC1.AG = as.numeric(df$VKORC1..1639.consensus == "A/G"),
    VKORC1.AA = as.numeric(df$VKORC1..1639.consensus == "A/A"),
    CYP2C9.12 = as.numeric(df$CYP2C9.consensus == "*1/*2"),
    CYP2C9.13 = as.numeric(df$CYP2C9.consensus == "*1/*3"),
    CYP2C9.other = as.numeric(
        !df$CYP2C9.consensus %in% c("*1/*1", "*1/*2", "*1/*3")
    ),
    Dose = df$Therapeutic.Dose.of.Warfarin,
    INR = df$INR.on.Reported.Therapeutic.Dose.of.Warfarin
)

arrow::write_parquet(traindata, output_path)
