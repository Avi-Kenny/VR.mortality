# Title: "Estimating under-five mortality from vital registration data"
# Author: Avi Kenny
# Date: 2020-10-16



################.
##### INFO #####
################.

# The "Setup" sections create and save the following objects:
#   obj_name: obj description
#   obj_name: obj description
#   obj_name: obj description

# List of data objects
robjs <- c("national_data", "national_data_trans", "region_data",
           "region_data_trans", "region_data_trans_svy",
           "region_data_trans_vr", "mat", "years_3", "geo2", "reg_names")



#################.
##### SETUP #####
#################.

# Set working directory
if (Sys.getenv("USERDOMAIN")=="AVI-KENNY-T460") {
  setwd("C:/Users/avike/OneDrive/Desktop/Avi/Biostats + Research/Research/Jon Wakefield/Project - VR/z.vitalreg/R")
} else {
  setwd("z.stepped.wedge/R")
}

# Load packages
{
  library(rgdal)
  library(spdep)
  library(SUMMER)
  library(gridExtra)
  library(INLA)
  library(survey)
  library(readr)
  library(dplyr)
  library(magrittr)
  library(ggplot2)
  library(haven)
  library(tibble)
  library(tidyr)
  library(foreign)
  library(sqldf)
  library(rstan)
  library(raster)
  library(parallel)
}

# Load functions
{
  # source("my_func.R")
}

# Set code blocks to run
{
  run_setup_data <- FALSE
  run_load_data <- TRUE
  run_visualizations <- TRUE
  run_analysis <- FALSE
}



##########################################.
##### SETUP: Import Ecuador GIS data #####
##########################################.

if (run_setup_data) {
  
  # Read in and modify GIS data
  geo <- readOGR(
    dsn = "../data/Quito workshop/gadm36_ECU_shp",
    layer = "gadm36_ECU_1",
    verbose = FALSE
  )
  geo <- geo[geo$NAME_1 != "GalÃ¡pagos", ]
  levels(geo$NAME_1)[levels(geo$NAME_1)=="CaÃ±ar"] <- "Canar"
  geo$NAME_1 <- factor(geo$NAME_1)
  geo2 <- fortify(geo, region = "NAME_1")
  
  # Generate adjacency matrix "mat"
  nb.r <- poly2nb(geo, queen=F, row.names = geo$NAME_1)
  mat <- nb2mat(nb.r, style="B",zero.policy=TRUE)
  mat <- as.matrix(mat[1:dim(mat)[1], 1:dim(mat)[1]])
  
  # Save region names
  reg_names = rownames(mat)
  
}



###################################################.
##### SETUP: Import Ecuador ENSANUT 2012 data #####
###################################################.

# Source: faculty.washington.edu/jonno/UNICEF-WORKSHOPS/ensanut2012fbh.csv
# !!!!! Quality-check the entire data importing process

if (run_setup_data) {
  fbh <- read_csv("../data/Quito workshop/ensanut2012fbh.csv")
  fbh %<>% mutate(
    province = ifelse(province=="Cañar", "Canar", province)
  )
}



#####################################################.
##### SETUP: Import Ecuador VR data (1993-2010) #####
#####################################################.

# Source (1): www.ecuadorencifras.gob.ec/nacimientos-bases-de-datos/
# Source (2): www.ecuadorencifras.gob.ec/defunciones-generales-y-fetales-bases-de-datos/
# Creates temporary variables called "births_2000", "deaths_2000", etc.
# Note: defunciones_2009 has been modified (added header row)

if (run_setup_data) {
  
  for (i in 1993:2010) {
    
    # Read births datasets
    assign(
      paste0("births_",i),
      read_tsv(paste0("../data/Ecuador - VR (births)/nacimientos_",i,".dat"))
    )
    
    # Read deaths datasets
    assign(
      paste0("deaths_",i),
      read_tsv(paste0("../data/Ecuador - VR (deaths)/defunciones_",i,".dat"))
    )
    
  }
  
  # Parse VR data into a single data frame
  df_VR <- data.frame(
    "region" = character(),
    "year" = integer(),
    "births" = integer(),
    "deaths" = integer(),
    stringsAsFactors=FALSE
  )
  
  # Set up regions_vr list
  # Note this is ONLY used for importing data; use `reg_names` for IDs
  regions_vr = list(
    "01" = "Azuay",
    "02" = "Bolivar",
    "03" = "Canar",
    "04" = "Carchi",
    "05" = "Cotopaxi",
    "06" = "Chimborazo",
    "07" = "El Oro",
    "08" = "Esmeraldas",
    "09" = "Guayas",
    "10" = "Imbabura",
    "11" = "Loja",
    "12" = "Los Rios",
    "13" = "Manabi",
    "14" = "Morona Santiago",
    "15" = "Napo",
    "16" = "Pastaza",
    "17" = "Pichincha",
    "18" = "Tungurahua",
    "19" = "Zamora Chinchipe",
    "21" = "Sucumbios",
    "22" = "Orellana",
    "23" = "Santo Domingo de los Tsachilas",
    "24" = "Santa Elena"
  )
  
  # Function: Update data frame with VR data
  update_df_VR = function(year, df_VR, df_births, df_deaths) {
    
    # Different death codes for different VR years
    death_code = ifelse(year<=2007, "CODI", "COD_EDA")
    
    # Filter deaths table to only include under-five deaths
    if (year <= 2007) {
      df_deaths %<>% filter(
        CODI != 9 &
          (CODI<=3 | (CODI==4 & EDAD<5))
      )
    }
    if (year >= 2008) {
      df_deaths %<>% filter(
        COD_EDA != 9 &
          (COD_EDA<=3 | (COD_EDA==4 & EDAD<5))
      )
    }
    
    # Add regional rows
    for (j in 1:length(regions_vr)) {
      
      # Add row
      i = nrow(df_VR) + 1
      region_code = names(regions_vr[j])[1]
      region_name = regions_vr[[j]]
      births = df_births %>%
        filter(as.numeric(PRO_RES) == as.numeric(region_code)) %>%
        nrow()
      deaths = df_deaths %>%
        filter(as.numeric(PROVF) == as.numeric(region_code)) %>%
        nrow()
      df_VR[i,"region"] = region_name
      df_VR[i,"year"] = year
      df_VR[i,"births"] = births
      df_VR[i,"deaths"] = deaths
      
    }
    
    # Return updated table
    return (df_VR)
    
  }
  
  # Add VR data to df_VR
  for (i in 1993:2010) {
    
    df_VR = update_df_VR(
      i,
      df_VR,
      eval(parse(text=paste0("births_",i))),
      eval(parse(text=paste0("deaths_",i)))
    )
    
  }
  
  # Aggregate VR regional deaths to three-year blocks
  df_VR %<>% mutate(
    year = case_when(
      between(year, 1993, 1995) ~ "93-95",
      between(year, 1996, 1998) ~ "96-98",
      between(year, 1999, 2001) ~ "99-01",
      between(year, 2002, 2004) ~ "02-04",
      between(year, 2005, 2007) ~ "05-07",
      between(year, 2008, 2010) ~ "08-10",
      TRUE ~ "error"
    )
  ) %>% group_by(region, year) %>% summarize(
    births = sum(births),
    deaths = sum(deaths)
  )
  
  # Ungroup
  df_VR %<>% ungroup()
  
  # Clean up objects
  for (i in 1993:2010) {
    rm(list=paste0("births_",i))
    rm(list=paste0("deaths_",i))
  }
  rm(regions_vr)
  
}



########################################################.
##### SETUP: Calculate mortality rates for ENSANUT #####
########################################################.

# From Wakefield Quito workshop

if (run_setup_data) {
  
  # Calculate variables
  fbh %<>% mutate(
    u5_death = ifelse(fbh$alive=="no" & fbh$B7<60,1,0),
    birth = 1,
    dob_years = case_when(
      between(B3, 1117, 1152) ~ "93-95",
      between(B3, 1153, 1188) ~ "96-98",
      between(B3, 1189, 1224) ~ "99-01",
      between(B3, 1225, 1260) ~ "02-04",
      between(B3, 1261, 1296) ~ "05-07",
      between(B3, 1297, 1332) ~ "08-10",
      TRUE ~ "other"
    ),
    dod_years = case_when(
      between(B3+B7, 1117, 1152) ~ "93-95",
      between(B3+B7, 1153, 1188) ~ "96-98",
      between(B3+B7, 1189, 1224) ~ "99-01",
      between(B3+B7, 1225, 1260) ~ "02-04",
      between(B3+B7, 1261, 1296) ~ "05-07",
      between(B3+B7, 1297, 1332) ~ "08-10",
      TRUE ~ "other"
    )
  )
  
  # Set survey design
  # !!!!! double-check this
  design <- svydesign(
    ids = ~V021,
    weights = ~V005,
    strata = ~V024,
    data = fbh
  )
  
  # Calculate direct counts of births and deaths
  direct_births <- svyby(
    formula = ~birth,
    by = ~province+dob_years,
    design = design,
    FUN = svytotal
  )
  direct_deaths <- svyby(
    formula = ~u5_death,
    by = ~province+dod_years,
    design = design,
    FUN = svytotal
  )
  
  # Code below adapted from Wakefield Quito workshop code
  
  # Region-level estimates
  births_3_years <- getBirths(
    data = fbh,
    surveyyear = 2012,
    variables = c("idhog", "province", "V005", "V021", "V024", "V025"),
    strata = "V024",
    dob = "B3",
    alive = "alive",
    age = "B7",
    date.interview = "V008",
    year.cut = seq(1993, 2011, by = 3)
  )
  years_3 <- levels(births_3_years$time)
  df_survey <- getDirect(
    births = births_3_years,
    years = years_3,
    regionVar = "province",
    timeVar = "time",
    clusterVar = "~V021+idhog",
    ageVar = "age",
    weightsVar = "V005",
    geo.recode = NULL,
    national.only = FALSE
  )
  
  # National estimates
  df_survey_national <- getDirect(
    births = births_3_years,
    years = years_3,
    regionVar = "province",
    timeVar = "time",
    clusterVar = "~V021+idhog",
    ageVar = "age",
    weightsVar = "V005",
    geo.recode = NULL,
    national.only = TRUE
  )
  
}



#################################################.
##### SETUP: Create data frame for analysis #####
#################################################.

# Uses `df_survey` object and `df_VR` objects generated above
# !!!!! Import standard errors

if (run_setup_data) {
  
  # Set up variables
  years = df_VR$year %>% unique()
  
  # Set up data frames
  region_data <- data.frame(
    "id_reg" = double(),
    "region" = character(),
    "years" = integer(),
    "n_births_vr" = integer(),
    "n_deaths_vr" = integer(),
    "n_births_svy" = integer(),
    "n_deaths_svy" = integer(),
    "u5mr_svy" = double(),
    stringsAsFactors=FALSE
  )
  national_data <- df_survey_national %>%
    subset(select=c(years,mean)) %>%
    rename("u5mr_svy"=`mean`) %>%
    mutate("u5mr_vr"=NA)
  
  # Fill national data frame
  for (y in 1:length(years)) {
    yr <- years[y]
    d <- df_VR %>% filter(year==yr)
    u5mr <- round(sum(d$deaths)/sum(d$births),5)
    national_data %<>% mutate(
      u5mr_vr = ifelse(years==yr, u5mr, u5mr_vr)
    )
  }
  
  # Fill regional data frame
  for (r in 1:length(reg_names)) {
    for (y in 1:length(years)) {
      
      # Create alias for years
      yrs <- years
      
      # Get n_births_vr
      n_births_vr <- df_VR %>%
        filter(region==reg_names[[r]] &
                 years==yrs[y]) %>%
        subset(select=births) %>% as.numeric()
      
      # Get n_deaths_vr
      n_deaths_vr <- df_VR %>%
        filter(region==reg_names[[r]] &
                 years==yrs[y]) %>%
        subset(select=deaths) %>% as.numeric()
      
      # Get n_births_svy
      n_births_svy <- direct_births %>%
        filter(province==reg_names[[r]] & 
                 dob_years==yrs[y]) %>%
        subset(select=birth) %>% as.numeric() %>%
        round()
      
      # Get n_deaths_svy
      n_deaths_svy <- direct_deaths %>%
        filter(province==reg_names[[r]] &
                 dod_years==yrs[y]) %>%
        subset(select=u5_death) %>% as.numeric() %>%
        round()
      
      # Get u5mr_svy
      u5mr_svy <- df_survey %>%
        filter(region==reg_names[[r]] &
                 years==yrs[y]) %>%
        subset(select=mean) %>% as.numeric()
      
      # Add row to data frame
      region_data[nrow(region_data)+1,] <- list(
        r, reg_names[r],
        years[y], n_births_vr, n_deaths_vr,
        n_births_svy, n_deaths_svy, u5mr_svy
      )
      
    }
  }
  
  # Remove one row with no survey deaths (Tungurahua 08-10)
  region_data %<>% filter(!is.na(u5mr_svy))
  
  # Generate calculated variables
  # !!!!! Use more complicated formula to calculate U5MR from VR
  region_data %<>% mutate(
    u5mr_vr = round(n_deaths_vr/n_births_vr,5),
    time = case_when(
      years=="93-95" ~ 1,
      years=="96-98" ~ 2,
      years=="99-01" ~ 3,
      years=="02-04" ~ 4,
      years=="05-07" ~ 5,
      years=="08-10" ~ 6,
      TRUE ~ 999
    ),
    log_mbias = log(u5mr_vr/u5mr_svy)
  )
  national_data %<>% mutate(
    time = case_when(
      years=="93-95" ~ 1,
      years=="96-98" ~ 2,
      years=="99-01" ~ 3,
      years=="02-04" ~ 4,
      years=="05-07" ~ 5,
      years=="08-10" ~ 6,
      TRUE ~ 999
    ),
    log_mbias = log(u5mr_vr/u5mr_svy)
  )
  
  # Create transformed datasets
  # Source: 0=VR, 1=survey
  region_data_trans <- sqldf(
    "SELECT id_reg, region, years, time, 1 AS source,
   n_births_svy AS n_births, n_deaths_svy AS n_deaths,
   u5mr_svy AS u5mr, log_mbias FROM region_data
   
   UNION SELECT id_reg, region, years, time, 0,
   n_births_vr, n_deaths_vr, u5mr_vr, NULL FROM region_data
   
   ORDER BY region, years, source"
  )
  national_data_trans <- sqldf(
    "SELECT years, time, 'Survey' AS source, u5mr_svy AS u5mr,
   log_mbias FROM national_data
   
   UNION SELECT years, time, 'VR', u5mr_vr,
   NULL FROM national_data
   ORDER BY source, time"
  )
  
  # Print data structure
  head(region_data)
  head(region_data_trans)
  head(national_data)
  head(national_data_trans)
  
  # Create copies of id variables
  region_data_trans %<>% mutate(
    id_reg_copy = id_reg,
    time_copy = time
  )
  
  # Create filtered data sources
  region_data_trans_svy <- region_data_trans %>%
    filter(source==0) %>%
    rowid_to_column("id_row")
  
  region_data_trans_vr <- region_data_trans %>%
    filter(source==1) %>%
    rowid_to_column("id_row")
  
}



#################################################.
##### SETUP: Save objects for quick loading #####
#################################################.

if (run_setup_data) {
  
  for (i in 1:length(robjs)) {
    filename <- paste0("../robj/", robjs[i], ".robj")
    saveRDS(eval(as.name(robjs[i])), filename)
  }
  
}



################################################.
##### SETUP: Load previously saved objects #####
################################################.

if (run_load_data) {
  
  for (i in 1:length(robjs)) {
    filename <- paste0("../robj/", robjs[i], ".robj")
    assign(robjs[i], readRDS(filename))
  }
  
}



###########################################.
##### VISUAL: U5MR visualizations (1) #####
###########################################.

if (run_visualizations) {
  
  # Get survey 80% CIs
  svy_ci_low <- expit(
    logit(df_survey_national$mean) +
      qnorm(0.1) * sqrt(df_survey_national$var.est)
  )
  svy_ci_high <- expit(
    logit(df_survey_national$mean) +
      qnorm(0.9) * sqrt(df_survey_national$var.est)
  )
  
  # Get VR counts
  vr_births <- c()
  vr_deaths <- c()
  for (y in 1:length(years)) {
    yr <- years[y]
    d <- df_VR %>% filter(year==yr)
    vr_births <- c(vr_births, sum(d$births))
    vr_deaths <- c(vr_deaths, sum(d$deaths))
  }
  vr_births <- c(vr_births[4:6],vr_births[1:3])
  vr_deaths <- c(vr_deaths[4:6],vr_deaths[1:3])
  
  # Calculate VR CIs
  vr_ci_low <- c()
  vr_ci_high <- c()
  for (i in 1:6) {
    B <- rpois(n=1000, lambda=vr_births[i])
    D <- rbinom(n=1000, size=B, prob=(vr_deaths[i]/vr_births[i]))
    vr_ci_low <- c(vr_ci_low, quantile(D/B, 0.1))
    vr_ci_high <- c(vr_ci_high, quantile(D/B, 0.9))
  }
  
  # Temp "Measurement error" adjustment
  vr_ci_ranges <- 0.5*(vr_ci_high-vr_ci_low)
  vr_ci_low <- vr_ci_low - 5*vr_ci_ranges
  vr_ci_high <- vr_ci_high + 5*vr_ci_ranges
  
  # Line graph: U5MR vs. time, by source
  ggplot(
    data = cbind(
      national_data_trans,
      CI_low = c(svy_ci_low, vr_ci_low),
      CI_high = c(svy_ci_high, vr_ci_high)
    ),
    aes(x=time, y=(u5mr*1000),
        group = source,
        color = as.factor(source))
  ) +
    geom_ribbon(
      aes(
        ymin = (CI_low*1000),
        ymax = (CI_high*1000),
        fill = "80% CI",
        color = "80% CI",
      ),
      color = NA,
      alpha = 0.2
    ) +
    geom_line(size=1.5) +
    scale_fill_manual("",values="grey12") +
    scale_x_discrete(
      limits = 1:6,
      labels = years_3
    ) +
    labs(
      # title = "U5MR vs. time, by source",
      x = "Time period (years)",
      y = "U5MR",
      color = "Data source"
    ) +
    theme(
      text = element_text(size=20)
    )
  
}



###########################################.
##### VISUAL: U5MR visualizations (2) #####
###########################################.

if (run_visualizations) {
  
  # Map: U5MR (VR) by region, [year]
  for (i in 1:length(years_3)) {
    print(
      plot_map(
        df = region_data %>%
          filter(years==years_3[i]) %>%
          subset(select=c(id_reg,u5mr_vr)) %>%
          mutate(u5mr_vr=1000*u5mr_vr),
        geo2 = geo2,
        title = paste("U5MR (VR), years:",years_3[i]),
        limits = c(0,60)
      )
    )
  }
  
  # Map: U5MR (SVY) by region, [year]
  for (i in 1:length(years_3)) {
    print(
      plot_map(
        df = region_data %>%
          filter(years==years_3[i]) %>%
          subset(select=c(id_reg,u5mr_svy)) %>%
          mutate(u5mr_svy=1000*u5mr_svy),
        geo2 = geo2,
        title = paste("U5MR (SVY), years:",years_3[i]),
        limits = c(0,60)
      )
    )
  }
  
  # Line graph: U5MR (VR) vs. time, by region (sample)
  ggplot(
    data = region_data %>%
      filter(region %in% reg_names[1:10]) %>%
      mutate(u5mr_vr=1000*u5mr_vr),
    aes(x=time, y=u5mr_vr, group=region,
        color=as.factor(region))
  ) +
    geom_line() +
    scale_x_discrete(limits=1:6, labels=years_3) +
    labs(
      color = "Region",
      x = "Time period",
      y = "U5MR",
      title = "U5MR (VR) vs. time, by region (sample)"
    )
  
  # Line graph: U5MR (SVY) vs. time, by region (sample)
  ggplot(
    data = region_data %>%
      filter(region %in% reg_names[1:10]) %>%
      mutate(u5mr_svy=1000*u5mr_svy),
    aes(x=time, y=u5mr_svy, group=region,
        color=as.factor(region))
  ) +
    geom_line() +
    scale_x_discrete(limits=1:6, labels=years_3) +
    labs(
      color = "Region",
      x = "Time period",
      y = "U5MR",
      title = "U5MR (SVY) vs. time, by region (sample)"
    )
  
}



#############################################.
##### VISUAL: log(MBIAS) visualizations #####
#############################################.

if (run_visualizations) {
  
  # Line graph: MBIAS vs. time
  ggplot(national_data, aes(x=time, y=exp(log_mbias))) + 
    geom_line() +
    scale_x_discrete(
      limits = 1:6,
      labels = years_3
    ) +
    labs(
      x = "Time period",
      y = "MBIAS",
      title = "MBIAS vs. time"
    )
  
  # Plot maps of log(MBIAS) for all six time periods separately
  for (i in 1:length(years_3)) {
    print(
      plot_map(
        df = region_data %>%
          filter(years==years_3[i]) %>%
          subset(select=c(region,log_mbias)),
        geo2 = geo2,
        title = paste("log(MBIAS), years:",years_3[i]),
        limits = c(-2.2, 2.2),
        sfg = "rwg"
      )
    )
  }
  
  # Line graph: log(MBIAS) vs. time, by region (sample)
  ggplot(
    data = region_data %>%
      filter(region %in% reg_names[1:10]),
    aes(x=time, y=log_mbias, group=region,
        color=as.factor(region))
  ) +
    geom_line() +
    scale_x_discrete(limits=1:6, labels=years_3) +
    labs(
      color = "Region",
      x = "Time period",
      y = "log(MBIAS)",
      title = "log(MBIAS) vs. time, by region (sample)"
    )
  
}



################################.
##### ANALYSIS: Run models #####
################################.

if (run_analysis) {
  
  # Model 4
  model_4_vr <- inla(
    n_deaths ~ 1 +
      f(id_reg, model="iid"),
    family = "binomial",
    Ntrials = region_data_trans_vr$n_births,
    data = region_data_trans_vr,
    control.family = list(link="logit"),
    control.compute=list(dic=TRUE)
  )
  model_4_svy <- inla(
    n_deaths ~ 1 +
      f(id_reg, model="iid"),
    family = "binomial",
    Ntrials = region_data_trans_svy$n_births,
    data = region_data_trans_svy,
    control.family = list(link="logit"),
    control.compute=list(dic=TRUE)
  )
  
  # Model 5
  model_5_vr <- inla(
    n_deaths ~ 1 +
      f(time, model="rw2"),
    family = "binomial",
    Ntrials = region_data_trans_vr$n_births,
    data = region_data_trans_vr,
    control.family = list(link="logit"),
    control.compute=list(dic=TRUE)
  )
  model_5_svy <- inla(
    n_deaths ~ 1 +
      f(time, model="rw2"),
    family = "binomial",
    Ntrials = region_data_trans_svy$n_births,
    data = region_data_trans_svy,
    control.family = list(link="logit"),
    control.compute=list(dic=TRUE)
  )
  
  # Model 6
  model_6_vr <- inla(
    n_deaths ~ 1 +
      f(id_reg, model="iid") +
      f(time, model="rw2"),
    family = "binomial",
    Ntrials = region_data_trans_vr$n_births,
    data = region_data_trans_vr,
    control.family = list(link="logit"),
    control.compute=list(dic=TRUE)
  )
  model_6_svy <- inla(
    n_deaths ~ 1 +
      f(id_reg, model="iid") +
      f(time, model="rw2"),
    family = "binomial",
    Ntrials = region_data_trans_svy$n_births,
    data = region_data_trans_svy,
    control.family = list(link="logit"),
    control.compute=list(dic=TRUE)
  )
  
  # Model 7
  model_7_vr <- inla(
    n_deaths ~ 1 +
      f(id_reg, model="iid") +
      f(id_reg_copy, model="besag", graph=mat) +
      f(time, model="rw2") +
      f(time_copy, model="iid"),
    family = "binomial",
    Ntrials = region_data_trans_vr$n_births,
    data = region_data_trans_vr,
    control.family = list(link="logit"),
    control.compute=list(dic=TRUE)
  )
  model_7_svy <- inla(
    n_deaths ~ 1 +
      f(id_reg, model="besag", graph=mat) +
      f(id_reg_copy, model="iid") +
      f(time, model="rw2") +
      f(time_copy, model="iid"),
    family = "binomial",
    Ntrials = region_data_trans_svy$n_births,
    data = region_data_trans_svy,
    control.family = list(link="logit"),
    control.compute=list(dic=TRUE)
  )
  
  # Models 7a/7b are the same as above but with the BYM2 parameterization
  #   (and with the DIC computed)
  model_7b_vr <- inla(
    n_deaths ~ 1 +
      f(id_reg, model="bym2", graph=mat) +
      f(time, model="rw2") +
      f(time_copy, model="iid"),
    family = "binomial",
    Ntrials = region_data_trans_vr$n_births,
    data = region_data_trans_vr,
    control.family = list(link="logit"),
    control.compute=list(dic=TRUE)
  )
  model_7b_svy <- inla(
    n_deaths ~ 1 +
      f(id_reg, model="bym2", graph=mat) +
      f(time, model="rw2") +
      f(time_copy, model="iid"),
    family = "binomial",
    Ntrials = region_data_trans_svy$n_births,
    data = region_data_trans_svy,
    control.family = list(link="logit"),
    control.compute=list(dic=TRUE)
  )
  
  # Copy vectors for use in models 8a-d
  id_reg_vr <- region_data_trans_vr$id_reg
  time_vr <- region_data_trans_vr$time
  id_reg_svy <- region_data_trans_svy$id_reg
  time_svy <- region_data_trans_svy$time
  
  # Model 8a
  model_8a_vr <- inla(
    n_deaths ~ 1 +
      f(id_reg, model="bym2", graph=mat) +
      f(time, model="rw2") +
      f(time_copy, model="iid") +
      f(id_row, model="iid"),
    family = "binomial",
    Ntrials = region_data_trans_vr$n_births,
    data = region_data_trans_vr,
    control.family = list(link="logit"),
    control.compute=list(dic=TRUE)
  )
  model_8a_svy <- inla(
    n_deaths ~ 1 +
      f(id_reg, model="bym2", graph=mat) +
      f(time, model="rw2") +
      f(time_copy, model="iid") +
      f(id_row, model="iid"),
    family = "binomial",
    Ntrials = region_data_trans_svy$n_births,
    data = region_data_trans_svy,
    control.family = list(link="logit"),
    control.compute=list(dic=TRUE)
  )
  
  # Model 8b
  model_8b_vr <- inla(
    n_deaths ~ 1 +
      f(id_reg, model="bym2", graph=mat) +
      f(time, model="rw2") +
      f(time_copy, model="iid") +
      f(id_reg_vr, model="iid", group=time_vr,
        control.group=list(model="rw2")),
    family = "binomial",
    Ntrials = region_data_trans_vr$n_births,
    data = region_data_trans_vr,
    control.family = list(link="logit"),
    control.compute=list(dic=TRUE)
  )
  model_8b_svy <- inla(
    n_deaths ~ 1 +
      f(id_reg, model="bym2", graph=mat) +
      f(time, model="rw2") +
      f(time_copy, model="iid") +
      f(id_reg_svy, model="iid", group=time_svy,
        control.group=list(model="rw2")),
    family = "binomial",
    Ntrials = region_data_trans_svy$n_births,
    data = region_data_trans_svy,
    control.family = list(link="logit"),
    control.compute=list(dic=TRUE)
  )
  
  # Model 8c
  model_8c_vr <- inla(
    n_deaths ~ 1 +
      f(id_reg, model="bym2", graph=mat) +
      f(time, model="rw2") +
      f(time_copy, model="iid") +
      f(time_vr, model="iid", group=id_reg_vr,
        control.group=list(model="besag", graph=mat)),
    family = "binomial",
    Ntrials = region_data_trans_vr$n_births,
    data = region_data_trans_vr,
    control.family = list(link="logit"),
    control.compute=list(dic=TRUE)
  )
  model_8c_svy <- inla(
    n_deaths ~ 1 +
      f(id_reg, model="bym2", graph=mat) +
      f(time, model="rw2") +
      f(time_copy, model="iid") +
      f(time_svy, model="iid", group=id_reg_svy,
        control.group=list(model="besag", graph=mat)),
    family = "binomial",
    Ntrials = region_data_trans_svy$n_births,
    data = region_data_trans_svy,
    control.family = list(link="logit"),
    control.compute=list(dic=TRUE)
  )
  
  # Model 8d
  model_8d_vr <- inla(
    n_deaths ~ 1 +
      f(id_reg, model="bym2", graph=mat) +
      f(time, model="rw2") +
      f(time_copy, model="iid") +
      f(id_reg_vr, model="besag", graph=mat,
        group=time_vr, control.group=list(model="rw2")),
    family = "binomial",
    Ntrials = region_data_trans_vr$n_births,
    data = region_data_trans_vr,
    control.family = list(link="logit"),
    control.compute=list(dic=TRUE)
  )
  model_8d_svy <- inla(
    n_deaths ~ 1 +
      f(id_reg, model="bym2", graph=mat) +
      f(time, model="rw2") +
      f(time_copy, model="iid") +
      f(id_reg_svy, model="besag", graph=mat,
        group=time_svy, control.group=list(model="rw2")),
    family = "binomial",
    Ntrials = region_data_trans_svy$n_births,
    data = region_data_trans_svy,
    control.family = list(link="logit"),
    control.compute=list(dic=TRUE)
  )
  
}



#####################################################.
##### VISUAL: Random effects visualizations (1) #####
#####################################################.

if (run_visualizations) {
  
  # Plot time effects
  for (m in c(5,6,7)) {
    for (src in c("svy","vr")) {
      
      # Get reference to model
      model <- eval(as.name(paste("model",m,src,sep="_")))
      
      # Get time random effects
      assign(
        paste("time", src, sep="_"),
        model$summary.random$time[,2]
      )
      
    }
    
    print(
      ggplot(data.frame(x=time_vr, y=time_svy),
             aes(x=x, y=y)) +
        geom_point(color="turquoise") +
        labs(title=paste("RW2 effects, model",m)) +
        labs(y='VR', x='Survey') +
        geom_abline(slope=1, intercept=0)
    )
    
  }
  
  # Plot IID effects
  for (m in c(4,6,7)) {
    for (src in c("svy","vr")) {
      
      # Get reference to model
      model <- eval(as.name(paste("model",m,src,sep="_")))
      
      # Get IID random effects
      assign(
        paste("iid", src, sep="_"),
        model$summary.random$id_reg[,2]
      )
      
    }
    
    print(
      ggplot(data.frame(x=iid_vr, y=iid_svy),
             aes(x=x, y=y)) +
        geom_point(color="turquoise") +
        labs(title=paste("IID effects, model",m)) +
        labs(y='VR', x='Survey') +
        geom_abline(slope=1, intercept=0)
    )
    
  }
  
}



#####################################################.
##### VISUAL: Random effects visualizations (2) #####
#####################################################.

if (run_visualizations) {
  
  # Map IID effects
  for (m in c(4,6,7)) {
    for (src in c("svy","vr")) {
      
      # Get reference to model
      model <- eval(as.name(paste("model",m,src,")",sep="_")))
      
      # Get IID random effects
      assign(
        paste("iid", src, sep="_"),
        model$summary.random$id_reg[,2]
      )
      
    }
    
    print(plot_map(
      data.frame(id=(1:23), value=iid_svy),
      geo2 = geo2,
      title = paste0("IID effects: SVY (model ", m, ""),
      limits = c(-0.5,0.5),
      sfg = "rwg"
    ))
    print(plot_map(
      data.frame(id=(1:23), value=iid_vr),
      geo2 = geo2,
      title = paste0("IID effects: VR (model ", m, ""),
      limits = c(-0.5,0.5),
      sfg = "rwg"
    ))
    
  }
  
}



#########################################################.
##### ANALYSIS: Run models (separate VR+SVY models) #####
#########################################################.

if (run_analysis) {
  
  results <- data.frame(
    "model" = integer(),
    "source" = character(),
    "intercept" = double(),
    "var_space" = double(),
    "phi_space" = double(),
    "var_time_rw2" = double(),
    "var_time_iid" = double(),
    "var_int" = double(),
    "dic" = double(),
    stringsAsFactors=FALSE
  )
  
  for (m in c("7b","8a","8b","8c","8d")) {
    for (src in c("svy","vr")) {
      
      # Get reference to model
      model <- eval(as.name(paste("model",m,src,sep="_")))
      summ <- summary(model)
      
      results[nrow(results)+1,] <- list(
        model = m,
        source = src,
        intercept = summ$fixed[1,1],
        var_space = 1/summ$hyperpar[1,1],
        phi_space = summ$hyperpar[2,1],
        var_time_rw2 = 1/summ$hyperpar[3,1],
        var_time_iid = 1/summ$hyperpar[4,1],
        var_int = ifelse(m=="7b", NA, 1/summ$hyperpar[5,1]),
        dic = summ$dic$dic
      )
      
    }
  }
  
  results %>% subset(
    select=c(model, source, intercept, var_space, phi_space,
             var_time_rw2, var_time_iid)
  ) %>% print()
  
  results %>% subset(
    select=c(model, source, var_int, dic)
  ) %>% print()
  
}



#########################################################.
##### ANALYSIS: Run models (combined VR+SVY models) #####
#########################################################.

if (run_analysis) {
  
  # Model #1 (INLA)
  model_1 <- inla(
    n_deaths ~ 1 + source,
    family = "binomial",
    Ntrials = region_data_trans$n_births,
    data = region_data_trans,
    control.family = list(link="logit"),
    control.predictor = list(compute=TRUE)
  )
  summary_1 <- summary(model_1)
  summary_1
  
  # Actual mean, by source (0,1)
  sum(region_data_trans_svy$n_deaths, na.rm=TRUE) / 
    sum(region_data_trans_svy$n_births)
  sum(region_data_trans_vr$n_deaths, na.rm=TRUE) / 
    sum(region_data_trans_vr$n_births)
  
  # Predicted mean, by source (0,1)
  beta_intercept <- summary_1$fixed[1,"mean"]
  beta_covariate <- summary_1$fixed[2,"mean"]
  expit(beta_intercept)
  expit(beta_intercept + beta_covariate)
  
  # Set up STAN
  options(mc.cores = detectCores()-1)
  rstan_options(auto_write = TRUE)
  
  # Model #1 (STAN)
  stan_code <- quote("
    data {
      int N;
      int id_reg[N];
      int time[N];
      row_vector[N] source;
      int n_births[N];
      int n_deaths[N];
      real u5mr[N];
    }
    parameters {
      real beta_0;
      real beta_1;
    }
    model {
      row_vector[N] p;
      p = inv_logit(beta_0 + beta_1 * source);
      beta_0 ~ normal(0, 100);
      beta_1 ~ normal(0, 100);
      n_deaths ~ binomial(n_births, p);
    }
  ")
  
  stan_data <- region_data_trans %>% 
    subset( select=-c(region,years,log_mbias) ) %>%
    filter(!is.na(n_deaths) & !is.na(u5mr)) %>%
    as.list()
  stan_data$N = length(stan_data$id_reg)
  stan_data$T = length(unique(stan_data$time))
  stan_data$I = length(unique(stan_data$id_reg))
  fit <- stan(
    model_code = stan_code,
    data = stan_data
  )
  print(fit)
  
  # Model #2 (INLA)
  model_2 <- inla(
    n_deaths ~ 1 +
      source +
      f(time, model="ar1") +
      f(id_reg, model="bym2", graph=mat),
    family = "binomial",
    Ntrials = region_data_trans$n_births,
    data = region_data_trans,
    control.family = list(link="logit"),
    control.predictor = list(compute=TRUE)
  )
  summary_2 <- summary(model_2)
  summary_2
  
}



###########################.
##### MISC: BYM model #####
###########################.

# !!!!! Why do these two not give the same estimates?
#         - look at precision of spatial component

if (run_misc) {
  
  model_1 <- inla(
    n_deaths ~ 1 +
      f(id_reg, model="bym", graph=mat),
    family = "binomial",
    Ntrials = region_data_trans_vr$n_births,
    data = region_data_trans_vr,
    control.family = list(link="logit"),
    control.predictor = list(compute=TRUE),
    control.compute=list(dic=TRUE)
  )
  
  model_2 <- inla(
    n_deaths ~ 1 +
      f(id_reg_copy, model="besag", graph=mat) +
      f(id_reg, model="iid"),
    family = "binomial",
    Ntrials = region_data_trans_vr$n_births,
    data = region_data_trans_vr,
    control.family = list(link="logit"),
    control.predictor = list(compute=TRUE),
    control.compute=list(dic=TRUE)
  )
  
  summary(model_1)
  summary(model_2)
  
}



###################.
##### ARCHIVE #####
###################.

if (FALSE) {
  
  # # Plot basic map of Ecuador with regions labeled
  # ggplot(data=geo2) +
  #   geom_polygon(aes(x=long, y=lat, group=group, fill=id),
  #                color="black") +
  #   theme_void() +
  #   coord_map()
  
  
  
  # # STAN model #3
  # stan_code <- quote("
  # // saved as stan_test.stan
  # data {
  #   int N;
  #   int T;
  #   int id_reg[249];
  #   int time[249];
  #   row_vector[249] source;
  #   int n_births[249];
  #   int n_deaths[249];
  #   real u5mr[249];
  # }
  # parameters {
  #   real beta_0;
  #   real beta_1;
  #   real<lower=0> sigma;
  #   real gamma;
  #   real<lower=0> rho;
  # }
  # model {
  #   row_vector[249] p;
  #   p = inv_logit(beta_0 + beta_1 * source + gamma);
  #   beta_0 ~ normal(0, 100);
  #   beta_1 ~ normal(0, 100);
  #   gamma[2:6] ~ normal(rho * gamma[1:(T - 1)], sigma);
  #   n_deaths ~ binomial(n_births, p);
  # }")
  # stan_data <- region_data_trans %>% 
  #   subset( select=-c(region,years,log_mbias) ) %>%
  #   filter(!is.na(n_deaths) & !is.na(u5mr)) %>%
  #   as.list()
  # stan_data$N = length(stan_data$id_reg)
  # fit <- stan(
  #   model_code = stan_code,
  #   data = stan_data
  # )
  # print(fit)
  
  
  
  # # Old models for log(MBIAS)
  # 
  # # Model #1
  # model_1 <- inla(
  #   log_mbias ~ 1,
  #   family = "gaussian",
  #   data = region_data
  # )
  # summary(model_1)
  # 
  # # Model #2
  # model_2 <- inla(
  #   log_mbias ~ time,
  #   family = "gaussian",
  #   data = region_data
  # )
  # summary(model_2)
  # 
  # # Model #3
  # model_3 <- inla(
  #   log_mbias ~ f(time, model="ar1"),
  #   family = "gaussian",
  #   data = region_data
  # )
  # summary(model_3)
  # 
  # # Model #4
  # model_4 <- inla(
  #   log_mbias ~ f(id_reg,model="iid"),
  #   family = "gaussian",
  #   data = region_data
  # )
  # summary(model_4)
  # 
  # # Model #5
  # model_5 <- inla(
  #   log_mbias ~ f(id_reg,model="bym2",graph=mat),
  #   family = "gaussian",
  #   data = region_data
  # )
  # summary(model_5)
  # 
  # # Model #6
  # model_6 <- inla(
  #   log_mbias ~ f(time, model="ar1") +
  #     f(id_reg,model="bym2",graph=mat),
  #   family = "gaussian",
  #   data = region_data
  # )
  # summary(model_6)
  # 
  # 
  # 
  # # INLA testing
  # 
  # # Model #1
  # df_1 <- data.frame(
  #   "geo_id" = c(1,2,3,4,5,6,7,8,9,10),
  #   "n_births" = c(100,100,100,100,100,100,100,100,100,100),
  #   "n_deaths" = c(5,5,5,6,4,5,7,5,3,5)
  # )
  # model_1 <- inla(
  #   n_deaths ~ 1,
  #   family = "binomial",
  #   Ntrials = df_1$n_birth,
  #   data = df_1,
  #   control.family = list(link="logit"),
  #   control.predictor = list(compute=TRUE)
  # )
  # summary_1 <- summary(model_1)
  # 
  # # Actual mean
  # sum(df_1$n_deaths)/sum(df_1$n_births)
  # 
  # # Predicted mean
  # expit(summary_1$fixed[1,"mean"])
  # 
  # # Model #2
  # df_2 <- data.frame(
  #   "geo_id" = c(1,2,3,4,5,1,2,3,4,5),
  #   "n_births" = c(100,100,100,100,100,100,100,100,100,100),
  #   "n_deaths" = c(5,5,5,6,4,15,16,15,14,15),
  #   "source" = c(0,0,0,0,0,1,1,1,1,1) # 0=vr, 1=survey
  # )
  # model_2 <- inla(
  #   n_deaths ~ 1 + source,
  #   family = "binomial",
  #   Ntrials = df_2$n_birth,
  #   data = df_2,
  #   control.family = list(link="logit"),
  #   control.predictor = list(compute=TRUE)
  # )
  # summary_2 <- summary(model_2)
  # 
  # # Actual mean, by source (0,1)
  # sum(df_2[1:5,]$n_deaths)/sum(df_2[1:5,]$n_births)
  # sum(df_2[6:10,]$n_deaths)/sum(df_2[6:10,]$n_births)
  # 
  # # Predicted mean, by covariate group (0,1)
  # beta_intercept <- summary_2$fixed[1,"mean"]
  # beta_covariate <- summary_2$fixed[2,"mean"]
  # expit(beta_intercept)
  # expit(beta_intercept + beta_covariate)
  # 
  # 
  # 
  # # Plot Ecuador mortality rates over time by region with CIs
  # regions <- unique(u5m_regional$region)
  # u5m_regional$years <- factor(u5m_regional$years, levels = years_3)
  # u5m_regional$region <- factor(u5m_regional$region, levels = regions)
  # ggplot(u5m_regional, aes(x = years, y = u5m, ymin = lower,
  #                          ymax = upper, color = region)) +
  #   geom_point() +
  #   geom_errorbar(width = 0.2) +
  #   facet_wrap(~region, ncol = 6) +
  #   theme(legend.position="none")
  # 
  # 
  # 
  # # Adapted from Wakefield Quito workshop code
  # # National-level estimates
  # births_one_year <- getBirths(
  #   data = fbh,
  #   surveyyear = 2012,
  #   variables = c("idhog", "province", "V005", "V021", "V024", "V025"),
  #   strata = "V024",
  #   dob = "B3",
  #   alive = "alive",
  #   age = "B7",
  #   date.interview = "V008",
  #   year.cut = seq(1992, 2012, by = 1)
  # )
  # years_1 <- levels(births_one_year$time)
  # u5m_country <- getDirect(
  #   births = births_one_year,
  #   years = years_1,
  #   regionVar = "province",
  #   timeVar = "time",
  #   clusterVar = "~V021+idhog",
  #   ageVar = "age",
  #   weightsVar = "V005",
  #   geo.recode = NULL,
  #   national.only = TRUE
  # )
  # print(u5m_country$u5m)
  # print(u5m_country$lower)
  # print(u5m_country$upper)
  # 
  # 
  # 
  # # Function: Update data frame with survey data
  # update_df_SVY = function(df_VR, df_survey, type) {
  # 
  #   for (j in 1:nrow(df_survey)) {
  # 
  #     if (type=="national") {
  # 
  #       # Calculate next table index
  #       i = nrow(df_VR) + 1
  # 
  #       # Parse year
  #       year = as.numeric(substr(df_survey[j,"years"],1,2))+1900
  #       year = ifelse(year<1950,year+100,year)
  # 
  #       # Add national row
  #       df_VR[i,"region"] = "National"
  #       df_VR[i,"source"] = "Survey"
  #       df_VR[i,"year"] = year
  #       df_VR[i,"u5mr"] = df_survey[j,"u5m"]
  # 
  #     }
  # 
  #     if (type == "regional") {
  # 
  #       # Calculate next table index
  #       i = nrow(df_VR) + 1
  # 
  #       # Add regional rows
  #       region = levels(df_survey[j,"region"])[as.numeric(df_survey[j,"region"])]
  #       df_VR[i,"region"] = region
  #       df_VR[i,"source"] = "Survey"
  #       df_VR[i,"year"] = df_survey[j,"years"]
  #       df_VR[i,"u5mr"] = df_survey[j,"u5m"]
  # 
  #     }
  # 
  #   }
  # 
  #   # Return updated table
  #   return (df_VR)
  # 
  # }
  # 
  # 
  # 
  # # Add survey data to df_VR
  # df_VR = update_df_SVY(df_VR, u5m_country, "national")
  # u5m_regional %<>% filter(region != "All")
  # df_VR = update_df_SVY(df_VR, u5m_regional, "regional")
  # 
  # 
  # 
  # # Generate "groups" (note: different than dplyr "groups" above)
  # df_VR %<>% mutate(
  #   group = paste(region,source, sep=": ")
  # )
  # 
  # 
  # 
  # # Generate plots
  # df_n = df_VR %>% filter(region == "National")
  # ggplot(data=df_n, aes(x=year, y=u5mr, group=group, color=group)) +
  #   geom_line() + geom_point()
  # df_r = df_VR %>% filter(
  #   region != "National" & year %in% c("02-06","07-11")
  # )
  # ggplot(data=df_r, aes(
  #   x = year,
  #   y = u5mr,
  #   group = group,
  #   color = substr(group, nchar(group)-1, nchar(group))
  # )) +
  #   geom_line() + geom_point() +
  #   facet_wrap(~ region, ncol=2)
  
}
