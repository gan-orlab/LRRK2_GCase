# In this project we used various lineal regression models to test for association between _LRRK2_ variants and GCase activity in the Columbia and PPMI cohorts.
1. In the Columbia cohort we adjusted for age, sex, PD status, _GBA_ status and ethnicity
2. In the PPMI cohort we adjusted for age, sex, PD status, _GBA_ status, ethnicity and white blood cells count
3. For both cohorts we repeated analysis for	 1)PD cases+ controls 2)PD cases 3) Controls 4) Excluding carriers of LRRK2 p.G2019S 5) Excluding carriers of GBA variants and LRRK2 p.G2019S variant 6) Analyses stratifying the cohorts by sex

## Linear regression models used to test for association between _LRRK2_ variants and GCase activity in the Columbia

    final <- merge(covar, gcase, by.x = "FID2", by.y = "FID2" )
    final2 <- merge(final, var, by.x = "FID2", by.y = "FID2" )
    assoc <- function(x,n,covar_assoc) {
      x <- as.numeric(x)
      if(!(any(x %in% c(1,2)))){return(NULL)}
      fit <- lm(covar_assoc$Gcase ~ x + covar_assoc$Sex + covar_assoc$Age.y + covar_assoc$Status.y + covar_assoc$GBA_coding + covar_assoc$Ethn)
      output <- cbind(coef(summary(fit))[c(1,2),],confint.default(fit)[c(1,2),])
      no_na_covar <- covar_assoc[!is.na(covar_assoc[,n]),]
      output <- cbind(output,mean(no_na_covar[no_na_covar[,n] > 0,]$Gcase, na.rm = T))
      colnames(output)[7] <- "Gcase"
      row.names(output)[2] <- names(covar_assoc)[n]
      output <- cbind(output,length(which(no_na_covar[,n] > 0)))
      colnames(output)[8] <- "N of carriers"
      output <- cbind(output,sd(no_na_covar[no_na_covar[,n] > 0,]$Gcase, na.rm = T))
      colnames(output)[9] <- "Gcase_SD"
      return(output[2,,drop = FALSE])
    }
    
    final_output <- Reduce(rbind,lapply(9:ncol(final2), function(x){
     assoc(final2[,x],x,final2)
    }))
    final_output <- cbind(final_output,(final_output[,"Gcase"]-mean(final_output[,"Gcase"]))/sd(final_output[,"Gcase"]))
    colnames(final_output)[10] <- "Gcase_Zscore"
    write.csv(final_output,file="PD_controls.csv",quote = FALSE)
    
    #selecting only PD cases
    case_covar <- final2[final2$Status.y == 2,]
    
    assoc <- function(x,n,covar_assoc) {
      x <- as.numeric(x)
      if(!(any(x %in% c(1,2)))){return(NULL)}
      fit <- lm(covar_assoc$Gcase ~ x + covar_assoc$Sex + covar_assoc$Age.y + covar_assoc$Status.y + covar_assoc$GBA_coding + covar_assoc$Ethn)
      output <- cbind(coef(summary(fit))[c(1,2),],confint.default(fit)[c(1,2),])
      no_na_covar <- covar_assoc[!is.na(covar_assoc[,n]),]
      output <- cbind(output,mean(no_na_covar[no_na_covar[,n] > 0,]$Gcase, na.rm = T))
      colnames(output)[7] <- "Gcase"
      row.names(output)[2] <- names(covar_assoc)[n]
      output <- cbind(output,length(which(no_na_covar[,n] > 0)))
      colnames(output)[8] <- "N of carriers"
      output <- cbind(output,sd(no_na_covar[no_na_covar[,n] > 0,]$Gcase, na.rm = T))
      colnames(output)[9] <- "Gcase_SD"
      return(output[2,,drop = FALSE])
    }
    
    final_output <- Reduce(rbind,lapply(9:ncol(case_covar), function(x){
     assoc(case_covar[,x],x,case_covar)
    }))
    
    final_output <- cbind(final_output,(final_output[,"Gcase"]-mean(final_output[,"Gcase"]))/sd(final_output[,"Gcase"]))
    colnames(final_output)[10] <- "Gcase_Zscore"
    write.csv(final_output,file="PD.csv",quote = FALSE)
    
    
    #selecting only controls
    control_covar <- final2[final2$Status.y == 1,]
    assoc <- function(x,n,covar_assoc) {
      x <- as.numeric(x)
      if(!(any(x %in% c(1,2)))){return(NULL)}
      fit <- lm(covar_assoc$Gcase ~ x + covar_assoc$Sex + covar_assoc$Age.y + covar_assoc$Status.y + covar_assoc$GBA_coding + covar_assoc$Ethn)
      output <- cbind(coef(summary(fit))[c(1,2),],confint.default(fit)[c(1,2),])
      no_na_covar <- covar_assoc[!is.na(covar_assoc[,n]),]
      output <- cbind(output,mean(no_na_covar[no_na_covar[,n] > 0,]$Gcase, na.rm = T))
      colnames(output)[7] <- "Gcase"
      row.names(output)[2] <- names(covar_assoc)[n]
      output <- cbind(output,length(which(no_na_covar[,n] > 0)))
      colnames(output)[8] <- "N of carriers"
      output <- cbind(output,sd(no_na_covar[no_na_covar[,n] > 0,]$Gcase, na.rm = T))
      colnames(output)[9] <- "Gcase_SD"
      return(output[2,,drop = FALSE])
    }
    
    final_output <- Reduce(rbind,lapply(9:ncol(control_covar), function(x){
     assoc(control_covar[,x],x,control_covar)
    }))
    
    final_output <- cbind(final_output,(final_output[,"Gcase"]-mean(final_output[,"Gcase"]))/sd(final_output[,"Gcase"]))
    colnames(final_output)[10] <- "Gcase_Zscore"
    write.csv(final_output,file="controls.csv",quote = FALSE)
    
    #Excluding carriers of LRRK2 
    LRRK2_status <- final2[final2$G2019S == 0,]
    assoc <- function(x,n,covar_assoc) {
      x <- as.numeric(x)
      if(!(any(x %in% c(1,2)))){return(NULL)}
      fit <- lm(covar_assoc$Gcase ~ x + covar_assoc$Sex + covar_assoc$Age.y + covar_assoc$Status.y + covar_assoc$GBA_coding + covar_assoc$Ethn)
      output <- cbind(coef(summary(fit))[c(1,2),],confint.default(fit)[c(1,2),])
      no_na_covar <- covar_assoc[!is.na(covar_assoc[,n]),]
      output <- cbind(output,mean(no_na_covar[no_na_covar[,n] > 0,]$Gcase, na.rm = T))
      colnames(output)[7] <- "Gcase"
      row.names(output)[2] <- names(covar_assoc)[n]
      output <- cbind(output,length(which(no_na_covar[,n] > 0)))
      colnames(output)[8] <- "N of carriers"
      output <- cbind(output,sd(no_na_covar[no_na_covar[,n] > 0,]$Gcase, na.rm = T))
      colnames(output)[9] <- "Gcase_SD"
      return(output[2,,drop = FALSE])
    }
    
    final_output <- Reduce(rbind,lapply(9:ncol(LRRK2_status), function(x){
     assoc(LRRK2_status[,x],x,LRRK2_status)
    }))
    
    final_output <- cbind(final_output,(final_output[,"Gcase"]-mean(final_output[,"Gcase"]))/sd(final_output[,"Gcase"]))
    colnames(final_output)[10] <- "Gcase_Zscore"
    write.csv(final_output,file="excludingLRRK2G2019Scarriers.csv",quote = FALSE)
    
    
    #Excluding carriers of LRRK2 and GBA
    GBA_LRRK2_status <- final2[final2$GBA_coding == 0 & final2$G2019S == 0,]
    
    assoc <- function(x,n,covar_assoc) {
      x <- as.numeric(x)
      if(!(any(x %in% c(1,2)))){return(NULL)}
      fit <- lm(covar_assoc$Gcase ~ x + covar_assoc$Sex + covar_assoc$Age.y + covar_assoc$Status.y + covar_assoc$GBA_coding + covar_assoc$Ethn)
      output <- cbind(coef(summary(fit))[c(1,2),],confint.default(fit)[c(1,2),])
      no_na_covar <- covar_assoc[!is.na(covar_assoc[,n]),]
      output <- cbind(output,mean(no_na_covar[no_na_covar[,n] > 0,]$Gcase, na.rm = T))
      colnames(output)[7] <- "Gcase"
      row.names(output)[2] <- names(covar_assoc)[n]
      output <- cbind(output,length(which(no_na_covar[,n] > 0)))
      colnames(output)[8] <- "N of carriers"
      output <- cbind(output,sd(no_na_covar[no_na_covar[,n] > 0,]$Gcase, na.rm = T))
      colnames(output)[9] <- "Gcase_SD"
      return(output[2,,drop = FALSE])
    }
    
    final_output <- Reduce(rbind,lapply(9:ncol(GBA_LRRK2_status), function(x){
     assoc(GBA_LRRK2_status[,x],x,GBA_LRRK2_status)
    }))
    
    final_output <- cbind(final_output,(final_output[,"Gcase"]-mean(final_output[,"Gcase"]))/sd(final_output[,"Gcase"]))
    colnames(final_output)[10] <- "Gcase_Zscore"
    write.csv(final_output,file="excludingLRRK2andGBA.csv",quote = FALSE)
    
    
    #Selecting only males
    sex2 <- final2[final2$Sex == 2,]
    assoc <- function(x,n,covar_assoc) {
      x <- as.numeric(x)
      if(!(any(x %in% c(1,2)))){return(NULL)}
      fit <- lm(covar_assoc$Gcase ~ x + covar_assoc$Sex + covar_assoc$Age.y + covar_assoc$Status.y + covar_assoc$GBA_coding + covar_assoc$Ethn)
      output <- cbind(coef(summary(fit))[c(1,2),],confint.default(fit)[c(1,2),])
      no_na_covar <- covar_assoc[!is.na(covar_assoc[,n]),]
      output <- cbind(output,mean(no_na_covar[no_na_covar[,n] > 0,]$Gcase, na.rm = T))
      colnames(output)[7] <- "Gcase"
      row.names(output)[2] <- names(covar_assoc)[n]
      output <- cbind(output,length(which(no_na_covar[,n] > 0)))
      colnames(output)[8] <- "N of carriers"
      output <- cbind(output,sd(no_na_covar[no_na_covar[,n] > 0,]$Gcase, na.rm = T))
      colnames(output)[9] <- "Gcase_SD"
      return(output[2,,drop = FALSE])
    }
    
    final_output <- Reduce(rbind,lapply(9:ncol(sex2), function(x){
     assoc(sex2[,x],x,sex2)
    }))
    
    final_output <- cbind(final_output,(final_output[,"Gcase"]-mean(final_output[,"Gcase"]))/sd(final_output[,"Gcase"]))
    colnames(final_output)[10] <- "Gcase_Zscore"
    write.csv(final_output,file="sex2.csv",quote = FALSE)
    
    
    #Selecting only females
    sex1 <- final2[final2$Sex == 1,]
    assoc <- function(x,n,covar_assoc) {
      x <- as.numeric(x)
      if(!(any(x %in% c(1,2)))){return(NULL)}
      fit <- lm(covar_assoc$Gcase ~ x + covar_assoc$Sex + covar_assoc$Age.y + covar_assoc$Status.y + covar_assoc$GBA_coding + covar_assoc$Ethn)
      output <- cbind(coef(summary(fit))[c(1,2),],confint.default(fit)[c(1,2),])
      no_na_covar <- covar_assoc[!is.na(covar_assoc[,n]),]
      output <- cbind(output,mean(no_na_covar[no_na_covar[,n] > 0,]$Gcase, na.rm = T))
      colnames(output)[7] <- "Gcase"
      row.names(output)[2] <- names(covar_assoc)[n]
      output <- cbind(output,length(which(no_na_covar[,n] > 0)))
      colnames(output)[8] <- "N of carriers"
      output <- cbind(output,sd(no_na_covar[no_na_covar[,n] > 0,]$Gcase, na.rm = T))
      colnames(output)[9] <- "Gcase_SD"
      return(output[2,,drop = FALSE])
    }
    
    final_output <- Reduce(rbind,lapply(9:ncol(sex1), function(x){
     assoc(sex1[,x],x,sex1)
    }))
    
    final_output <- cbind(final_output,(final_output[,"Gcase"]-mean(final_output[,"Gcase"]))/sd(final_output[,"Gcase"]))
    colnames(final_output)[10] <- "Gcase_Zscore"
    write.csv(final_output,file="sex1.csv",quote = FALSE)

## Linear regression models used to test for association between _LRRK2_ variants and GCase activity in the PPMI cohorts

    covar <- read.csv("PPMI_covar_final.csv", header = T, sep=",")
    assoc <- function(x,n,covar) {
          x <- as.numeric(x)
          #  x[x==2] <- 1
          if(!(any(x %in% c(1,2)))){return(NULL)}
          fit <- lm(covar$GCase_average ~ x + covar$Sex + covar$Age + covar$Status + covar$GBA_coding+covar$Ethn+covar$WBC_mean)
          output <- cbind(coef(summary(fit))[c(1,2),],confint.default(fit)[c(1,2),])
          no_na_covar <- covar[!is.na(covar[,n]),]
          output <- cbind(output,mean(no_na_covar[no_na_covar[,n] > 0,]$GCase_average, na.rm = T))
          colnames(output)[7] <- "GCase_average"
          row.names(output)[2] <- names(covar)[n]
          output <- cbind(output,length(which(no_na_covar[,n] > 0)))
          colnames(output)[8] <- "N of carriers"
          output <- cbind(output,sd(no_na_covar[no_na_covar[,n] > 0,]$GCase_average, na.rm = T))
          colnames(output)[9] <- "Gcase_SD"
          #row.names(output) <- names(var)[n]
          return(output[2,,drop = FALSE])
    }
    
    final_output <- Reduce(rbind,lapply(10:ncol(covar ), function(x){
          assoc(covar [,x],x,covar )
    }))
    final_output <- cbind(final_output,(final_output[,"GCase_average"]-mean(final_output[,"GCase_average"]))/sd(final_output[,"GCase_average"]))
    colnames(final_output)[9] <- "Gcase_Zscore"
    write.csv(final_output,file="PD+controls.csv",quote = FALSE)
    
    #Selecting only PD cases
    case_covar <- covar[covar$Pheno == 2,]
    assoc <- function(x,n,covar) {
          x <- as.numeric(x)
          #  x[x==2] <- 1
          if(!(any(x %in% c(1,2)))){return(NULL)}
          fit <- lm(covar$GCase_average ~ x + covar$Sex + covar$Age + covar$Status + covar$GBA_coding+covar$Ethn+covar$WBC_mean)
          output <- cbind(coef(summary(fit))[c(1,2),],confint.default(fit)[c(1,2),])
          no_na_covar <- covar[!is.na(covar[,n]),]
          output <- cbind(output,mean(no_na_covar[no_na_covar[,n] > 0,]$GCase_average, na.rm = T))
          colnames(output)[7] <- "GCase_average"
          row.names(output)[2] <- names(covar)[n]
          output <- cbind(output,length(which(no_na_covar[,n] > 0)))
          colnames(output)[8] <- "N of carriers"
          output <- cbind(output,sd(no_na_covar[no_na_covar[,n] > 0,]$GCase_average, na.rm = T))
          colnames(output)[9] <- "Gcase_SD"
          #row.names(output) <- names(var)[n]
          return(output[2,,drop = FALSE])
    }
    
    final_output <- Reduce(rbind,lapply(10:ncol(case_covar), function(x){
          assoc(case_covar[,x],x,case_covar)
    }))
    
    final_output <- cbind(final_output,(final_output[,"GCase_average"]-mean(final_output[,"GCase_average"]))/sd(final_output[,"GCase_average"]))
    colnames(final_output)[9] <- "Gcase_Zscore"
    write.csv(final_output,file="PD.csv",quote = FALSE)
    
    #Selecting only controls
    control_covar <- covar[covar$Pheno == 1,]
    
    assoc <- function(x,n,covar) {
          x <- as.numeric(x)
          #  x[x==2] <- 1
          if(!(any(x %in% c(1,2)))){return(NULL)}
          fit <- lm(covar$GCase_average ~ x + covar$Sex + covar$Age + covar$Status + covar$GBA_coding+covar$Ethn+covar$WBC_mean)
          output <- cbind(coef(summary(fit))[c(1,2),],confint.default(fit)[c(1,2),])
          no_na_covar <- covar[!is.na(covar[,n]),]
          output <- cbind(output,mean(no_na_covar[no_na_covar[,n] > 0,]$GCase_average, na.rm = T))
          colnames(output)[7] <- "GCase_average"
          row.names(output)[2] <- names(covar)[n]
          output <- cbind(output,length(which(no_na_covar[,n] > 0)))
          colnames(output)[8] <- "N of carriers"
          output <- cbind(output,sd(no_na_covar[no_na_covar[,n] > 0,]$GCase_average, na.rm = T))
          colnames(output)[9] <- "Gcase_SD"
          #row.names(output) <- names(var)[n]
          return(output[2,,drop = FALSE])
    }
    
    final_output <- Reduce(rbind,lapply(10:ncol(covar ), function(x){
          assoc(covar [,x],x,covar )
    }))
    final_output <- cbind(final_output,(final_output[,"GCase_average"]-mean(final_output[,"GCase_average"]))/sd(final_output[,"GCase_average"]))
    colnames(final_output)[11] <- "Gcase_Zscore"
    write.csv(final_output,file="controls.csv",quote = FALSE)
    
    #excluding carriers of LRRK2G2019S
    LRRK2_status <- covar[covar$G2019S == 0,]
    
    assoc <- function(x,n,covar) {
          x <- as.numeric(x)
          #  x[x==2] <- 1
          if(!(any(x %in% c(1,2)))){return(NULL)}
          fit <- lm(covar$GCase_average ~ x + covar$Sex + covar$Age + covar$Status + covar$GBA_coding+covar$Ethn+covar$WBC_mean)
          output <- cbind(coef(summary(fit))[c(1,2),],confint.default(fit)[c(1,2),])
          no_na_covar <- covar[!is.na(covar[,n]),]
          output <- cbind(output,mean(no_na_covar[no_na_covar[,n] > 0,]$GCase_average, na.rm = T))
          colnames(output)[7] <- "GCase_average"
          row.names(output)[2] <- names(covar)[n]
          output <- cbind(output,length(which(no_na_covar[,n] > 0)))
          colnames(output)[8] <- "N of carriers"
          output <- cbind(output,sd(no_na_covar[no_na_covar[,n] > 0,]$GCase_average, na.rm = T))
          colnames(output)[9] <- "Gcase_SD"
          #row.names(output) <- names(var)[n]
          return(output[2,,drop = FALSE])
    }
    
    final_output <- Reduce(rbind,lapply(10:ncol(covar ), function(x){
          assoc(covar [,x],x,covar )
    }))
    final_output <- cbind(final_output,(final_output[,"GCase_average"]-mean(final_output[,"GCase_average"]))/sd(final_output[,"GCase_average"]))
    colnames(final_output)[11] <- "Gcase_Zscore"
    write.csv(final_output,file="PD+controls_excludingLRRK2.csv",quote = FALSE)
    #excluding carriers of GBA and LRRK2 p.G2019S
    GBA_LRRK2_status <- covar[covar$GBA_coding == 0 & covar$G2019S == 0,]
    assoc <- function(x,n,covar) {
          x <- as.numeric(x)
          #  x[x==2] <- 1
          if(!(any(x %in% c(1,2)))){return(NULL)}
          fit <- lm(covar$GCase_average ~ x + covar$Sex + covar$Age + covar$Status + covar$GBA_coding+covar$Ethn+covar$WBC_mean)
          output <- cbind(coef(summary(fit))[c(1,2),],confint.default(fit)[c(1,2),])
          no_na_covar <- covar[!is.na(covar[,n]),]
          output <- cbind(output,mean(no_na_covar[no_na_covar[,n] > 0,]$GCase_average, na.rm = T))
          colnames(output)[7] <- "GCase_average"
          row.names(output)[2] <- names(covar)[n]
          output <- cbind(output,length(which(no_na_covar[,n] > 0)))
          colnames(output)[8] <- "N of carriers"
          output <- cbind(output,sd(no_na_covar[no_na_covar[,n] > 0,]$GCase_average, na.rm = T))
          colnames(output)[9] <- "Gcase_SD"
          #row.names(output) <- names(var)[n]
          return(output[2,,drop = FALSE])
    }
    
    final_output <- Reduce(rbind,lapply(10:ncol(covar ), function(x){
          assoc(covar [,x],x,covar )
    }))
    final_output <- cbind(final_output,(final_output[,"GCase_average"]-mean(final_output[,"GCase_average"]))/sd(final_output[,"GCase_average"]))
    colnames(final_output)[11] <- "Gcase_Zscore"
    write.csv(final_output,file="PD+controls_excludingLRRK2andGBA.csv",quote = FALSE)
    
    #select only males
    sex2 <- final2[final2$Sex == 2,]
    assoc <- function(x,n,covar) {
          x <- as.numeric(x)
          #  x[x==2] <- 1
          if(!(any(x %in% c(1,2)))){return(NULL)}
          fit <- lm(covar$GCase_average ~ x + covar$Sex + covar$Age + covar$Status + covar$GBA_coding+covar$Ethn+covar$WBC_mean)
          output <- cbind(coef(summary(fit))[c(1,2),],confint.default(fit)[c(1,2),])
          no_na_covar <- covar[!is.na(covar[,n]),]
          output <- cbind(output,mean(no_na_covar[no_na_covar[,n] > 0,]$GCase_average, na.rm = T))
          colnames(output)[7] <- "GCase_average"
          row.names(output)[2] <- names(covar)[n]
          output <- cbind(output,length(which(no_na_covar[,n] > 0)))
          colnames(output)[8] <- "N of carriers"
          output <- cbind(output,sd(no_na_covar[no_na_covar[,n] > 0,]$GCase_average, na.rm = T))
          colnames(output)[9] <- "Gcase_SD"
          #row.names(output) <- names(var)[n]
          return(output[2,,drop = FALSE])
    }
    
    final_output <- Reduce(rbind,lapply(10:ncol(covar ), function(x){
          assoc(covar [,x],x,covar )
    }))
    final_output <- cbind(final_output,(final_output[,"GCase_average"]-mean(final_output[,"GCase_average"]))/sd(final_output[,"GCase_average"]))
    colnames(final_output)[11] <- "Gcase_Zscore"
    write.csv(final_output,file="only_man.csv",quote = FALSE)
    
    #select females only
    sex1 <- final2[final2$Sex == 1,]
    
    assoc <- function(x,n,covar) {
          x <- as.numeric(x)
          #  x[x==2] <- 1
          if(!(any(x %in% c(1,2)))){return(NULL)}
          fit <- lm(covar$GCase_average ~ x + covar$Sex + covar$Age + covar$Status + covar$GBA_coding+covar$Ethn+covar$WBC_mean)
          output <- cbind(coef(summary(fit))[c(1,2),],confint.default(fit)[c(1,2),])
          no_na_covar <- covar[!is.na(covar[,n]),]
          output <- cbind(output,mean(no_na_covar[no_na_covar[,n] > 0,]$GCase_average, na.rm = T))
          colnames(output)[7] <- "GCase_average"
          row.names(output)[2] <- names(covar)[n]
          output <- cbind(output,length(which(no_na_covar[,n] > 0)))
          colnames(output)[8] <- "N of carriers"
          output <- cbind(output,sd(no_na_covar[no_na_covar[,n] > 0,]$GCase_average, na.rm = T))
          colnames(output)[9] <- "Gcase_SD"
          #row.names(output) <- names(var)[n]
          return(output[2,,drop = FALSE])
    }
    
    final_output <- Reduce(rbind,lapply(10:ncol(covar ), function(x){
          assoc(covar [,x],x,covar )
    }))
    final_output <- cbind(final_output,(final_output[,"GCase_average"]-mean(final_output[,"GCase_average"]))/sd(final_output[,"GCase_average"]))
    colnames(final_output)[11] <- "Gcase_Zscore"
    write.csv(final_output,file="onlyfemale.csv",quote = FALSE)


