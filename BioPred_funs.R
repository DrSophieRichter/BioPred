#BioPred functions ##############

###########################################
# logit helper function
###########################################

logit <- LaplacesDemon::logit

###########################################
#rendering functions for table1
###########################################

#set up display of continuous variables
my.render.cont <- function(x) 
{
  #display 2 sig figs:
  with(stats.apply.rounding(stats.default(x), digits=2), 
       #in first column display variable name and "Median (IQR)":
       c("", "Median (Q1-Q3)" = 
           #in second column calculate and display:
           sprintf("%s (%s %s %s)", 
                   MEDIAN,  
                   round(quantile(x, 0.25, na.rm = T),1) , 
                   "-", 
                   round(quantile(x, 0.75, na.rm = T),1))))   #calculations requested
}


#set up display of categorical variables
my.render.cat <- function(x) 
{
  #in first column display variable name:
  c("", 
    #in second column calculate and display:
    sapply(stats.default(x), function(y) with(y, sprintf("%d (%0.0f %%)", 
                                                         FREQ, 
                                                         round(PCT, 0)))))
}


#########################################################
#
# Adjust biomarker concentrations for time of sampling
#
########################################################

adjust_bio <- function(imp, imp_list, mybio){
  
  #Fit model to multiply imputed data
  #that predicts biomarkers had they all been sampled at the same time
  fit <- 
    imp %>%
    with(lm(get(mybio) ~ 
              poly(Time_to_bio,3) +
              Age +
              Sex +
              Latest_GCS +
              Motor_score +
              Reactive_pupils +
              Marshall_score +
              CT_SAH +
              CT_EDH +
              CT_TAI +
              CT_Cistern +
              CT_Mids +
              CT_Haematoma +
              MajorEI +
              Hypoxia +
              Hypotension +
              Admission_status +
              Alcohol_intoxication))
  pooled_fit <- pool(fit)
  
  #put pooled coefficients into an lm object to make predictions
  pooled_lm <- fit$analyses[[1]]
  # Replace the fitted coefficients with the pooled
  #   estimates (need to check they are replaced in
  #   the correct order)
  pooled_lm$coefficients <- summary(pooled_fit)$estimate
  
  
  # WITHIN EACH IMPUTED DATASET
  # Predict the difference between 2h and actual prediction time
  # then add this delta (which may be negative for late samples)
  # to the observed concentration to obtain the time adjusted
  # biomarker result
  for (i in 1:length(imp_list)){
    
    mydata_orig <- imp_list[[i]]
    mydata_2h   <- imp_list[[i]]
    mydata_2h$Time_to_bio <- 2
    
    bio_pred    <- predict(pooled_lm, newdata = mydata_orig)
    bio_2h      <- predict(pooled_lm, newdata = mydata_2h)
    bio_delta   <- bio_2h - bio_pred
    newname     <- gsub("_unadj", "", mybio)
    imp_list[[i]][,newname] <- imp_list[[i]][, paste0(mybio)] + bio_delta
    
  }
  
  return(imp_list)
}



##########################################################
#
# Function to select the optimal panel of biomarkers
# to predict chosen outcome
# considering multiply imputed data
#
##########################################################

select_panel <- function(imputed_list = imp_list, imputations = m){
  
  myoutcome <- "Outcome_5plus"
  
  #initialize vector to store lambdas from m imputed datasets
  my_lambdas <- vector()
  
  #variables to consider for selection
  cols = c('GFAP', 'NFL', 'NSE', 'S100B', "Tau", "UCH.L1")
  
  #1.	Find one lambda per imputed dataset (using 10 x 10fold cross-validation). 
  #This will yield 10 lambdas
  for (imp_num in 1:imputations) {
    #select one imputed dataset
    myimp        <- 
      imputed_list[[imp_num]] %>% 
      select(all_of(myoutcome), 
             GFAP:UCH.L1)
    
    #scale and center data
    pre_proc_val <- preProcess(myimp[,cols], method = c("center", "scale"))
    myimp[,cols] <- predict(pre_proc_val, myimp[,cols])
    
    #make model matrix
    x <-  as.matrix(myimp[, -1])
    y <-  myimp[,1]
    
    #determine lambda using 10x repeated 10fold cross-validation
    #note glmnet standardizes variables by default
    MSEs <- NULL
    for (i in 1:10){
      cv <- cv.glmnet(x = x, y = y, alpha=1, nfolds=10)  
      MSEs <- cbind(MSEs, cv$cvm)
    }
    rownames(MSEs) <- cv$lambda
    lambda_min <- as.numeric(names(which.min(rowMeans(MSEs))))
    my_lambdas[[imp_num]] <- lambda_min
  }
  
  
  #2.	I take the mean lambda as my optimal_lambda
  optimal_lambda <- mean(my_lambdas) 
  
  #3.	I stack all 10 imputed datasets
  mystack        <- 
    bind_rows(imputed_list) %>% 
    select(all_of(myoutcome), 
           GFAP:UCH.L1)
  
  #scale and center stacked data
  pre_proc_val <- preProcess(mystack[,cols], method = c("center", "scale"))
  mystack[,cols] <- predict(pre_proc_val, mystack[,cols])
  
  #make model matrix
  x <-  as.matrix(mystack[, -1])
  y <-  mystack[,1]
  
  #4.	I fit one lasso model to the stacked data, using 
  #    a.	Weights = 0.1
  #    b.	Lambda = optimal_lambda
  myweights <- rep((1/imputations), nrow(mystack))
  lasso_model <- glmnet(x, y, 
                        alpha = 1, #for lasso
                        lambda = optimal_lambda, #mean of 10 imputed datasets
                        weights = myweights #weighted to not inflate sample size
  )
  #Extract coefficients
  mycoefs <- 
    coef(lasso_model) %>%
    as.matrix() %>%
    as.data.frame()
  
  colnames(mycoefs) <- myoutcome
  
  return(mycoefs)
  
}

##########################################################
#
# Format lasso results table
# and save to word
#
##########################################################

format_and_save_panel_table <- function(temp){
  
  colnames(temp) <- c("all","m12", "m3456")
  temp <- 
    temp %>% 
    mutate_if(is.numeric, round, digits = 2)
  
  mymax <- max(temp[-1,])
  myhalf <<- mymax/2
  
  
  temp$Variable <- rownames(temp)
  temp <- temp %>% select(Variable, everything())
  temp$Variable <- gsub(".", "-", temp$Variable, fixed = TRUE)
  
  colourer <- scales::col_numeric(
    palette = "BrBG",
    domain = c(-mymax, mymax))
  
  ft_panel <-
    temp %>%
    regulartable() %>% 
    set_header_labels(all = "Overall",
                      m12 = "Marshall <3", 
                      m3456 = "Marshall >=3") %>%
    bg(
      bg = colourer,
      i = c(rownames(temp)[-1]),
      j = c(colnames(temp)[-1]),
      part = "body") %>%
    color(i = ~ all >= myhalf, j = "all", color = "white") %>%
    color(i = ~ m12 >= myhalf, j = "m12", color = "white") %>%
    color(i = ~ m3456 >= myhalf, j = "m3456", color = "white") %>%
    color(i = 1, color = "black") %>%
    vline(j = 1, border = fp_border_default()) %>%
    vline_right(part = "all") %>%
    vline_left(part = "all") %>%
    fix_border_issues(part = "all") 
  
  save_as_docx(ft_panel, path = "Figures/Panel_selection.docx")
  
}

###########################################################
#
#Pool coefficients across imputed datasets
#
###########################################################

pool_my_coefs <- function(estimate_list, my_imp_list) {
  mycols <- colnames(estimate_list[[1]])
  myrows <- rownames(as.data.frame(estimate_list[[1]]))
  arr_estimate <- array(unlist(estimate_list) , 
                        c(length(myrows),
                          length(mycols),
                          length(my_imp_list)) )
  pooled_estimate <- apply(arr_estimate, 1:2, mean)
  colnames(pooled_estimate) <- mycols
  rownames(pooled_estimate) <- myrows
  
  return(pooled_estimate)
  
}

############################################################
#
# Compare populations with Marshall 1/2 vs Marshall >2
# by checking for interaction effects
#
############################################################

check_for_interactions <- function(my_imp_list, mygroup){
  
  coef_bio <- vector(mode = "list", length = length(my_imp_list))
  coef_crash <- vector(mode = "list", length = length(my_imp_list))
  coef_impact <- vector(mode = "list", length = length(my_imp_list))
  
  #Within each imputed dataset
  
  for (i in 1:length(my_imp_list)) {
    # fit 3 models using biomarkers only / crash / impact variables

    
    fit_crash <- 
      my_imp_list[[i]] %>%
      with(glm(Outcome_5plus ~ 
                 Age_crash + 
                 Latest_GCS + 
                 Reactive_pupils +
                 MajorEI +
                 CT_TAI +
                 CT_SAH +
                 CT_Mids +
                 CT_Haematoma,
               family = "binomial"))
    
    fit_impact <- 
      my_imp_list[[i]] %>%
      with(glm(Outcome_5plus ~  
                 Age +
                 Motor_score +
                 Reactive_pupils +
                 Hypoxia +
                 Hypotension +
                 Marshall_score +
                 CT_SAH +
                 CT_EDH,
               family = "binomial"))
    
    # and use them to predict outcome
    
    my_imp_list[[i]]$pred_crash <- 
      predict(fit_crash, 
              newdata = my_imp_list[[i]])
    
    my_imp_list[[i]]$pred_impact <- 
      predict(fit_impact, 
              newdata = my_imp_list[[i]])
    
    # use predictions to check for interaction effects
    # i.e. whether bio/crash/impact work better in
    # one vs the other patient group
    fit_1 <- glm(Outcome_5plus ~ get(mygroup) * pred_bio, data = my_imp_list[[i]], family = "binomial")
    fit_2 <- glm(Outcome_5plus ~ get(mygroup) * pred_crash, data = my_imp_list[[i]], family = "binomial")
    fit_3 <- glm(Outcome_5plus ~ get(mygroup) * pred_impact, data = my_imp_list[[i]], family = "binomial")
    
    
    coef_bio[[i]]    <- summary(fit_1)$coefficients
    coef_crash[[i]]  <- summary(fit_2)$coefficients
    coef_impact[[i]] <- summary(fit_3)$coefficients
    
  }
  
  all_coef <- list(coef_bio,
                   coef_crash,
                   coef_impact)
  names(all_coef) <- c("bio", "crash", "impact")
  
  return(all_coef)
  
}

###########################################################
#
# Format pooled coefficients from interaction function
#
###########################################################

format_pooled_coefs <- function(pooled_coefs, model){
  
  temp <- format(round(pooled_coefs, 3), nsmall = 3)
  temp <- temp %>% as.data.frame()
  
  temp$Variable <- rownames(temp)
  temp$Model    <- model
  temp <- 
    temp %>%
    select(Model, Variable, everything())
  
  colnames(temp) <- gsub("Pr(>|z|)", "p-value", colnames(temp), fixed = TRUE)
  temp[, "p-value"] <- ifelse(as.numeric(temp[, "p-value"]) < 0.001, "<0.001", temp[, "p-value"])
  temp[, "p-value"] <- as.character(temp[, "p-value"])
  
  return(temp)
}

####################################
# Get the mean and se of pred_bio 
# stratified by outcome and group
# pooled across imputed datasets
###################################

summarise_pred_bio <- function(imp_list_all) {
  #summarise pred_bio score (on logit scale) in each imputed dataset
  mean_list <- vector(mode = "list", length = m)
  se_list   <- vector(mode = "list", length = m)
  n_list    <- vector(mode = "list", length = m)
  
  for (i in 1:m) {
    se_list[[i]] <-
      imp_list_all[[i]] %>%
      group_by(Group2, Outcome_5plus) %>%
      summarise(se = sd(pred_bio)/n()) %>%
      ungroup() %>%
      select(se)
    mean_list[[i]] <-
      imp_list_all[[i]] %>%
      group_by(Group2, Outcome_5plus) %>%
      summarise(mean = mean(pred_bio)) %>%
      ungroup() %>%
      select(mean)
    n_list[[i]] <-
      imp_list_all[[i]] %>%
      group_by(Group2, Outcome_5plus) %>%
      summarise(n = n())
  }
  
  ####Pool results via Rubin's rules###
  
  #pooled mean estimates
  pooled_estimate <- pool_my_coefs(mean_list, imp_list_all)
  pooled_n        <- pool_my_coefs(n_list, imp_list_all)
  
  #pooled se
  #first calculate within imputation variance
  squ_fun <- function(x){
    x^2
  }
  
  temp <- lapply(se_list, squ_fun)
  arr_within <- array(unlist(temp) , c(nrow(se_list[[1]]),ncol(se_list[[1]]),m) )
  within_variance <- apply(arr_within, 1:2, mean)
  
  #then calculate between imputation variance
  sub_fun <- function(x){
    x - pooled_estimate
  }
  
  delta_estimate <- lapply(mean_list, sub_fun)
  delta_estimate_squ <- lapply(delta_estimate, squ_fun)
  arr_between <- array(unlist(delta_estimate_squ) , c(nrow(se_list[[1]]),ncol(se_list[[1]]),m))
  
  spec_fun <- function(x){
    sum(x)/(m-1)
  }
  between_variance <- apply(arr_between, 1:2, spec_fun)
  
  #Then calculate total variance and pooled standard error
  total_variance <- within_variance + between_variance + I(between_variance/m)
  pooled_se      <- sqrt(total_variance)
  
  #calculate the limits of the 95% confidence interval
  pooled_ll <- pooled_estimate - I(1.96*pooled_se)
  pooled_ul <- pooled_estimate + I(1.96*pooled_se)
  
  #combine results into a table
  pooled_estimate <- format(round(pooled_estimate, 3), nsmall = 3)
  pooled_ll <- format(round(pooled_ll, 3), nsmall = 3)
  pooled_ul <- format(round(pooled_ul, 3), nsmall = 3)
  
  res <- n_list[[1]]
  res$n <- round(res$n, 0)
  res$Bio <- paste0(pooled_estimate, " (", pooled_ll, "-", pooled_ul, ")")
  colnames(res) <- c("Marshall score", "Outcome", "N", "Biomarker score")
  res$`Marshall score` <- gsub("Marshall", "", res$`Marshall score`)
  
  ftpred_bio <-
    res %>%
    regulartable()%>%
    vline_left() %>%
    vline_right() %>%
    hline(i = 2) %>%
    merge_v(j = 1) %>%
    fix_border_issues(part = "all") %>%
    bold(part = "header")
  
  #save to word document
  save_as_docx(ftpred_bio, path = "Figures/Biomarker_score.docx")
}


###########################################################
#
# Adds a column that contains the probability of a poor
# outcome based on the lasso biomarker coefficients
#
###########################################################

get_bio_prob <- function(data, panel){
  
  # standardize biomarker concentrations
  # to allow usage of standardized coeefficients
  cols     <-  c("GFAP", "NFL", "NSE", "S100B", "Tau", "UCH.L1")
  stancols <-  paste0("stan_", cols)
  pre_proc_val     <- preProcess(data[,cols], method = c("center", "scale"))
  data[ ,stancols] <- predict(pre_proc_val, data[ ,cols])
  
  data$B_GFAP   <- panel["GFAP", 1]   * data$stan_GFAP
  data$B_NFL    <- panel["NFL", 1]    * data$stan_NFL
  data$B_NSE    <- panel["NSE", 1]    * data$stan_NSE
  data$B_S100B  <- panel["S100B", 1]  * data$stan_S100B
  data$B_Tau    <- panel["Tau", 1]    * data$stan_Tau
  data$B_UCH.L1 <- panel["UCH.L1", 1] * data$stan_UCH.L1
  
  #getting the logodds of a poor outcome
  bcols <- 
    data %>% 
    select(starts_with("B_")) %>% 
    colnames()
  data$pred_bio <- 
    rowSums(data[, bcols])
  data <-
    data %>%
    select(-contains("B_"))
  
  return(data)
  
}


#############################################################
#
# Fit the 8 models to the chosen dataset
#
#############################################################

fit_mymodels <- function(df, mymodels){
  
  models_imp <- vector(mode = "list", length = length(mymodels))
  names(models_imp) <- mymodels
  
  #midline shift and cisternal obstruction is extremely rare in my dataset
  #in some bootstrap samples there will be no patient in whom they will be present
  #in these cases I will need to omit the variable from the model
  
  if (df[1, "Label"] == "m1"){ #for true CT-occult cohort
    
    models_imp[["plain_refit_crash_5"]] <- 
      glm(Outcome_5plus ~ 
            Age_crash + 
            Latest_GCS + 
            Reactive_pupils +
            MajorEI, 
          data = df,
          family = "binomial")
    
    models_imp[["plain_refit_impact_5"]] <- 
      glm(Outcome_5plus ~ 
            Age +
            Motor_score +
            Reactive_pupils +
            Hypoxia +
            Hypotension, 
          data = df,
          family = "binomial")
    
    #model with biomarkers
    
    models_imp[["withbio_refit_crash_5"]] <- 
      glm(Outcome_5plus ~ 
            Age_crash + 
            Latest_GCS + 
            Reactive_pupils +
            MajorEI +
            pred_bio, 
          data = df,
          family = "binomial")
    
    models_imp[["withbio_refit_impact_5"]] <- 
      glm(Outcome_5plus ~ 
            Age +
            Motor_score +
            Reactive_pupils +
            Hypoxia +
            Hypotension +
            pred_bio, 
          data = df,
          family = "binomial")
  
  } else if (df[1, "Label"] == "m2" || #Marshall 2
             df[1, "Label"] == "m56"  || #Marshall 5 & 6
             df[1, "Label"] == "m34") { #Marshall 3 & 4 as too few within each subcategory
    
    if (summary(df$CT_Mids)[[2]] > 0 & summary(df$CT_Cistern)[[2]] > 0){
      
      #plain models
      
      models_imp[["plain_refit_crash_5"]] <- 
        glm(Outcome_5plus ~ 
              Age_crash + 
              Latest_GCS + 
              Reactive_pupils +
              MajorEI +
              CT_TAI +
              CT_Cistern +
              CT_SAH +
              CT_Mids +
              CT_Haematoma, 
            data = df,
            family = "binomial")
      
      models_imp[["plain_refit_impact_5"]] <- 
        glm(Outcome_5plus ~ 
              Age +
              Motor_score +
              Reactive_pupils +
              Hypoxia +
              Hypotension +
              CT_SAH +
              CT_EDH, 
            data = df,
            family = "binomial")
      
      #model with biomarkers
      
      models_imp[["withbio_refit_crash_5"]] <- 
        glm(Outcome_5plus ~ 
              Age_crash + 
              Latest_GCS + 
              Reactive_pupils +
              MajorEI +
              CT_TAI +
              CT_Cistern +
              CT_SAH +
              CT_Mids +
              CT_Haematoma +
              pred_bio, 
            data = df,
            family = "binomial")
      
      models_imp[["withbio_refit_impact_5"]] <- 
        glm(Outcome_5plus ~ 
              Age +
              Motor_score +
              Reactive_pupils +
              Hypoxia +
              Hypotension +
              CT_SAH +
              CT_EDH +
              pred_bio, 
            data = df,
            family = "binomial")
      
      
    } else if (summary(df$CT_Mids)[[2]] > 0 & summary(df$CT_Cistern)[[2]] == 0){
      
      models_imp[["plain_refit_crash_5"]] <- 
        glm(Outcome_5plus ~ 
              Age_crash + 
              Latest_GCS + 
              Reactive_pupils +
              MajorEI +
              CT_TAI +
              CT_SAH +
              CT_Mids +
              CT_Haematoma, 
            data = df,
            family = "binomial")
      
      models_imp[["plain_refit_impact_5"]] <- 
        glm(Outcome_5plus ~ 
              Age +
              Motor_score +
              Reactive_pupils +
              Hypoxia +
              Hypotension +
              CT_SAH +
              CT_EDH, 
            data = df,
            family = "binomial")
      
      #model with biomarkers
      
      models_imp[["withbio_refit_crash_5"]] <- 
        glm(Outcome_5plus ~ 
              Age_crash + 
              Latest_GCS + 
              Reactive_pupils +
              MajorEI +
              CT_TAI +
              CT_SAH +
              CT_Mids +
              CT_Haematoma +
              pred_bio, 
            data = df,
            family = "binomial")
      
      models_imp[["withbio_refit_impact_5"]] <- 
        glm(Outcome_5plus ~ 
              Age +
              Motor_score +
              Reactive_pupils +
              Hypoxia +
              Hypotension +
              CT_SAH +
              CT_EDH +
              pred_bio, 
            data = df,
            family = "binomial")
      
      
    } else if (summary(df$CT_Mids)[[2]] == 0 & summary(df$CT_Cistern)[[2]] > 0){
      
      models_imp[["plain_refit_crash_5"]] <- 
        glm(Outcome_5plus ~ 
              Age_crash + 
              Latest_GCS + 
              Reactive_pupils +
              MajorEI +
              CT_TAI +
              CT_Cistern +
              CT_SAH +
              CT_Haematoma, 
            data = df,
            family = "binomial")
      
      models_imp[["plain_refit_impact_5"]] <- 
        glm(Outcome_5plus ~ 
              Age +
              Motor_score +
              Reactive_pupils +
              Hypoxia +
              Hypotension +
              CT_SAH +
              CT_EDH, 
            data = df,
            family = "binomial")
      
      #model with biomarkers
      
      models_imp[["withbio_refit_crash_5"]] <- 
        glm(Outcome_5plus ~ 
              Age_crash + 
              Latest_GCS + 
              Reactive_pupils +
              MajorEI +
              CT_TAI +
              CT_Cistern +
              CT_SAH +
              CT_Haematoma +
              pred_bio, 
            data = df,
            family = "binomial")
      
      models_imp[["withbio_refit_impact_5"]] <- 
        glm(Outcome_5plus ~ 
              Age +
              Motor_score +
              Reactive_pupils +
              Hypoxia +
              Hypotension +
              CT_SAH +
              CT_EDH +
              pred_bio, 
            data = df,
            family = "binomial")
      
    } else {
      
      models_imp[["plain_refit_crash_5"]] <- 
        glm(Outcome_5plus ~ 
              Age_crash + 
              Latest_GCS + 
              Reactive_pupils +
              MajorEI +
              CT_TAI +
              CT_SAH +
              CT_Haematoma, 
            data = df,
            family = "binomial")
      
      models_imp[["plain_refit_impact_5"]] <- 
        glm(Outcome_5plus ~ 
              Age +
              Motor_score +
              Reactive_pupils +
              Hypoxia +
              Hypotension +
              CT_SAH +
              CT_EDH, 
            data = df,
            family = "binomial")
      
      #model with biomarkers
      
      models_imp[["withbio_refit_crash_5"]] <- 
        glm(Outcome_5plus ~ 
              Age_crash + 
              Latest_GCS + 
              Reactive_pupils +
              MajorEI +
              CT_TAI +
              CT_SAH +
              CT_Haematoma +
              pred_bio, 
            data = df,
            family = "binomial")
      
      models_imp[["withbio_refit_impact_5"]] <- 
        glm(Outcome_5plus ~ 
              Age +
              Motor_score +
              Reactive_pupils +
              Hypoxia +
              Hypotension +
              CT_SAH +
              CT_EDH +
              pred_bio, 
            data = df,
            family = "binomial")
      
    }
    } else {
      if (summary(df$CT_Mids)[[2]] > 0 & summary(df$CT_Cistern)[[2]] > 0){
        
        #plain models
        
        models_imp[["plain_refit_crash_5"]] <- 
          glm(Outcome_5plus ~ 
                Age_crash + 
                Latest_GCS + 
                Reactive_pupils +
                MajorEI +
                CT_TAI +
                CT_Cistern +
                CT_SAH +
                CT_Mids +
                CT_Haematoma, 
              data = df,
              family = "binomial")
        
        models_imp[["plain_refit_impact_5"]] <- 
          glm(Outcome_5plus ~ 
                Age +
                Motor_score +
                Reactive_pupils +
                Hypoxia +
                Hypotension +
                Marshall_score +
                CT_SAH +
                CT_EDH, 
              data = df,
              family = "binomial")
        
        #model with biomarkers
        
        models_imp[["withbio_refit_crash_5"]] <- 
          glm(Outcome_5plus ~ 
                Age_crash + 
                Latest_GCS + 
                Reactive_pupils +
                MajorEI +
                CT_TAI +
                CT_Cistern +
                CT_SAH +
                CT_Mids +
                CT_Haematoma +
                pred_bio, 
              data = df,
              family = "binomial")
        
        models_imp[["withbio_refit_impact_5"]] <- 
          glm(Outcome_5plus ~ 
                Age +
                Motor_score +
                Reactive_pupils +
                Hypoxia +
                Hypotension +
                Marshall_score +
                CT_SAH +
                CT_EDH +
                pred_bio, 
              data = df,
              family = "binomial")
        
        
      } else if (summary(df$CT_Mids)[[2]] > 0 & summary(df$CT_Cistern)[[2]] == 0){
        
        models_imp[["plain_refit_crash_5"]] <- 
          glm(Outcome_5plus ~ 
                Age_crash + 
                Latest_GCS + 
                Reactive_pupils +
                MajorEI +
                CT_TAI +
                CT_SAH +
                CT_Mids +
                CT_Haematoma, 
              data = df,
              family = "binomial")
        
        models_imp[["plain_refit_impact_5"]] <- 
          glm(Outcome_5plus ~ 
                Age +
                Motor_score +
                Reactive_pupils +
                Hypoxia +
                Hypotension +
                Marshall_score +
                CT_SAH +
                CT_EDH, 
              data = df,
              family = "binomial")
        
        #model with biomarkers
        
        models_imp[["withbio_refit_crash_5"]] <- 
          glm(Outcome_5plus ~ 
                Age_crash + 
                Latest_GCS + 
                Reactive_pupils +
                MajorEI +
                CT_TAI +
                CT_SAH +
                CT_Mids +
                CT_Haematoma +
                pred_bio, 
              data = df,
              family = "binomial")
        
        models_imp[["withbio_refit_impact_5"]] <- 
          glm(Outcome_5plus ~ 
                Age +
                Motor_score +
                Reactive_pupils +
                Hypoxia +
                Hypotension +
                Marshall_score +
                CT_SAH +
                CT_EDH +
                pred_bio, 
              data = df,
              family = "binomial")
        
        
      } else if (summary(df$CT_Mids)[[2]] == 0 & summary(df$CT_Cistern)[[2]] > 0){
        
        models_imp[["plain_refit_crash_5"]] <- 
          glm(Outcome_5plus ~ 
                Age_crash + 
                Latest_GCS + 
                Reactive_pupils +
                MajorEI +
                CT_TAI +
                CT_Cistern +
                CT_SAH +
                CT_Haematoma, 
              data = df,
              family = "binomial")
        
        models_imp[["plain_refit_impact_5"]] <- 
          glm(Outcome_5plus ~ 
                Age +
                Motor_score +
                Reactive_pupils +
                Hypoxia +
                Hypotension +
                Marshall_score +
                CT_SAH +
                CT_EDH, 
              data = df,
              family = "binomial")
        
        #model with biomarkers
        
        models_imp[["withbio_refit_crash_5"]] <- 
          glm(Outcome_5plus ~ 
                Age_crash + 
                Latest_GCS + 
                Reactive_pupils +
                MajorEI +
                CT_TAI +
                CT_Cistern +
                CT_SAH +
                CT_Haematoma +
                pred_bio, 
              data = df,
              family = "binomial")
        
        models_imp[["withbio_refit_impact_5"]] <- 
          glm(Outcome_5plus ~ 
                Age +
                Motor_score +
                Reactive_pupils +
                Hypoxia +
                Hypotension +
                Marshall_score +
                CT_SAH +
                CT_EDH +
                pred_bio, 
              data = df,
              family = "binomial")
        
      } else {
        
        models_imp[["plain_refit_crash_5"]] <- 
          glm(Outcome_5plus ~ 
                Age_crash + 
                Latest_GCS + 
                Reactive_pupils +
                MajorEI +
                CT_TAI +
                CT_SAH +
                CT_Haematoma, 
              data = df,
              family = "binomial")
        
        models_imp[["plain_refit_impact_5"]] <- 
          glm(Outcome_5plus ~ 
                Age +
                Motor_score +
                Reactive_pupils +
                Hypoxia +
                Hypotension +
                Marshall_score +
                CT_SAH +
                CT_EDH, 
              data = df,
              family = "binomial")
        
        #model with biomarkers
        
        models_imp[["withbio_refit_crash_5"]] <- 
          glm(Outcome_5plus ~ 
                Age_crash + 
                Latest_GCS + 
                Reactive_pupils +
                MajorEI +
                CT_TAI +
                CT_SAH +
                CT_Haematoma +
                pred_bio, 
              data = df,
              family = "binomial")
        
        models_imp[["withbio_refit_impact_5"]] <- 
          glm(Outcome_5plus ~ 
                Age +
                Motor_score +
                Reactive_pupils +
                Hypoxia +
                Hypotension +
                Marshall_score +
                CT_SAH +
                CT_EDH +
                pred_bio, 
              data = df,
              family = "binomial")
        
      }
    
  }
  
  
  return(models_imp)
  
}

##############################################################
#
# extract coefficients from a model object
#
##############################################################

coef_fun <- function(x){
  coefs <- summary(x)$coefficients
  return(coefs)
}

##############################################################
#
# for a chosen cohort, extract pooled model coefficients
# and save as a flextable in word
#
##############################################################

save_model_coefficients <- function(imp_list, imp_models_list, mytitle){
  
  #for each model, extract coefficients per imputed dataset
  #then pool coefficients across imputed datasets
  mymodels <- names(imp_models_list[[1]])
  pooled_model_coefs <- vector(mode = "list", length = length(mymodels))
  names(pooled_model_coefs) <- mymodels
  
  for (modnum in 1:length(mymodels)){
    certain_model_list <- 
      lapply(imp_models_list,'[[', modnum)
    certain_coef_list  <- 
      lapply(certain_model_list, coef_fun)
    pooled_model_coefs[[modnum]] <- 
      pool_my_coefs(certain_coef_list, imp_list)
  }
  
  #format pooled coefficients for each model
  for (modnum in 1:length(mymodels)){
    pooled_model_coefs[[modnum]] <- 
      format_pooled_coefs(pooled_coefs = pooled_model_coefs[[modnum]], 
                          model = names(pooled_model_coefs)[modnum])
  }
  
  #combine coefficients from multiple models into one table
  pooled_model_coefs <- 
    bind_rows(pooled_model_coefs)
  
  #format this table
  pooled_model_coefs$Model <- gsub("plain_refit_crash_5", "Refitted CRASH-CT", pooled_model_coefs$Model)
  pooled_model_coefs$Model <- gsub("plain_refit_impact_5", "Refitted IMPACT-CT", pooled_model_coefs$Model)
  pooled_model_coefs$Model <- gsub("withbio_refit_crash_5", "Refitted CRASH-CT plus biomarkers", pooled_model_coefs$Model)
  pooled_model_coefs$Model <- gsub("withbio_refit_impact_5", "Refitted IMPACT-CT plus biomarkers", pooled_model_coefs$Model)
  
  pooled_model_coefs$Variable <- gsub("Age", "Age (years)", pooled_model_coefs$Variable)
  pooled_model_coefs$Variable <- gsub("Age_crash", "Age over 40 (years)", pooled_model_coefs$Variable)
  pooled_model_coefs$Variable <- gsub("pred_bio", "Biomarker panel logit(probability)", pooled_model_coefs$Variable)
  pooled_model_coefs$Variable <- gsub("cprob", "CRASH-CT probability", pooled_model_coefs$Variable)
  pooled_model_coefs$Variable <- gsub("CT_Cisternpresent", "Cisternal compression (present)", pooled_model_coefs$Variable)
  pooled_model_coefs$Variable <- gsub("CT_EDHpresent", "Extradural haematoma (present)", pooled_model_coefs$Variable)
  pooled_model_coefs$Variable <- gsub("CT_Haematomapresent", "Unevacuated haematoma (present)", pooled_model_coefs$Variable)
  pooled_model_coefs$Variable <- gsub("CT_Midspresent", "Midline shift (present)", pooled_model_coefs$Variable)
  pooled_model_coefs$Variable <- gsub("CT_SAHpresent", "Subarachnoid haemorrhage (present)", pooled_model_coefs$Variable)
  pooled_model_coefs$Variable <- gsub("CT_TAIpresent", "Petechial haemorrhages (present)", pooled_model_coefs$Variable)
  pooled_model_coefs$Variable <- gsub("Hypotensionpresent or suspected", "Hypotension (present or suspected)", pooled_model_coefs$Variable)
  pooled_model_coefs$Variable <- gsub("Hypoxiapresent or suspected", "Hypoxia (present or suspected)", pooled_model_coefs$Variable)
  pooled_model_coefs$Variable <- gsub("iprob", "IMPACT-CT probability", pooled_model_coefs$Variable)
  pooled_model_coefs$Variable <- gsub("Latest_GCS", "Glasgow Coma Score", pooled_model_coefs$Variable)
  pooled_model_coefs$Variable <- gsub("MajorEIpresent", "Major extra-cranial injury (present)", pooled_model_coefs$Variable)
  pooled_model_coefs$Variable <- gsub("Marshall_score2", "Marshall score (2)", pooled_model_coefs$Variable)
  pooled_model_coefs$Variable <- gsub("Marshall_score3", "Marshall score (3)", pooled_model_coefs$Variable)
  pooled_model_coefs$Variable <- gsub("Marshall_score4", "Marshall score (4)", pooled_model_coefs$Variable)
  pooled_model_coefs$Variable <- gsub("Marshall_score5", "Marshall score (5)", pooled_model_coefs$Variable)
  pooled_model_coefs$Variable <- gsub("Marshall_score6", "Marshall score (6)", pooled_model_coefs$Variable)
  pooled_model_coefs$Variable <- gsub("Marshall_scoreunknown", "Marshall score (unknown)", pooled_model_coefs$Variable)
  pooled_model_coefs$Variable <- gsub("Motor_score2", "Motor score (2)", pooled_model_coefs$Variable)
  pooled_model_coefs$Variable <- gsub("Motor_score3", "Motor score (3)", pooled_model_coefs$Variable)
  pooled_model_coefs$Variable <- gsub("Motor_score4", "Motor score (4)", pooled_model_coefs$Variable)
  pooled_model_coefs$Variable <- gsub("Motor_score5_or_6", "Motor score (5 or 6)", pooled_model_coefs$Variable)
  pooled_model_coefs$Variable <- gsub("Reactive_pupils1", "Reactive pupils (1)", pooled_model_coefs$Variable)
  pooled_model_coefs$Variable <- gsub("Reactive_pupils2", "Reactive pupils (2)", pooled_model_coefs$Variable)
  pooled_model_coefs$Variable <- gsub("UCH.L1", "UCH-L1 (log pg/ml)", pooled_model_coefs$Variable)
  pooled_model_coefs$Variable <- gsub("NFL", "NFL (log pg/ml)", pooled_model_coefs$Variable)
  pooled_model_coefs$Variable <- gsub("Tau", "Tau (log pg/ml)", pooled_model_coefs$Variable)
  pooled_model_coefs$Variable <- gsub("GFAP", "GFAP (log ng/ml)", pooled_model_coefs$Variable)
  pooled_model_coefs$Variable <- gsub("NSE", "NSE (log ng/ml)", pooled_model_coefs$Variable)
  pooled_model_coefs$Variable <- gsub("S100B", "S100B (log ng/ml)", pooled_model_coefs$Variable)
  pooled_model_coefs$Variable <- gsub("(Intercept)", "Intercept", pooled_model_coefs$Variable, fixed = TRUE)
  
  #turn table into flextable object
  ft_model_coefs <- 
    pooled_model_coefs %>%
    regulartable() %>%
    vline_left() %>%
    vline_right() %>%
    vline(j = 1) %>%
    merge_v(j = 1) %>%
    fix_border_issues(part = "all") %>%
    bold(part = "header") %>%
    align(j = c(3,4,5,6), align = "right") %>%
    padding(padding = 1, padding.top = 3) %>%
    fontsize(size = 10, part = "all") %>%
    valign(valign = "bottom", j = c(2:6))
  
  #add horizontal line after merged cell
  row_loc <- rle(cumsum(ft_model_coefs$body$spans$columns[,1] ))$values
  ft_model_coefs <- 
    ft_model_coefs %>% 
    border(border.bottom = fp_border_default(),
           i=row_loc, 
           j = 1:6, 
           part="body") 
  ft_model_coefs <- 
    ft_model_coefs %>% 
    border(border.bottom = fp_border_default(),
           i = ft_model_coefs$body$spans$columns[,1] > 1, 
           j = 1, 
           part="body") %>% 
    border(border.bottom = fp_border_default(), 
           border.top = fp_border_default(),
           part = "header") 
  
  #save to word document
  save_as_docx(ft_model_coefs, path = paste0("Figures/modelcoefs_", mytitle ,".docx"))
  
}



##############################################################
#
# Test the chosen models on the chosen dataset
#
##############################################################

test_mymodels <- function(df, models_imp){
  
  #Make predictions on imp data
  
  prob_imp        <- vector(mode = "list", length = length(mymodels))
  names(prob_imp) <- mymodels
  
  for (j in 1:length(mymodels)){
    prob_imp[[j]]  <- predict(object = models_imp[[j]],
                              newdata = df,
                              type = "response")
  }
  
  
  #Fit calibration models
  
  cal_imp <- vector(mode = "list", length = length(mymodels))
  names(cal_imp) <- mymodels
  
  for (j in 1:length(mymodels)){
    
    y    <- df$Outcome_5plus
    prob <- prob_imp[[j]]
    cal_imp[[j]] <- 
      glm(y ~ log(prob/(1-prob)), family = "binomial")
  }
  
  
  #Now get all the performance measures I want to report
  # so that they can be corrected for optimism
  # AUC
  # Nahelkerke R2
  # calibration slope
  # calibration intercept
  # p-value LRT
  
  mycols <- c("auc", "hl", "cal_slope", "cal_int", "brier" ,"nag", "lrt", "delta_auc", "delta_nag")
  
  res <- matrix(, ncol = length(mycols), nrow = length(mymodels))
  colnames(res) <- mycols
  rownames(res) <- mymodels
  
  for (j in 1:length(mymodels)){
    
    #AUC
    res[j, "auc"] <- ci.auc(roc(df$Outcome_5plus ~ prob_imp[[j]], quiet = TRUE))[2]
    
    #Hosmer-Lemeshow test
    res[j, "hl"]  <- hoslem.test(x = df$Outcome_5plus, 
                                 y = prob_imp[[j]])$p.value
    
    #Calibration slope & intercept
    coef <- summary(cal_imp[[j]])$coefficients
    res[j, "cal_slope"]     <- coef[2, 1]
    res[j, "cal_int"]       <- coef[1, 1]
    
    #Brier score
    res[j, "brier"]     <- BrierScore(resp = df$Outcome_5plus,
                                      pred = prob_imp[[j]])
    
    #Nagelkerke R2
    nag           <- val.prob(p = prob_imp[[j]], 
                              y = df$Outcome_5plus, pl = FALSE)
    res[j, "nag"] <- nag["R2"][[1]]*100
    
    #LRT (comparing nested models)
    c <- 0.0000000001 # constant in case AUC = 1, logit of which would be Inf
    if (j <= 2){
      res[j, "lrt"] <- NA
      res[j, "delta_auc"] <- NA
      res[j, "delta_nag"] <- NA
    } else {
      temp          <- lmtest::lrtest(models_imp[[j-2]], models_imp[[j]])
      res[j, "lrt"] <- temp[2, "Pr(>Chisq)"]
      res[j, "delta_auc"] <- logit(I(res[j, "auc"] - c)) - logit(I(res[j-2, "auc"] - c))
      res[j, "delta_nag"] <- res[j, "nag"] - res[j-2, "nag"]
    }
    
  }
  
  prob            <- bind_rows(prob_imp)
  res_list        <- list(res, prob)
  names(res_list) <- c("res", "prob")
  return(res_list)
  
}



#######################################################################
#
# Get results for one imputed dataset
# B specified the number of bootstrap samples
#
#######################################################################

results_for_one_imp <- function(imputed_set,mymodels, B){
  
  #Fit models to first imputed dataset and test them on it
  imp_models <- fit_mymodels(imputed_set, mymodels)
  imp_res    <- test_mymodels(imputed_set, imp_models)[["res"]]
  
  
  #For this imputed dataset create B bootstrap samples
  boot_list <- list()
  for (b in 1:B) {
    set.seed(b)
    boot_list[[b]]=imputed_set[sample(x       = rownames(imputed_set), 
                                      size    = nrow(imputed_set), 
                                      replace = TRUE),]
  }
  
  #For each bootstrap sample calculate the amount of optimism
  optimism_list   <- vector(mode = "list", length = B) #will store optimism
  orig_res_list   <- vector(mode = "list", length = B) #will store performance of bootstrap model in original imputed data
  orig_prob_list  <- rep(0, nrow(imputed_set)) #will store mean predictions per patient across all bootstrapped models
  boot_res_list   <- vector(mode = "list", length = B) #will store performance of bootstrap model on the same data it was fit (for R2)
  
  for (b in 1:B){
    #get performance of bootstrap model on original imputed dataset
    boot_mymodels       <- fit_mymodels(boot_list[[b]], mymodels)
    boot_res_list[[b]]  <- test_mymodels(boot_list[[b]], boot_mymodels)[["res"]]
    temp                <- test_mymodels(imputed_set, boot_mymodels)
    orig_res_list[[b]]  <- temp[["res"]]
    optimism_list[[b]]  <- boot_res_list[[b]] - orig_res_list[[b]]
    
    #use predictions of each bootstrap model on original imputed dataset
    #to update the mean prediction per patient that will be used for plotting
    existing_weight     <- (b-1)/b
    new_weight          <- 1/b
    orig_prob_list      <- (existing_weight * orig_prob_list) + (new_weight * temp[["prob"]])
  }
  
  #calculate mean optimism for current imputed dataset across B bootstrapped samples
  mycols <- c("auc","hl", "cal_slope", "cal_int", "brier", "nag", "lrt", "delta_auc", "delta_nag")
  arr_opt <- array( unlist(optimism_list) , c(length(mymodels),length(mycols),B) )
  mean_optimism <- apply(arr_opt, 1:2, mean)
  colnames(mean_optimism) <- mycols
  rownames(mean_optimism) <- mymodels
  
  #adjust initial model performance estimates for optimism
  #but do not adjust R2 as already adjusted by rms::val.prob
  adj_res <- imp_res - mean_optimism
  adj_res[ , "nag"] <- imp_res[ , "nag"]
  adj_res[ , "delta_nag"] <- imp_res[ , "delta_nag"]
  
  
  #get standard errors across all bootstrap samples of this imputed dataset
  #as a conservative approach I am using the standard errors
  # of the models fitted on the bootstrap but tested on the imp data
  arr_orig <- array( unlist(orig_res_list) , c(length(mymodels),length(mycols),B) )
  standard_error <- function(x) sd(x) / sqrt(B)
  se_orig  <- apply(arr_orig, 1:2, standard_error)
  colnames(se_orig) <- mycols
  rownames(se_orig) <- mymodels
  
  #For R2 use standard error for models tested on same data it was fit on
  arr_boot <- array( unlist(boot_res_list) , c(length(mymodels),length(mycols),B) )
  se_boot  <- apply(arr_boot, 1:2, standard_error)
  colnames(se_boot) <- mycols
  rownames(se_boot) <- mymodels
  se_orig[ , "nag"]       <-  se_boot[ , "nag"]
  se_orig[ , "delta_nag"] <-  se_boot[ , "delta_nag"]

  
  #keep mean predictions of bootstrapped model to create plots later
  prob_df                <- suppressMessages(bind_cols(imputed_set[,c("Outcome_5plus")], t(orig_prob_list)))
  colnames(prob_df)      <- c("Outcome_5plus", mymodels)
  imp_result             <- list(adj_res, se_orig, prob_df, imp_models)
  names(imp_result)      <- c("estimate", "se", "prob", "imp_models")
  
  return(imp_result)
  
}

##################################################################################
#
# Use Rubin's rules to pool results
# across m imputed datasets
#
##################################################################################

pool_myresults <- function(estimate_list, se_list, m, mymodels){
  #Pool results across imputed datasets
  #using Rubin's rules
  
  #####For point estimates that is the same as taking the mean
  mycols <- c("auc", "hl", "cal_slope", "cal_int", "brier", "nag", "lrt", "delta_auc", "delta_nag")
  arr_estimate <- array(unlist(estimate_list) , c(length(mymodels),length(mycols),m) )
  pooled_estimate <- apply(arr_estimate, 1:2, mean)
  colnames(pooled_estimate) <- mycols
  rownames(pooled_estimate) <- mymodels
  
  #####For SE Rubin's rules are more complicated
  
  #first calculate within imputation variance
  squ_fun <- function(x){
    x^2
  }
  temp <- lapply(se_list, squ_fun)
  
  mycols <- c("auc", "hl", "cal_slope", "cal_int", "brier","nag", "lrt", "delta_auc", "delta_nag")
  arr_within <- array(unlist(temp) , c(length(mymodels),length(mycols),m) )
  within_variance <- apply(arr_within, 1:2, mean)
  colnames(within_variance) <- mycols
  rownames(within_variance) <- mymodels
  
  
  #then calculate between imputation variance
  sub_fun <- function(x){
    x - pooled_estimate
  }
  
  delta_estimate <- lapply(estimate_list, sub_fun)
  
  delta_estimate_squ <- lapply(delta_estimate, squ_fun)
  mycols <- c("auc", "hl", "cal_slope", "cal_int", "brier","nag", "lrt", "delta_auc", "delta_nag")
  arr_between <- array(unlist(delta_estimate_squ) , c(length(mymodels),length(mycols),m) )
  
  spec_fun <- function(x){
    sum(x)/(m-1)
  }
  between_variance <- apply(arr_between, 1:2, spec_fun)
  colnames(between_variance) <- mycols
  rownames(between_variance) <- mymodels
  
  
  #Then calculate total variance and pooled standard error
  total_variance <- within_variance + between_variance + I(between_variance/m)
  pooled_se      <- sqrt(total_variance)
  
  
  #calculate the limits of the 95% confidence interval
  pooled_ll <- pooled_estimate - I(1.96*pooled_se)
  pooled_ul <- pooled_estimate + I(1.96*pooled_se)
  
  #set sensible limits for confidence intervals of
  #auc (0-1)
  pooled_ll[,"auc"] <- ifelse(pooled_ll[,"auc"] < 0, 0, pooled_ll[,"auc"])
  pooled_ul[,"auc"] <- ifelse(pooled_ul[,"auc"] > 1, 1, pooled_ul[,"auc"])
  
  #collect results in a list
  pooled_results <- list(pooled_estimate, pooled_ll, pooled_ul, pooled_se)
  names(pooled_results) <- c("pooled_estimate", "pooled_ll", "pooled_ul", "pooled_se")
  
  return(pooled_results)
}



#####################################################
#
# Round figures in results table
#
#####################################################

round_fun <- function(df){
  #round p-values to 3 digits
  #auc, slope, int 2
  df[,c("hl", "lrt")] <- 
    format(round(df[,c("hl", "lrt")],3), 
           nsmall = 3)
  df[, c("auc", "cal_slope", "cal_int", "brier", "delta_auc")] <- 
    format(round(as.numeric(df[, c("auc", "cal_slope", "cal_int", "brier", "delta_auc")]), 2),
           nsmall = 2)
  df[, c("nag", "delta_nag")] <-
    format(round(as.numeric(df[, c("nag", "delta_nag")]), 0),
           nsmall = 0)
  df <- 
    df %>%
    trimws()
  return(df)
}

####################################################
#
# save myres as a pretty table in word
#
####################################################

save_result_as_pretty_table <- function(myres, mytitle) {
  #paste point estimates with corresponding 95% confidence interval into the same cell
  temp <- paste0(round_fun(myres[["pooled_estimate"]]), 
                 "\n(", 
                 round_fun(myres[["pooled_ll"]]), 
                 "-", 
                 round_fun(myres[["pooled_ul"]]), 
                 ")"
  )
  dim(temp) <- c(nrow(myres[["pooled_estimate"]]), ncol(myres[["pooled_estimate"]]))
  colnames(temp) <- colnames(myres[["pooled_estimate"]])
  rownames(temp) <- rownames(myres[["pooled_estimate"]])
  
  #format model name display
  temp <- 
    temp %>%
    as.data.frame() 
  
  temp$model <-
    row.names(temp)
  
  temp <- 
    temp %>%
    select(auc,delta_auc, nag, delta_nag, brier, cal_int, cal_slope, hl, lrt, model)
  
  #adjust p-value (note only one p per two models, so NAs present
  temp[c(1, 2), "lrt"] <- p.adjust(temp[c(1, 2), "lrt"], method = "fdr")
  
  
  #reorder rows
  temp$model <-
    factor(temp$model,
           levels = c("plain_refit_crash_5",
                      "withbio_refit_crash_5",
                      "plain_refit_impact_5",
                      "withbio_refit_impact_5"))
  temp <- 
    temp %>%
    arrange(model)
  
  #remove CI for H-L p-value as non-sensical
  temp <-
    temp %>%
    separate(col = "hl", 
             into = "hl",
             sep = "\n") %>%
    separate(col = "lrt", 
             into = "lrt",
             sep = "\n")
  
  #for p-values (HL and lrt) display as "<0.001" if appropriate
  temp$hl  <- ifelse(temp$hl < 0.001, "<0.001", temp$hl)
  temp$lrt <- ifelse(temp$lrt < 0.001, "<0.001", temp$lrt)
  temp$hl  <- ifelse(temp$hl > 0.999, ">0.999", temp$hl)
  
  temp$lrt <-
    ifelse(temp$lrt == "NA", "", temp$lrt)
  temp$delta_auc <-
    ifelse(temp$delta_auc == "NA\n(NA-NA)", "", temp$delta_auc)
  temp$delta_nag <-
    ifelse(temp$delta_nag == "NA\n(NA-NA)", "", temp$delta_nag)
  
  row.names(temp) <- temp$model
  temp <- temp %>% select(-model)
  temp <- t(temp) %>% as.data.frame()
  
  temp$Metric <- row.names(temp)
  temp <- 
    temp %>%
    select(Metric, everything())
  
  #select only those metrics I want to display
  temp <- 
    temp %>%
    filter(Metric %in% c("auc", "nag", "cal_int", "cal_slope", "lrt"))
  
  #change order and label of metrics displayed
  temp$Metric <-
    factor(temp$Metric,
           levels = c("auc", "delta_auc", "nag", "delta_nag", "brier", "cal_int", "cal_slope", "hl", "lrt"),
           labels = c("Area under\nthe curve", 
                      "Incremental AUC",
                      "Variation\nexplained (%)", 
                      "Incremental variation\nexplained (%)",
                      "Brier score", 
                      "Calibration\nintercept", 
                      "Calibration\nslope", 
                      "Hosmer-\nLemeshow", 
                      "Likelihood ratio\ntest p-value"))
  
  gose5 <-
    temp %>% 
    select(Metric, contains("_5"))
  
  
  #create flextable object for GOSE < 5 results
  ft5 <-
    gose5 %>%
    regulartable() %>%
    set_header_labels(Metric = "Proteins added",
                      "plain_refit_crash_5" = "none",
                      "withbio_refit_crash_5" = "panel",
                      "plain_refit_impact_5" = "none",
                      "withbio_refit_impact_5" = "panel") %>%
    add_header_row(values = c(" ","CRASH-CT", "IMPACT-CT"),
                   colwidths = c(1, 2, 2))%>%
    vline(j = c(1, 3, 5))%>%
    vline_left() %>%
    vline_right() %>%
    align(align = "right", part = "all")%>%
    align(align = "center", part = "header", i = c(1:2)) %>%
    align(align = "left", j = c(1), part = "all") %>% 
    padding(padding = 0, part = "all", padding.right = 3, padding.left = 3, padding.top = 3, padding.bottom = 3) %>%
    fontsize(size = 10, part = "all") %>%
    fix_border_issues(part = "all") 
  
  save_as_docx(ft5, path = paste0("Figures/", mytitle, "_performance.docx"))
  
  
}


#####################################################
#
# Save model coefficients to word document
#
#####################################################

save_coefs <- function(temp, mytitle){
  
  #temp is the list of pooled model estimates
  # as generated from summary(mice::pool())
  
  for (i in 1:length(temp)){
    temp[[i]]$model <- names(temp)[i]
  }
  
  temp <- 
    bind_rows(temp) %>%
    as.data.frame()
  
  summary(as.factor(temp$model))
  temp$model <- gsub("_5", " (GOSE<5)", temp$model)
  temp$model <- gsub("_crash", " CRASH-CT", temp$model)
  temp$model <- gsub("_impact", " IMPACT-CT", temp$model)
  temp$model <- gsub("_refit", " refitted,", temp$model)
  temp$model <- gsub("plain", "no proteins,", temp$model)
  temp$model <- gsub("withbio", "with proteins,", temp$model)
  
  ft <-
    temp[-c(1:2),] %>%
    regulartable() %>%
    merge_v(j = "model") %>%
    theme_box() %>%
    valign( j = 1, valign = "top", part = "all") %>%
    fix_border_issues(part = "all") 
  
  save_as_docx(ft, path = paste0("Figures/", mytitle, "_coefs.docx"))
  
}

###########################################################
#
# Make calibration plot
#
###########################################################

plot_cal <- function(mydata, predvar, myoutcome, mycolor){
  temp <- mutate(mydata, bin = ntile(get(predvar), 10)) %>% 
    # Bin prediction into 10ths
    group_by(bin) %>%
    mutate(n = n(), # Get ests and CIs
           bin_pred = mean(get(predvar)), 
           bin_prob = mean(as.numeric(get(myoutcome))), 
           se = sqrt((bin_prob * (1 - bin_prob)) / n), 
           ul = bin_prob + 1.96 * se, 
           ll = bin_prob - 1.96 * se) %>%
    ungroup()
  
  g1 <- ggplot(data = temp, aes(x = bin_pred, y = bin_prob, ymin = ll, ymax = ul)) +
    geom_pointrange(size = 0.5, color = "black") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    geom_abline() + # 45 degree line indicating perfect calibration
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed", 
                color = "black", formula = y~-1 + x) + 
    # straight line fit through estimates
    geom_smooth(data = temp,
                aes(x = get(predvar), y = as.numeric(get(myoutcome))), 
                color = mycolor, se = FALSE, method = "loess") + 
    # loess fit through estimates
    xlab("") +
    ylab("Observed Probability") +
    theme_minimal() 
  
  # The distribution plot  
 g2 <- 
   ggplot(data = temp, aes(x = get(predvar), y = after_stat(density))) +
   geom_histogram(fill = "black", bins = 30) +
   scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
   scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, by = 2)) +
   xlab("Predicted Probability") +
   ylab("%") +
   theme_minimal() +
   theme(panel.grid.minor = element_blank())
  
  # Combine them    
  myplot <- cowplot::plot_grid(g1, NULL, g2, align = "hv", nrow = 3, rel_heights = c(3/4, -0.06, 1/4))
  
  return(myplot)
}


###########################################################
#
# Make ROC plot
#
###########################################################

make_one_roc_plot <- function(roc_with, roc_without, color_with, color_without, auc_with, auc_without){
  roclist <- list("roc_with" =  roc_with, "roc_without" =  roc_without)
  
  myplot <- 
    ggroc(roclist, legacy.axes= T, size=0.6, aes = "colour") +
    theme_classic(base_size = 6) +
    theme(legend.justification=c(1,0), 
          legend.position=c(0.98,0.08), 
          legend.text=element_text(size=6),
          legend.title=element_text(size=6),
          legend.key.size = unit(0.6,"line"),
          plot.title = element_text(size = 9, hjust = 0.5),
          axis.title.x = element_text(size = 6),
          axis.title.y = element_text(size = 6),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6)) +
    geom_abline(intercept = 0, 
                slope = 1, 
                color = "darkgrey", 
                linetype = "dashed") +
    scale_colour_manual(name   = "AUC (95% CI)",
                        values = c("roc_with" =  color_with, "roc_without" = color_without),
                        labels = c(paste0("with proteins:", auc_with),
                                   paste0("without proteins:", auc_without))) 
  
  return(myplot)
}


###############################################################
# 
# Use prob_list to make a panel
# of calibration plots
#
###############################################################

make_calibration_plots <- function(prob_list, mytitle){
  arr_prob <- 
    array( unlist(prob_list) , 
           c(nrow(prob_list[[1]]), ncol(prob_list[[1]]), m) )
  myprobs  <- 
    apply(arr_prob, 1:2, mean) %>% 
    as.data.frame()
  colnames(myprobs) <- colnames(prob_list[[1]])
  
  #make individual plots and store them in a list
  calplot_list <- lapply(c(2:ncol(myprobs)), FUN = function(x) {
    
    predvar    <- names(myprobs)[x]
    myoutcome  <- "Outcome_5plus"
    
    if(str_detect(predvar, "crash") == TRUE){
      if(str_detect(predvar, "plain") == TRUE){
        mycolor = "orange"
      } else {
        mycolor = "red"
      }
    } else {
      if(str_detect(predvar, "plain") == TRUE){
        mycolor = "lightblue"
      } else {
        mycolor = "purple"
      }
    }
    
    plot_cal(mydata = myprobs,
             predvar = predvar,
             myoutcome = myoutcome,
             mycolor = mycolor)
  })
  
  names(calplot_list) <- colnames(myprobs)[-1]
  
  
  #arrange outcome 5 plots in a panel and save image
  c <- calplot_list[c(1, 3, 2, 4)]
  calplot_origrefit_5 <- 
    grid.arrange(
      arrangeGrob(c[[1]], top="Model without proteins",left= "CRASH-CT"), 
      arrangeGrob(c[[2]], top="Model with proteins",   left= ""),
      arrangeGrob(c[[3]], top="", left="IMPACT-CT"), 
      arrangeGrob(c[[4]], top = "", left= ""),
      ncol=2, nrow=2, 
      widths = c(4,4),
      heights = c(4,4))
  
  ggsave(filename = paste0("Figures/calplot_", mytitle, ".tiff"), 
         plot = calplot_origrefit_5,
         dpi = 300,
         width = 20,
         height = 20,
         units = "cm",
         limitsize = FALSE)
  
}


####################################################################
#
# Create panels of ROC curves and save as images
#
####################################################################

make_roc_plots <- function(prob_list, myres, mytitle) {
  
  myprobs <- bind_rows(prob_list)
  
  #Make individual ROC curves and store in a list
  roc_list <- vector(mode = "list", length = length(mymodels))
  names(roc_list) <- mymodels
  
  for (i in 1:length(mymodels)){
    roc_list[[i]] <- roc(unlist(myprobs[, "Outcome_5plus"]), 
                         unlist(myprobs[, mymodels[i]]))
  }
  
  aucs <- formatC(myres[["pooled_estimate"]][, "auc"], digits = 2, format = "f")
  lls  <- formatC(myres[["pooled_ll"]][, "auc"], digits = 2, format = "f")
  uls  <- formatC(myres[["pooled_ul"]][, "auc"], digits = 2, format = "f")
  aucs <- paste0(aucs, " (", lls, "-", uls, ")")
  
  
  rocplot_list <- vector(mode = "list", length = length(mymodels)/2)
  names(rocplot_list) <- mymodels[str_detect(mymodels, "plain")==FALSE]
  
  for (i in 1:length(rocplot_list)){
    
    if(str_detect(names(rocplot_list)[i], "crash") == TRUE){
      color_without = "orange"
      color_with    = "red"
    } else {
      color_without = "lightblue"
      color_with    = "purple"
    }
    
      with    <- i + 2
      without <- i
      rocplot_list[[i]] <- 
        make_one_roc_plot(
          roc_with = roc_list[[with]],
          roc_without = roc_list[[without]],
          color_with = color_with,
          color_without = color_without,
          auc_with = aucs[with],
          auc_without = aucs[without])
    
  }
  
  #arrange Outcome 5 plots in a panel and save image
  c <- rocplot_list[c(1:2)]
  rocplot_origrefit_5 <- 
    grid.arrange(
      arrangeGrob(c[[1]], top = "CRASH-CT"), 
      arrangeGrob(c[[2]], top = "IMPACT-CT"),
      ncol=2, nrow=1)
  
  ggsave(filename = paste0("Figures/rocplot_", mytitle, ".tiff"), 
         plot = rocplot_origrefit_5,
         dpi = 300,
         width = 16,
         height = 8,
         units = "cm")
  
}

#################################
# Helper functions for rounding
#################################

myfun3 <- function(x){
  format(round(x,3), nsmall = 3)
}
myfun0 <- function(x){
  format(round(x,0), nsmall = 0)
}
myfun2 <- function(x){
  format(round(x,2), nsmall = 2)
}


###########################################
#
# Compare incremental value of biomarkers
# between two patient groups
#
###########################################

compare_increments <- function(suffix_A, suffix_B) {
  
  imp_list_A <- get(paste0("imp_list_", suffix_A))
  imp_list_B <- get(paste0("imp_list_", suffix_B))
  myres_A    <- get(paste0("myres_", suffix_A))
  myres_B    <- get(paste0("myres_", suffix_B))
  
  #######
  # Compare incremental AUC between groups
  #######
  
  metric <- "delta_auc"
  
  #initialize results table
  chosen_models <- c("withbio_refit_crash_5",
                     "withbio_refit_impact_5")
  mycols <- c("n_A", "mean_A", "se_A", "n_B", "mean_B", "se_B", "mean_se", "df", "t", "p")
  df <- as.data.frame(matrix(,
                             nrow = length(chosen_models),
                             ncol = length(mycols)))
  colnames(df) <- mycols
  rownames(df) <- chosen_models
  
  #sample size and degrees of freedom
  df$n_A <- nrow(imp_list_A[[1]])
  df$n_B <- nrow(imp_list_B[[1]])
  df$df  <- df$n_A + df$n_B -2
  
  #pooled mean estimates of logit iAUC
  df$mean_A <- as.data.frame(myres_A$pooled_estimate)[chosen_models, metric]
  df$mean_B <- as.data.frame(myres_B$pooled_estimate)[chosen_models, metric]
  
  #pooled standard error of logit iAUC
  class(myres_A$pooled_se) <- "matrix"
  class(myres_B$pooled_se) <- "matrix"
  df$se_A <- as.data.frame(myres_A$pooled_se)[chosen_models, metric]
  df$se_B <- as.data.frame(myres_B$pooled_se)[chosen_models, metric]
  df$mean_se <- (df$se_A + df$se_B)/2
  
  #t and p values
  df$t <- (df$mean_A - df$mean_B)/df$mean_se
  df$p <- pt(df$t, df$df, lower.tail=FALSE)
  
  auc_df <- df
  
  auc_df$Model <- rownames(auc_df)
  auc_df$Model <- ifelse(auc_df$Model == "withbio_refit_crash_5", "CRASH-CT",
                                       ifelse(auc_df$Model == "withbio_refit_impact_5", "IMPACT-CT", "none"))
  
  auc_df <-
    auc_df %>%
    mutate(across(c("mean_A", "mean_B", "se_A", "se_B", "mean_se"), myfun2))
  
  ######
  #compare incremental Nagelkerke R2 between groups
  ######
  
  metric <- "delta_nag"
  
  #initialize results table
  chosen_models <- c("withbio_refit_crash_5",
                     "withbio_refit_impact_5")
  mycols <- c("n_A", "mean_A", "se_A", "n_B", "mean_B", "se_B", "mean_se", "df", "t", "p")
  df <- as.data.frame(matrix(,
                             nrow = length(chosen_models),
                             ncol = length(mycols)))
  colnames(df) <- mycols
  rownames(df) <- chosen_models
  
  #sample size and degrees of freedom
  df$n_A <- nrow(imp_list_A[[1]])
  df$n_B <- nrow(imp_list_B[[1]])
  df$df  <- df$n_A + df$n_B -2
  
  #pooled mean estimates of logit iAUC
  df$mean_A <- as.data.frame(myres_A$pooled_estimate)[chosen_models, metric]
  df$mean_B <- as.data.frame(myres_B$pooled_estimate)[chosen_models, metric]
  
  #pooled standard error of logit iAUC
  class(myres_A$pooled_se) <- "matrix"
  class(myres_B$pooled_se) <- "matrix"
  df$se_A <- as.data.frame(myres_A$pooled_se)[chosen_models, metric]
  df$se_B <- as.data.frame(myres_B$pooled_se)[chosen_models, metric]
  df$mean_se <- (df$se_A + df$se_B)/2
  
  #t and p values
  df$t <- (df$mean_A - df$mean_B)/df$mean_se
  df$p <- pt(df$t, df$df, lower.tail=FALSE)
  
  nag_df <- df
  nag_df$Model <- rownames(nag_df)
  
  nag_df$Model <- ifelse(nag_df$Model == "withbio_refit_crash_5", "CRASH-CT",
                                       ifelse(nag_df$Model == "withbio_refit_impact_5", "IMPACT-CT", "none"))
  nag_df <-
    nag_df %>%
    mutate(across(c("mean_A", "mean_B", "se_A", "se_B", "mean_se"), myfun0))
  
  #########
  #Combine AUC and NAG dataframes
  ########
  
  auc_df$Metric <- "delta(logit(AUC))"
  nag_df$Metric <- "delta(R^2)"
  both_df <- rbind(auc_df, nag_df) %>% as.data.frame()
  comp <- both_df
  
  #format table and save to word
  comp <- 
    comp %>%
    select(Metric, Model, everything())
  
  return(comp)
  
}

#####################################################
#
# Prepare table to create increment plot
# using 2 Marshall score categories
#
######################################################

table_for_increment_plot2 <- function(m3456,  metric) {
  
  if (metric == "auc") {
    vec <- c(1:2)
  } else if (metric == "nag"){
    vec <- c(3:4)
  } else {
    print("Please choose valid metric (auc or nag)")
  }
  
  mycols <- c("Marshall", "Model", "Mean", "LL", "UL", "P" )
  df <- as.data.frame(matrix(,nrow = 4, ncol = length(mycols)))
  names(df) <- mycols
  df[1:2, "Marshall"] <- "1-2"
  df[1:2, "Model"] <- m3456[vec, "Model"]
  df[1:2, "Mean"] <- as.numeric(m3456[vec, "mean_A"])
  df[1:2, "LL"] <-as.numeric(m3456[vec, "mean_A"]) - 1.96*as.numeric(m3456[vec, "se_A"])
  df[1:2, "UL"] <-as.numeric(m3456[vec, "mean_A"]) + 1.96*as.numeric(m3456[vec, "se_A"])
  
  df[3:4, "Marshall"] <- "3-6"
  df[3:4, "Model"] <-m3456[vec, "Model"]
  df[3:4, "Mean"] <- as.numeric(m3456[vec, "mean_B"])
  df[3:4, "LL"] <-as.numeric(m3456[vec, "mean_B"]) - 1.96*as.numeric(m3456[vec, "se_B"])
  df[3:4, "UL"] <-as.numeric(m3456[vec, "mean_B"]) + 1.96*as.numeric(m3456[vec, "se_B"])
  df[3:4, "P"] <-as.numeric(m3456[vec, "p"])
  
  df[-c(1:2), "FDR"] <- p.adjust(df[-c(1:2), "P"])
  df$Sig <- ifelse(df$FDR < 0.05, "significant", "ns")
  df$Label1 <- ifelse(df$FDR < 0.001, "< 0.001", paste0("= ", format(round(df$FDR, 3), nsmall = 3)))
  df$Label2 <- ifelse(df$FDR < 0.001, "***", 
                      ifelse(df$FDR < 0.01, "**",
                             ifelse(df$FDR < 0.05, "*", "ns")))
  df$group1 <- "1-2"
  df$group2 <- "3-6"
  
  return(df)
  
}


#####################################################
#
# Prepare table to create increment plot
# using 4 Marshall score categories
#
######################################################

table_for_increment_plot4 <- function(m2, m34, m56, metric) {
  
  if (metric == "auc") {
    vec <- c(1:2)
  } else if (metric == "nag"){
    vec <- c(3:4)
  } else {
    print("Please choose valid metric (auc or nag)")
  }
  
  mycols <- c("Marshall", "Model", "Mean", "LL", "UL", "P" )
  df <- as.data.frame(matrix(,nrow = 8, ncol = length(mycols)))
  names(df) <- mycols
  df[1:2, "Marshall"] <- "1"
  df[1:2, "Model"] <-m2[vec, "Model"]
  df[1:2, "Mean"] <- as.numeric(m2[vec, "mean_A"])
  df[1:2, "LL"] <-as.numeric(m2[vec, "mean_A"]) - 1.96*as.numeric(m2[vec, "se_A"])
  df[1:2, "UL"] <-as.numeric(m2[vec, "mean_A"]) + 1.96*as.numeric(m2[vec, "se_A"])
  
  df[3:4, "Marshall"] <- "2"
  df[3:4, "Model"] <-m2[vec, "Model"]
  df[3:4, "Mean"] <- as.numeric(m2[vec, "mean_B"])
  df[3:4, "LL"] <-as.numeric(m2[vec, "mean_B"]) - 1.96*as.numeric(m2[vec, "se_B"])
  df[3:4, "UL"] <-as.numeric(m2[vec, "mean_B"]) + 1.96*as.numeric(m2[vec, "se_B"])
  df[3:4, "P"] <-as.numeric(m2[vec, "p"])
  
  df[5:6, "Marshall"] <- "3-4"
  df[5:6, "Model"] <-m34[vec, "Model"]
  df[5:6, "Mean"] <- as.numeric(m34[vec, "mean_B"])
  df[5:6, "LL"] <-as.numeric(m34[vec, "mean_B"]) - 1.96*as.numeric(m34[vec, "se_B"])
  df[5:6, "UL"] <-as.numeric(m34[vec, "mean_B"]) + 1.96*as.numeric(m34[vec, "se_B"])
  df[5:6, "P"] <-as.numeric(m34[vec, "p"])
  
  df[7:8, "Marshall"] <- "5-6"
  df[7:8, "Model"] <-m56[vec, "Model"]
  df[7:8, "Mean"] <- as.numeric(m56[vec, "mean_B"])
  df[7:8, "LL"] <-as.numeric(m56[vec, "mean_B"]) - 1.96*as.numeric(m56[vec, "se_B"])
  df[7:8, "UL"] <-as.numeric(m56[vec, "mean_B"]) + 1.96*as.numeric(m56[vec, "se_B"])
  df[7:8, "P"] <-as.numeric(m56[vec, "p"])
  
  df[-c(1:2), "FDR"] <- p.adjust(df[-c(1:2), "P"])
  df$Sig <- ifelse(df$FDR < 0.05, "significant", "ns")
  df$Label1 <- ifelse(df$FDR < 0.001, "< 0.001", paste0("= ", format(round(df$FDR, 3), nsmall = 3)))
  df$Label2 <- ifelse(df$FDR < 0.001, "***", 
                      ifelse(df$FDR < 0.01, "**",
                             ifelse(df$FDR < 0.05, "*", "ns")))
  df$group1 <- c(rep(NA, 2), rep("1", 2), rep("2", 2), rep("3-4", 2))
  df$group2 <- c(rep(NA, 2), rep("2", 2), rep("3-4", 2), rep("5-6", 2))
  
  return(df)

}

############################################
#
# Helper function for increment plot
#
############################################

get_legend  <- function(myggplot)
{
  tmp       <- ggplot_gtable(ggplot_build(myggplot))
  leg       <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend    <- tmp$grobs[[leg]]
  return(legend)
}


##############################################
#
# Plot incremental AUC or R2 for
# different patient cohorts (Marshall scores)
#
##############################################

plot_increments <- function(df, metric, mywidth, ll, ul, ypos, label){
  
  if (metric == "auc"){
    mytitle <- "Incremental logit(AUC)"
  } else if (metric == "nag") {
    mytitle <- "Incremental R2"
  } else {
    print("Please provide valid metric (auc or nag)")
  }
  
  myplot <- 
    ggplot(df, aes(x=Model, y=Mean, ymin=LL, ymax=UL, col=Marshall, fill=Marshall)) + 
    geom_hline(yintercept=0, size = 1, colour = "grey80") +
    #specify position here
    geom_linerange(size=2,position=position_dodge(width = mywidth)) +
    #specify position here too
    geom_point(size=2, shape=21, colour="white", stroke = 0.5 ,position=position_dodge(width = mywidth)) +
    scale_y_continuous(limits=c(ll, ul),name = mytitle, oob=scales::rescale_none) +
    theme_bw() +
    theme(legend.position = "bottom") + 
    labs(fill = "Marshall score", color = "Marshall score")
  
  if (label == "full"){
    myplot <-
      myplot +
      ggpubr::stat_pvalue_manual(
        df[-c(1:2), ], 
        x = "Model", 
        y.position = ypos,
        label = "q {Label1} {Label2}",
        position = position_dodge(mywidth/1.4), 
        size = 3)
  } else if (label == "short") {
    myplot <-
      myplot +
      ggpubr::stat_pvalue_manual(
        df[-c(1:2), ], 
        x = "Model", 
        y.position = ypos,
        label = "Label2",
        position = position_dodge(mywidth/1.4), 
        size = 3)
  } else {
    print("Please provide valid label: full or short")
  }
    
  
  return(myplot)
  
}
