#### READ IN DATA SETS and LOAD IN PACKAGES ####

db.1 <- read.csv("humvar.csv") 
db.2 <- read.csv("lab_variants.csv")
db.3 <- read.csv("variants_annotated.csv")

library(randomForest)
library(tidyverse)
library(pROC)



                                                  #### CREATE THE TRAINING AND TEST SETS FOR HUMVAR ####


humvar_shuffled <- db.1[sample(1:nrow(db.1)),] # shuffles the data set so it's random and unbiased

humvar_bound <- nrow(humvar_shuffled) * 0.75 # selects 75% of the data

humvar_train <- humvar_shuffled[1:humvar_bound,] # get first 75% of rows for training
humvar_test <- humvar_shuffled[(humvar_bound+1):nrow(humvar_shuffled),] # get remaining 25% of rows for test  




                                                            #### RUN RANDOM FOREST ####


humvar_train_rf <- randomForest(as.factor(labels) ~ ., data = humvar_train) 
humvar_train_rf


                                                      
                                                #### MAKE HUMVAR PREDICTION AND FIND ACC, TPR, FPR ####


humvar_pred <- predict(humvar_train_rf, humvar_test)
humvar_pred

humvar_cm <- table(humvar_pred,humvar_test$labels)

humvar_pred_accuracy <- (humvar_cm[1,1]+humvar_cm[2,2]) / (humvar_cm[1,1] + humvar_cm[1,2] + humvar_cm[2,1] + humvar_cm[2,2])
humvar_pred_accuracy 
humvar_pred_TPR <- (humvar_cm[2,2]) / (humvar_cm[2,1] + humvar_cm[2,2]) 
humvar_pred_TPR 
humvar_pred_FPR <- (humvar_cm[1,2]) / (humvar_cm[1,2] + humvar_cm[1,1]) 
humvar_pred_FPR 



                                            #### Join db.2 & db.3 to show ID, label, and phenotype ####


pheno_patho_joined <- inner_join(db.2, db.3, by = "ids")



                                      #### FILTER NEW TABLE TO ONLY SHOW EPILEPSY OR MUSCULAR CONDITIONS ####


epilepsy_no_na <- pheno_patho_joined %>% 
  filter(phenotype == "epilepsy") %>%
  na.omit() 

muscular_no_na <- pheno_patho_joined %>%
  filter(phenotype == "muscular_conditions") %>%
  na.omit()



                                          #### CREATE SHUFFLED SET FOR PHENOTYPE TABLE ####


epilepsy_shuffled <- epilepsy_no_na[sample(1:nrow(epilepsy_no_na)),] # shuffles the data sets so it's random and unbiased
muscular_shuffled <- muscular_no_na[sample(1:nrow(muscular_no_na)),]

epilepsy_pred <- predict(humvar_train_rf, epilepsy_shuffled)
muscular_pred <- predict(humvar_train_rf, muscular_shuffled)

epilepsy_cm <- table(epilepsy_pred,epilepsy_shuffled$labels)
muscular_cm <- table(muscular_pred,muscular_shuffled$labels)



                                              #### EPILEPSY & MUSCULAR ACC, TPR, FPR ####

## EPILEPSY ##

epilepsy_pred_accuracy <- (epilepsy_cm[1,1] + epilepsy_cm[2,2]) / (epilepsy_cm[1,1] + epilepsy_cm[1,2] + epilepsy_cm[2,1] + epilepsy_cm[2,2])
epilepsy_pred_TPR <- (epilepsy_cm[2,2]) / (epilepsy_cm[2,1] + epilepsy_cm[2,2]) epilepsy_pred_TPR 
epilepsy_pred_FPR <- (epilepsy_cm[1,2]) / (epilepsy_cm[1,2] + epilepsy_cm[1,1])  


## MUSCULAR CONDITIONS ##

muscular_pred_accuracy <- (muscular_cm[1,1] + muscular_cm[2,2]) / (muscular_cm[1,1] + muscular_cm[1,2] + muscular_cm[2,1] + muscular_cm[2,2]) 
muscular_pred_TPR <- (muscular_cm[2,2]) / (muscular_cm[2,1] + muscular_cm[2,2])
muscular_pred_FPR <- (muscular_cm[1,2]) / (muscular_cm[1,2] + muscular_cm[1,1]) 




                                        #### PLOTTING THE RECEIVER OPERATING CURVE #### 


epilepsy_probs <- as.data.frame(predict(humvar_train_rf, epilepsy_shuffled, type = "prob"))
muscular_probs <- as.data.frame(predict(humvar_train_rf, muscular_shuffled, type = "prob"))

epilepsy_num_lab <- epilepsy_shuffled %>% 
  mutate(num_label = case_when(
    labels == "Benign" ~ 1,
    labels == "Pathogenic" ~ 2
  ))

muscular_num_lab <- muscular_shuffled %>% 
  mutate(num_label = case_when(
    labels == "Benign" ~ 1,
    labels == "Pathogenic" ~ 2
  ))

epilepsy_labels = c(epilepsy_num_lab$num_label)
epilepsy_benign_predictor = c(epilepsy_probs$Benign)

muscular_labels = c(muscular_num_lab$num_label)
muscular_benign_predictor = c(muscular_probs$Benign)

epilepsy_roc <- roc(response = epilepsy_labels, predictor = epilepsy_benign_predictor)
muscular_roc <- roc(response = muscular_labels, predictor = muscular_benign_predictor)

plot(epilepsy_roc, col = "blue", main = "ROC Curves for Epilepsy & Muscular Conditions")
lines(muscular_roc, col = "red")
legend("bottomright", legend = c("Epilepsy", "Muscular Conditions"),
      col = c("blue", "red"), lwd = 2)
