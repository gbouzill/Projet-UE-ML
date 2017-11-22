#############
# Chargement des données
######
library(RSQLite)

n <- dbDriver("SQLite", max.con=25)
con <- dbConnect(n, dbname="data/diabete.db")

# Liste des tables disponibles
req_tab <- dbSendQuery(con, "SELECT name FROM sqlite_master WHERE type='table'")
tab <- fetch(req_tab, n = -1)
dbHasCompleted(req_tab)
dbClearResult(req_tab)

# Table patient
req_patient <- dbSendQuery(con, "select * from training_patient")
patient <- fetch(req_patient, n = -1)
dbHasCompleted(req_patient)
dbClearResult(req_patient)

# Table diagnostics
req_diagnosis <- dbSendQuery(con, "select * from training_diagnosis")
diagnosis <- fetch(req_diagnosis, n = -1)
dbHasCompleted(req_diagnosis)
dbClearResult(req_diagnosis)

# Table transcript (= sejours)
# Il peut y avoir une erreur qui n'impacte a priori pas les données
req_transcript <- dbSendQuery(con, "select * from training_transcript")
transcript <- fetch(req_transcript, n = -1)
dbHasCompleted(req_diagnosis)
dbClearResult(req_diagnosis)

# Table Médicaments
req_med <- dbSendQuery(con, "select * from training_medication")
med <- fetch(req_med, n = -1)
dbHasCompleted(req_med)
dbClearResult(req_med)

dbDisconnect(con)

##########
# Data preparation
#######
library(plyr)
library(caret)
library(dplyr)
library(tidyr)
library(FactoMineR)
library(pROC)

# Fonction pour réaliser un capping sur une variable
# Cf. Cours
cap <- function(x) {
  th <- quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
  x[x < th[1]] <- th[1]
  x[x > th[2]] <- th[2]
  x
}

# Nettoyage de la table patient
patient_clean <- patient %>%
  select(PatientGuid,dmIndicator,Gender,YearOfBirth) %>%
  mutate(
    dmIndicator = factor(dmIndicator, levels = c(0,1), labels = c("Non", "Oui")),
    Gender_M = Gender == "M") %>%
  select(-Gender)

# Premier partionnement en apprentissage / test
set.seed(12345)
train_id <- createDataPartition(patient_clean$dmIndicator, p = 0.8)
train_guid <- patient$PatientGuid[train_id$Resample1]
# train <- complete_df[train_id$Resample1, ]
# test <- complete_df[-train_id$Resample1, ]

# Effectif des patients dans la partition apprentissage
train_nb_patient <- length(train_guid)

# Nettoyage de la table transcript

# Transformation des valeurs nulles ou à 0 en manquantes
transcript_clean <- transcript %>%
  mutate_at(.vars = vars(Height, SystolicBP, DiastolicBP, RespiratoryRate, HeartRate, Temperature),
            .funs = funs(ifelse(. == "NULL", NA, .) %>% as.numeric))

transcript_clean$Height[transcript_clean$Height == 0] <- NA
transcript_clean$Weight[transcript_clean$Weight == 0] <- NA
transcript_clean$BMI[transcript_clean$BMI == 0] <- NA

# Suppression des valeurs aberrantes
transcript_clean$BMI <- cap(transcript_clean$BMI)
transcript_clean$Height <- cap(transcript_clean$Height)
transcript_clean$Weight <- cap(transcript_clean$Weight)
transcript_clean$SystolicBP <- cap(transcript_clean$SystolicBP)
transcript_clean$DiastolicBP <- cap(transcript_clean$DiastolicBP)
transcript_clean$RespiratoryRate <- cap(transcript_clean$RespiratoryRate)
transcript_clean$Temperature <- cap(transcript_clean$Temperature)

# On ne conserve que la valeur maximale de chaque variable pour chaque patient
# On supprime les variables avec trop peu de données
# On modifie les valeurs à -Inf qui sont apparues quand les patients n'avaient
# aucune donnée.
transcript_clean <- transcript_clean %>%
  group_by(PatientGuid) %>%
  summarise_at(
    .vars = vars(Height, Weight, BMI, SystolicBP, DiastolicBP, RespiratoryRate, HeartRate, Temperature),
    .funs = funs(max(., na.rm = TRUE))
  ) %>%
  select(-HeartRate, -RespiratoryRate, -Temperature) %>%
  mutate_at(
    .vars = vars(SystolicBP, DiastolicBP),
    .funs = funs(ifelse(. == "-Inf", NA, .))
  )

# Nettoyage de la table medicament
# On crée une nouvelle variable basée uniquement sur le premier
# terme du médicament
med$molecule <- med$MedicationName %>% sapply(function(i) {
  str_split(i, " ")[[1]][[1]]
})

# On ne conserve que les médicament apparaissant que chez aux moins 1% des patients
# On transforme la variable en booléen dès lors que le patient a eu au moins
# une fois la molécule
med_filt <- med %>%
  filter(PatientGuid %in% train_guid) %>%
  group_by(molecule) %>%
  mutate(
    nb = n() / train_nb_patient * 100
  ) %>%
  filter(nb >= 1) %>%
  group_by(PatientGuid, molecule) %>%
  summarise(nb = n() > 0)

# On met les données en format large
med_wide <- med_filt %>% spread(molecule, nb, fill = FALSE) %>%
  left_join(patient_clean %>% select(PatientGuid, dmIndicator), by = "PatientGuid") %>%
  ungroup()

# Préparation des données de diagnosic

# Comme pour la cim-10, les chaque caractère de code correspond à un niveau
# dans la hierarchie de la terminologie.
# On choisit de regrouper les codes selon leur troisème niveau
diagnosis <- diagnosis %>%
  mutate(short_code = str_sub(ICD9Code, 1, 3))

# Même principe que pour les médicaments, on ne garde que les codes fréquent
# et on ne re conserve que l'information de la présence du code au moins
# une fois chez les patients
diag_filt <- diagnosis %>%
  group_by(short_code) %>%
  mutate(
    nb = n() / train_nb_patient * 100
  ) %>%
  filter(nb >= 10) %>%
  group_by(PatientGuid, short_code) %>%
  summarise(nb = n() > 0)

# Mise en format large
diag_wide <- diag_filt %>% spread(short_code, nb, fill = FALSE) %>%
  left_join(patient_clean %>% select(PatientGuid, dmIndicator), by = "PatientGuid") %>%
  ungroup()

# Jointure de l'ensemble des table sous forme d'un unique jeu de données
complete_df <- patient_clean %>%
  left_join(med_wide, by = c("PatientGuid", "dmIndicator")) %>%
  left_join(diag_wide, by = c("PatientGuid", "dmIndicator")) %>%
  left_join(transcript_clean, by = c("PatientGuid"))

# Ces deux variables n'ont qu'une donnée manquante donc on l'impute par la médiane
complete_df$SystolicBP[is.na(complete_df$SystolicBP)] <- median(complete_df$SystolicBP, na.rm = T)
complete_df$DiastolicBP[is.na(complete_df$DiastolicBP)] <- median(complete_df$DiastolicBP, na.rm = T)

# On considère les données manquantes comment la non présence de données et
# on les code en FALSE.
# On standardise es variables continues
complete_df <- complete_df %>%
  mutate_if(
    .predicate = funs(any(is.na(.))),
    .funs = funs(ifelse(is.na(.), FALSE, .))
  ) %>%
  mutate_if(
    .predicate = funs(is.numeric(.) | is.integer(.)),
    .funs = funs(scale(.))
  )

# On crée des variables binaire 0/1 pour toutes les variables qualitatives sauf pour
# la variable à prédire (dmIndicator)
dd <- dummyVars(~.-PatientGuid, data = complete_df, fullRank = TRUE)
complete_df <- predict(dd, complete_df) %>% data.frame() %>%
  dplyr::rename(dmIndicator = dmIndicator.Oui) %>%
  mutate(dmIndicator = factor(dmIndicator, levels = c(0,1), labels = c("Non", "Oui")))

# On partition notre jeu de données avec le partitionnement défini plus haut.
train <- complete_df[train_id$Resample1, ]
test <- complete_df[-train_id$Resample1, ]

# On subdivise notre jeu de données d'apprentissage en deux partitions :
# Apprentissage2 et Validation
train_id2 <- createDataPartition(train$dmIndicator, p = 0.5)
train2 <- train[train_id2$Resample1, ]
valid <- train[-train_id2$Resample1, ]

#############
# Modélisation
#####

# 1. Non supervisé
# Clustering ascendant hiérarchique sur la première partition de test
train_to_clust <- train %>% select(-dmIndicator)
train_to_clust <- sapply(train_to_clust, scale)

# Peut mettre quelques secondes
dist_mat <- dist(train_to_clust, method = "euclidean")
clust <- hclust(dist_mat, method = "ward.D2")
plot(clust)

# On crée 3 cluster
groups <- cutree(clust, 3)

# Répartition des cas de diabète selon les groups
prop.table(table(groups, train$dmIndicator), 1)

# Parallélisation sous Windows
# library(doParallel)
# cl <- makeCluster(detectCores() - 1)
# registerDoParallel(cl)

# Paramétrage pour l'optimisation des modèles (cf. Caret)
#  - validation croisée adaptative 10 blocs pour éviter
#    de tester toutes les combinaisons de paramètres
#  - on se base sur les rappels précisions pour évaluer les modèles (prSummary)
ctrl <- trainControl(
  method = "adaptive_cv",
  number = 10,
  repeats = 1,
  classProbs = TRUE,
  summaryFunction = prSummary,
  adaptive = list(
    min = 5, alpha = 0.05, method = "gls", complete = TRUE
  ),
  verboseIter = TRUE,
  search = "random",
  allowParallel = FALSE
)

# On attribue des poids différents à nos individus selon qu'ils ont le diabète
# ou non
model_weights <- ifelse(train2$dmIndicator == "Non",
(1/table(train2$dmIndicator)[1]) * 0.5,
(1/table(train2$dmIndicator)[2]) * 0.5)

# premier modèle randomForest sur l'ensemble de la partition apprentissage 2.
# On fournit les poids, on standardise les variables et on supprimes les
# variables qui auraient une variance quasi nulle.
# Attention l'AUC est ici l'aire sous la courbe rappel/précision et non sous
# la courbe ROC.
fit_all <- train(
  dmIndicator ~ .,
  data = train2,
  weights = model_weights,
  method = "rf",
  preProc = c("center", "scale", "nzv"),
  metric = "AUC",
  trControl = ctrl
)

# Aires sous la courbe ROC à parition des trois partition à disposition
roc(train2$dmIndicator, predict(fit_all, train2, type = "prob")[,2])
roc(valid$dmIndicator, predict(fit_all, valid, type = "prob")[,2])
roc(test$dmIndicator, predict(fit_all, test, type = "prob")[,2])

# Evaluation des faux positifs et faux négatifs en fonction de nos clusters
valid_pred <- predict(fit_all, valid, type = "prob")[,2] > 0.5
fp <- valid$dmIndicator == "Non" & valid_pred
fn <- valid$dmIndicator == "Oui" & !valid_pred

valid_groups <- groups[-train_id2$Resample1]

prop.table(table(valid_groups, fp), 1)
prop.table(table(valid_groups, fn), 1)

# Modèle basé uniquement sur les patients du premier cluster
model_weights <- ifelse(
  train2$dmIndicator[groups[train_id2$Resample1] == 1] == "Non",
  (1/table(train2$dmIndicator[groups[train_id2$Resample1] == 1])[1]) * 0.5,
  (1/table(train2$dmIndicator[groups[train_id2$Resample1] == 1])[2]) * 0.5
)

fit_g1 <- train(
  dmIndicator ~ .,
  data = train2 %>% filter(groups[train_id2$Resample1] == 1),
  weights = model_weights,
  method = "rf",
  preProc = c("center", "scale", "nzv"),
  metric = "AUC",
  trControl = ctrl
)

roc(train2$dmIndicator, predict(fit_g1, train2, type = "prob")[,2])
roc(valid$dmIndicator, predict(fit_g1, valid, type = "prob")[,2])
roc(test$dmIndicator, predict(fit_g1, test, type = "prob")[,2])

valid_pred <- predict(fit_g1, valid, type ="prob")[,2] > 0.5
fp <- valid$dmIndicator == "Non" & valid_pred
fn <- valid$dmIndicator == "Oui" & !valid_pred

valid_groups <- groups[-train_id2$Resample1]

prop.table(table(valid_groups, fp), 1)
prop.table(table(valid_groups, fn), 1)

# Modèle sur le deuxième cluster
model_weights <- ifelse(
  train2$dmIndicator[groups[train_id2$Resample1] == 2] == "Non",
  (1/table(train2$dmIndicator[groups[train_id2$Resample1] == 2])[1]) * 0.5,
  (1/table(train2$dmIndicator[groups[train_id2$Resample1] == 2])[2]) * 0.5
)

fit_g2 <- train(
  dmIndicator ~ .,
  data = train2 %>% filter(groups[train_id2$Resample1] == 2),
  weights = model_weights,
  method = "rf",
  preProc = c("center", "scale", "nzv"),
  metric = "AUC",
  trControl = ctrl
)

roc(train2$dmIndicator, predict(fit_g2, train2, type = "prob")[,2])
roc(valid$dmIndicator, predict(fit_g2, valid, type = "prob")[,2])
roc(test$dmIndicator, predict(fit_g2, test, type = "prob")[,2])

valid_pred <- predict(fit_g2, valid, type ="prob")[,2] > 0.5
fp <- valid$dmIndicator == "Non" & valid_pred
fn <- valid$dmIndicator == "Oui" & !valid_pred

valid_groups <- groups[-train_id2$Resample1]

prop.table(table(valid_groups, fp), 1)
prop.table(table(valid_groups, fn), 1)

# Modèle sur le troisième cluster
model_weights <- ifelse(
  train2$dmIndicator[groups[train_id2$Resample1] == 3] == "Non",
  (1/table(train2$dmIndicator[groups[train_id2$Resample1] == 3])[1]) * 0.5,
  (1/table(train2$dmIndicator[groups[train_id2$Resample1] == 3])[2]) * 0.5
)

fit_g3 <- train(
  dmIndicator ~ .,
  data = train2 %>% filter(groups[train_id2$Resample1] == 3),
  weights = model_weights,
  method = "rf",
  preProc = c("center", "scale", "nzv"),
  metric = "AUC",
  trControl = ctrl
)

roc(train2$dmIndicator, predict(fit_g3, train2, type = "prob")[,2])
roc(valid$dmIndicator, predict(fit_g3, valid, type = "prob")[,2])
roc(test$dmIndicator, predict(fit_g3, test, type = "prob")[,2])

valid_pred <- predict(fit_g3, valid, type ="prob")[,2] > 0.5
fp <- valid$dmIndicator == "Non" & valid_pred
fn <- valid$dmIndicator == "Oui" & !valid_pred

valid_groups <- groups[-train_id2$Resample1]

prop.table(table(valid_groups, fp), 1)
prop.table(table(valid_groups, fn), 1)


# Evaluation de la corrélation entre nos différents modèles (cf. Caret)
results <- resamples(list(mod = fit_all, mod1 = fit_g1, mod2 = fit_g2, mod3 = fit_g3))
modelCor(results)

# On réalise les prédictions avec nos différents modèles sur les partitions
# validation et test.
st <- data.frame(
  dmIndicator = valid$dmIndicator,
  xall = predict(fit_all, valid, type ="prob")[,2],
  x1 = predict(fit_g1, valid, type = "prob")[,2],
  x2 = predict(fit_g2, valid, type = "prob")[,2],
  x3 = predict(fit_g3, valid, type = "prob")[,2]
)

st2 <- data.frame(
  dmIndicator = test$dmIndicator,
  xall = predict(fit_all, test, type ="prob")[,2],
  x1 = predict(fit_g1, test, type = "prob")[,2],
  x2 = predict(fit_g2, test, type = "prob")[,2],
  x3 = predict(fit_g3, test, type = "prob")[,2]
)

# On recalcule les poids pour notre partition, de validation
model_weights <- ifelse(valid$dmIndicator == "Non",
                        (1/table(valid$dmIndicator)[1]) * 0.5,
                        (1/table(valid$dmIndicator)[2]) * 0.5)

# On realise un modèle final à partir des probabilités des modèles précédents
fit_final <- glm(dmIndicator~x1+x2+x3, data = st, weights = model_weights, family = binomial)

stopCluster(cl)

# On évalue notre modèle final sur la partition test
def <- predict(fit_final, st2)
roc(test$dmIndicator, def)
table(test$dmIndicator, def > 0.5)

