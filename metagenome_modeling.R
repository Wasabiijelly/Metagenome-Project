#########################
## project: CI final term in ABC
## Hyeyoon
#########################

## Set Working Env. ##
setRepositories(ind=1:8)
WD <- '/disk3/bihy/metagenomeData'
setwd(WD)

library(caret)
library(Rtsne)
library(reshape)
library(kknn)
library(caretEnsemble)
library(ggplot2)
library(rlist)

## Load data ##
data_phylum_orig <- fread(paste0(WD,'/final_phylum.csv'),stringsAsFactors = FALSE)
rowN <- data_phylum_orig$Phylum
data_phylum <- data_phylum_orig[,-c(1,2)]
rownames(data_phylum) <- rowN

data_class_orig <- fread(paste0(WD,'/final_class.csv'))
rowN <- data_class_orig$Class
data_class <- data_class_orig[,-c(1,2)]
rownames(data_class) <- rowN

data_order_orig <- fread(paste0(WD,'/final_order.csv'),stringsAsFactors = FALSE)
rowN <- data_order_orig$Order
data_order <- data_order_orig[,-c(1,2)]
rownames(data_order) <- rowN

data_family_orig <- fread(paste0(WD,'/final_family.csv'))
rowN <- data_family_orig$Family
data_family <- data_family_orig[,-c(1,2)]
rownames(data_family) <- rowN

data_genus_orig <- fread(paste0(WD,'/final_genus.csv'))
rowN <- data_genus_orig$Genus
data_genus <- data_genus_orig[,-c(1,2)]
rownames(data_genus) <- rowN

meta_final <- fread(paste0(WD,'/meta_final_v2.csv'))

## phylum level PCA classification ##
# Matching sample name in metagenome data
dim(meta_final)
dim(data_phylum)
dim(data_class)
dim(data_order)
dim(data_family)
dim(data_genus)
tmp<- intersect(colnames(data_phylum),meta_final$run_accession)
tmp2<-intersect(colnames(data_phylum),meta_final$analysis_accession)
length(tmp2)+length(tmp)
ncol(data_class)
nrow(meta_final)
accession <- c(tmp,tmp2)

################################info file###################################
info <- meta_final %>%
  filter(run_accession %in% accession) %>% 
  select(c(run_accession,host_disease)) %>% 
  dplyr::rename(accession = run_accession)
tmp_info <- meta_final %>%
  filter(analysis_accession %in% accession) %>% 
  select(c(analysis_accession,host_disease)) %>% 
  dplyr::rename(accession=analysis_accession)

info <- rbind(info,tmp_info)
info$host_disease <- factor(info$host_disease)
################################################################
# select sample matching with metadata
data_phylum_inter <- data_phylum %>% 
  select(intersect(info$accession,colnames(data_phylum)))
data_class_inter <- data_class %>% 
  select(intersect(info$accession,colnames(data_class)))
data_order_inter <- data_order %>% 
  select(intersect(info$accession,colnames(data_order)))
data_family_inter <- data_family %>% 
  select(intersect(info$accession,colnames(data_family)))
data_genus_inter <- data_genus %>% 
  select(intersect(info$accession,colnames(data_genus)))

# transpose
data_phylum_inter_t <- data.frame(t(data_phylum_inter))
data_class_inter_t <- data.frame(t(data_class_inter))
data_order_inter_t <- data.frame(t(data_order_inter))
data_family_inter_t <- data.frame(t(data_family_inter))
data_genus_inter_t <- data.frame(t(data_genus_inter))

# zero variance
nzv_phylum <- nearZeroVar(data_phylum_inter_t, saveMetrics = T, allowParallel = T)
nzv_class <- nearZeroVar(data_class_inter_t, saveMetrics = T, allowParallel = T)
nzv_order <- nearZeroVar(data_order_inter_t, saveMetrics = T, allowParallel = T)
nzv_family <- nearZeroVar(data_family_inter_t, saveMetrics = T, allowParallel = T)
nzv_genus <- nearZeroVar(data_genus_inter_t, saveMetrics = T, allowParallel = T)

data_phylum_inter_t <- data_phylum_inter_t[,-which(nzv_phylum$zeroVar==T)]
data_class_inter_t <- data_class_inter_t[,-which(nzv_class$zeroVar==T)]
data_order_inter_t <- data_order_inter_t[,-which(nzv_order$zeroVar==T)]
data_family_inter_t <- data_family_inter_t[,-which(nzv_family$zeroVar==T)]
data_genus_inter_t <- data_genus_inter_t[,-which(nzv_genus$zeroVar==T)]

# Scaling - Maximum Absolute Scaling for sparse data
data_phylum_inter_t <- data.frame(apply(data_phylum_inter_t,2,function(x) x/max(abs(x))))
data_class_inter_t <- data.frame(apply(data_class_inter_t,2,function(x) x/max(abs(x))))
data_order_inter_t <- data.frame(apply(data_order_inter_t,2,function(x) x/max(abs(x))))
data_family_inter_t<- data.frame(apply(data_family_inter_t,2,function(x) x/max(abs(x))))
data_genus_inter_t <- data.frame(apply(data_genus_inter_t,2,function(x) x/max(abs(x))))

# Set column name
colnames(data_phylum_inter_t) <- rownames(data_phylum)[-which(nzv_phylum$zeroVar==T)]
colnames(data_class_inter_t) <- rownames(data_class)[-which(nzv_class$zeroVar==T)]
colnames(data_order_inter_t) <- rownames(data_order)[-which(nzv_order$zeroVar==T)]
colnames(data_family_inter_t) <- rownames(data_family)[-which(nzv_family$zeroVar==T)]
colnames(data_genus_inter_t) <- rownames(data_genus)[-which(nzv_genus$zeroVar==T)]

data_phylum_with_disease <- cbind(data_phylum_inter_t,host_disease=info$host_disease)
data_class_with_disease <- cbind(data_class_inter_t,host_disease=info$host_disease)
data_order_with_disease <- cbind(data_order_inter_t,host_disease=info$host_disease)
data_family_with_disease <- cbind(data_family_inter_t,host_disease=info$host_disease)
data_genus_with_disease <- cbind(data_genus_inter_t,host_disease=info$host_disease)

# PCA
results_p <- prcomp(data_phylum_inter_t, scale = F)
fviz_pca_ind(results_p, habillage=data_phylum_with_disease$host_disease, label='none')
# autoplot(results_p ,data = pca_p, colour = 'host_disease', label.size=3, 
#          shape='circle') + theme_classic() + ggtitle('Phylum PCA Result')
summary(results_p)
results_p %>% biplot(cex=.5)
########## pca ADDITIONAL ##########
idx_outlier <- which(rownames(data_phylum_inter_t) %in% c('ERR1598260','ERR1598312','ERR1598236'))
tmp2 <- data_phylum_inter_t[-idx_outlier,]
nzv_tmp2 <- nearZeroVar(tmp2, saveMetrics = T)
tmp2 <- tmp2[,-which(nzv_tmp2$zeroVar==T)]
results_p3 <- prcomp(tmp2, scale = TRUE)
pca_p3 <- cbind(tmp2,host_disease=info$host_disease[-idx_outlier])
autoplot(results_p3 ,data = pca_p3,colour = 'host_disease', label.size=3, 
         shape=F) + theme_classic() + ggtitle('Phylum PCA Result3')

tmp <- data.frame(results_p$x[,1:49])
results_p2 <- prcomp(tmp, scale=T)
pca_p2 <- cbind(tmp,host_disease=info$host_disease)
autoplot(results_p2 ,data = pca_p2,colour = 'host_disease', label.size=3, 
         shape=F) + theme_classic() + ggtitle('Phylum PCA Result2')
fviz(results_p3,"var")
#ggsave('phylumPcaResult.pdf')
############################################################
results_c <- prcomp(data_class_inter_t, scale = F)
pca_c <- cbind(data_class_inter_t,info$host_disease)
autoplot(results_c ,data = pca_c,colour = 'info$host_disease', label.size=3, 
         shape='circle') + theme_classic()+ ggtitle('Class PCA Result')
#ggsave('classPcaResult.pdf')

results_o <- prcomp(data_order_inter_t, scale = F)
pca_o <- cbind(data_order_inter_t,info$host_disease)
autoplot(results_o ,data = pca_o,colour = 'info$host_disease', label.size=3, 
         shape='circle') + theme_classic()+ ggtitle('Order PCA Result')
#ggsave('orderPcaResult.pdf')

results_f <- prcomp(data_family_inter_t, scale = F)
pca_f <- cbind(data_family_inter_t,info$host_disease)
autoplot(results_f ,data = pca_f,colour = 'info$host_disease', label.size=3, 
         shape='circle') + theme_classic()+ ggtitle('Family PCA Result')
#ggsave('familyPcaResult.pdf')

results_g <- prcomp(data_genus_inter_t, scale = F)
pca_g <- cbind(data_genus_inter_t,info$host_disease)
autoplot(results_g ,data = pca_g,colour = 'info$host_disease', label.size=3, 
         shape=F) + theme_classic()+ ggtitle('Genus PCA Result')
#ggsave('genusPcaResult.pdf')

pca_g <- pca_g %>% clean_names()
levels(pca_g$info_host_disease)
test_idx <- which(pca_g$info_host_disease != 'Healthy')
dim(data_genus_inter_t)
dim(pca_g)
test <- data_genus_inter_t[test_idx,]
pca_g2 <- pca_g[test_idx,]

nzv <- nearZeroVar(test,saveMetrics = T)
test <- test[,-which(nzv$zeroVar==T)]
pca_g2 <- pca_g2[,-which(nzv$zeroVar==T)]
results_g2 <- prcomp(test, scale = TRUE)

autoplot(results_g2 ,data = pca_g2,colour = 'info_host_disease', label.size=3, 
         shape='circle') + theme_classic()+ ggtitle('Genus - Chronic rhinosinusitis  vs Acute Respiratory Infection')
#ggsave('genus_crVSar.pdf')

################# TSNE #################
tsne_phylum <- Rtsne(data_phylum_inter_t, perplexity=30,check_duplicates=F)
tsne_class <- Rtsne(data_class_inter_t, perplexity=30,check_duplicates=F)
tsne_order <- Rtsne(data_order_inter_t, perplexity=30,check_duplicates=F)
tsne_family <- Rtsne(data_family_inter_t, perplexity=30,check_duplicates=F)
tsne_genus <- Rtsne(data_genus_inter_t, perplexity=30,check_duplicates=F)

g1 <- ggplot(tsne_phylum$Y, aes(x = tsne_phylum$Y[, 1], y = tsne_phylum$Y[,2], colour = info$host_disease)) + 
  geom_point(alpha = 0.3) + theme_bw() + ggtitle('Phylum') ;g1
#ggsave('tsne1.pdf', width = 1900, height = 1080, units = 'px')
g2 <- ggplot(tsne_class$Y, aes(x = tsne_class$Y[, 1], y = tsne_class$Y[,2], colour = info$host_disease)) + 
  geom_point(alpha = 0.3) + theme_bw() + ggtitle('Class');g2
#ggsave('tsne2.pdf', width = 1900, height = 1080, units = 'px')
g3 <- ggplot(tsne_order$Y, aes(x = tsne_order$Y[, 1], y = tsne_order$Y[,2], colour = info$host_disease)) + 
  geom_point(alpha = 0.3) + theme_bw() + ggtitle('Order');g3
#ggsave('tsne3.pdf', width = 1900, height = 1080, units = 'px')
g4 <- ggplot(tsne_family$Y, aes(x = tsne_family$Y[, 1], y = tsne_family$Y[,2], colour = info$host_disease)) + 
  geom_point(alpha = 0.3) + theme_bw() + ggtitle('Family');g4
#ggsave('tsne4.pdf', width = 1900, height = 1080, units = 'px')
g5 <- ggplot(tsne_genus$Y, aes(x = tsne_genus$Y[, 1], y = tsne_genus$Y[,2], colour = info$host_disease)) + 
  geom_point(alpha = 0.3) + theme_bw() + ggtitle('Genus');g5
#ggsave('tsne5.pdf', width = 1900, height = 1080, units = 'px')

# tsne_figure <- ggarrange(g1, g2, g3, g4, g5,
#                             labels = c("Pylum","Class","Order","Family","Genus"),
#                             ncol = 3, nrow = 2);tsne_figure

############### k means clustering Phylum without outlier ##############
test_phylum <- scale(tmp2)
k_result <- kmeans(test_phylum,3)
table(k_result$cluster)
d_can <-get_dist(test_phylum, method = 'canberra') 
sil_result <- silhouette(k_result$cluster, d_can)
fviz_silhouette(sil_result)
# genus #
genusOutler <- fread('genusOutlier.txt',header = F)
idx_outlier_g <- which(rownames(data_genus_inter_t) %in% genusOutler$V1)
data_genus_with_disease <- data_genus_with_disease[-idx_outlier_g,]
test_genus <- data_genus_inter_t[-idx_outlier_g,]
nzv_test_g <- nearZeroVar(test_genus,saveMetrics = T)
test_genus <- test_genus[,-which(nzv_test_g$zeroVar==T)]
test_genus <- scale(test_genus,scale = T,center = T)
k_result_g <- kmeans(test_genus,3)
d_can_g <-parDist(test_genus, method = 'canberra',threads = 8) 
sil_result_g <- silhouette(k_result_g$cluster, d_can_g)
fviz_silhouette(sil_result_g)

meta_final$host_age_1
################# Modeling ###############
ctrl_phylum <- trainControl(method="repeatedcv", number=10, repeats = 3,
                               savePredictions = 'final',
                               allowParallel=T,
                               index = createFolds(data_phylum_with_disease$host_disease, 10))
ctrl_class <- trainControl(method="repeatedcv", number=10, repeats = 3,
                            savePredictions = 'final',
                            allowParallel=T,
                            index = createFolds(data_class_with_disease$host_disease, 10))
ctrl_order <- trainControl(method="repeatedcv", number=10, repeats = 3,
                            savePredictions = 'final',
                            allowParallel=T,
                            index = createFolds(data_order_with_disease$host_disease, 10))
ctrl_family <- trainControl(method="repeatedcv", number=10, repeats = 3,
                            savePredictions = 'final',
                            allowParallel=T,
                            index = createFolds(data_family_with_disease$host_disease, 10))
ctrl_genus <- trainControl(method="repeatedcv", number=10, repeats = 3,
                            savePredictions = 'final',
                            allowParallel=T,
                            index = createFolds(data_genus_with_disease$host_disease, 10))

# plot(phylum_model_pam)
# plot(phylum_model_protoclass)
# confusion_5nn <- table(pred=fitted(phylum_model_knn),real=test_phylum$host_disease)
# algorithmList <- c('knn', 'LogitBoost', 'svmRadial', 'earth','svmRadialWeights','mda','protoclass','pam'.'dda')
modelTypes <- list(knn = caretModelSpec(method = 'knn',tuneGrid = expand.grid(k=2:4)),
                   LogitBoost = caretModelSpec(method = 'LogitBooost'),
                   svmRadial = caretModelSpec(method = 'svmRadial'),
                   earth = caretModelSpec(method='earth'),
                   svmRadialWeights=caretModelSpec(method = 'svmRadialWeights'))
models_phylum <- caretList(host_disease~., 
                           data = data_phylum_with_disease,
                           trControl = ctrl_phylum, 
                           tuneList = modelTypes)
models_class <- caretList(host_disease~., 
                          data = data_class_with_disease,
                          trControl = ctrl_class, 
                          tuneList = modelTypes)
models_order <- caretList(host_disease~., 
                          data = data_order_with_disease,
                          trControl = ctrl_order, 
                          tuneList = modelTypes)
models_family <- caretList(host_disease~., 
                           data = data_family_with_disease,
                           trControl = ctrl_family, 
                           tuneList = modelTypes)
models_genus <- caretList(host_disease~., 
                          data = data_genus_with_disease,
                          trControl = ctrl_genus, 
                          tuneList = modelTypes)

results_phylum <- resamples(models_phylum)
results_class <- resamples(models_class)
results_order <- resamples(models_order)
results_family <- resamples(models_family)
results_genus <- resamples(models_genus)
summary(results_phylum)
summary(results_class)
summary(results_order)
summary(results_family)
summary(results_genus)
scales <- list(x=list(relation="free"), y=list(relation="free"))
bwplot(results_phylum, scales=scales, main = 'Performance by models (Phylum)')
bwplot(results_class, scales=scales, main = 'Performance by models (Class)')
bwplot(results_order, scales=scales, main = 'Performance by models (Order)')
bwplot(results_family, scales=scales, main = 'Performance by models (Family)')
bwplot(results_genus, scales=scales, main = 'Performance by models (Genus)')

models_phylum$LogitBooost$resample %>% 
  filter(max(Accuracy) == Accuracy) %>% 
  select(Resample) %>% 
  as.character() ## fold6
models_class$LogitBooost$resample %>% 
  filter(max(Accuracy) == Accuracy) %>% 
  select(Resample) %>% 
  as.character() ## fold6
models_order$LogitBooost$resample %>% 
  filter(max(Accuracy) == Accuracy) %>% 
  select(Resample) %>% 
  as.character() ## fold4
models_family$LogitBooost$resample %>% 
  filter(max(Accuracy) == Accuracy) %>% 
  select(Resample) %>% 
  as.character() ## fold8
models_genus$LogitBooost$resample %>% 
  filter(max(Accuracy) == Accuracy) %>% 
  select(Resample) %>% 
  as.character() ## fold8

idx_test <- models_phylum$LogitBooost$call$trControl$index$Fold06
predict_phylum <- data.frame(predict(models_phylum, newdata=data_phylum_inter_t[idx_test,]))
which(is.na(predict_phylum$LogitBooost))
confusion_phylum <- confusionMatrix(data = predict_phylum$LogitBooost %>% as.factor(),
                                    reference = data_phylum_with_disease[idx_test,]$host_disease) # Or using table()
ari_info <- confusion_phylum$byClass %>% t() %>% as.data.frame() %>% clean_names() %>% 
  select(class_acute_respiratory_infection) %>% 
  filter(rownames(.) %in% c('Sensitivity','Specificity')) %>% 
  dplyr::rename(phylum_ari = class_acute_respiratory_infection)

idx_test <- models_class$LogitBooost$call$trControl$index$Fold06
predict_class <- data.frame(predict(models_class, newdata=data_class_inter_t[idx_test,]))
confusion_class <- confusionMatrix(data = predict_class$LogitBooost %>% as.factor(),
                                    reference = data_class_with_disease[idx_test,]$host_disease)
tmp <- confusion_class$byClass %>% t() %>% as.data.frame() %>% clean_names() %>% 
  select(class_acute_respiratory_infection) %>% 
  filter(rownames(.) %in% c('Sensitivity','Specificity')) %>% 
  dplyr::rename(class_ari = class_acute_respiratory_infection)
ari_info <- cbind(ari_info,tmp)

idx_test <- models_order$LogitBooost$call$trControl$index$Fold04
predict_order <- data.frame(predict(models_order, newdata=data_order_inter_t[idx_test,]))
confusion_order <- confusionMatrix(data = predict_order$LogitBooost %>% as.factor(),
                                    reference = data_order_with_disease[idx_test,]$host_disease)
tmp <- confusion_order$byClass %>% t() %>% as.data.frame() %>% clean_names() %>% 
  select(class_acute_respiratory_infection) %>% 
  filter(rownames(.) %in% c('Sensitivity','Specificity')) %>% 
  dplyr::rename(order_ari = class_acute_respiratory_infection)
ari_info <- cbind(ari_info,tmp)

idx_test <- models_family$LogitBooost$call$trControl$index$Fold08
predict_family <- data.frame(predict(models_family, newdata=data_family_inter_t[idx_test,]))
confusion_family <- confusionMatrix(data = predict_family$LogitBooost %>% as.factor(),
                                    reference = data_family_with_disease[idx_test,]$host_disease)
tmp <- confusion_family$byClass %>% t() %>% as.data.frame() %>% clean_names() %>% 
  select(class_acute_respiratory_infection) %>% 
  filter(rownames(.) %in% c('Sensitivity','Specificity')) %>% 
  dplyr::rename(family_ari = class_acute_respiratory_infection)
ari_info <- cbind(ari_info,tmp)

idx_test <- models_genus$LogitBooost$call$trControl$index$Fold08
predict_genus <- data.frame(predict(models_genus, newdata=data_genus_inter_t[idx_test,]))
confusion_genus <- confusionMatrix(data = predict_genus$LogitBooost %>% as.factor(),
                                    reference = data_genus_with_disease[idx_test,]$host_disease)
tmp <- confusion_genus$byClass %>% t() %>% as.data.frame() %>% clean_names() %>% 
  select(class_acute_respiratory_infection) %>% 
  filter(rownames(.) %in% c('Sensitivity','Specificity')) %>% 
  dplyr::rename(genus_ari = class_acute_respiratory_infection)
ari_info <- cbind(ari_info,tmp)
ari_info <- ari_info %>% t() %>% as.data.frame() %>% 
  mutate(classes = c('Phylum','Class','Order','Family','Genus'))

ggplot(ari_info, aes(x = classes, y = Sensitivity, fill=classes)) +
  geom_bar(stat = 'identity')+
  scale_fill_brewer(palette="YlGnBu")+
  theme_classic()
ggsave('sensitivity.pdf')

ggplot(ari_info, aes(x = classes, y = Specificity, fill=classes)) +
  geom_bar(stat = 'identity')+
  scale_fill_brewer(palette="YlGnBu")+
  theme_classic()
ggsave('specificity.pdf')

###################### Modeling Comparision ####################
models_phylum$knn$resample %>% 
  filter(max(Accuracy) == Accuracy) %>% 
  select(Resample) %>% 
  as.character()
idx_test <- models_phylum$knn$call$trControl$index$Fold06
confusion_phylum <- list()
for(i in 1:ncol(predict_phylum)){
  tmp <- confusionMatrix(data = predict_phylum[,i] %>% as.factor(),
                                      reference = data_phylum_with_disease[idx_test,]$host_disease) 
  confusion_phylum <- list.append(confusion_phylum, tmp)
  print(1)
}
predict_phylum[,1]
