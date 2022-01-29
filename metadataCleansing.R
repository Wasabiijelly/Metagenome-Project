#########################
## project: CI final term in ABC
## Hyeyoon
#########################

## Set Working Env. ##
setRepositories(ind=1:8)
WD <- '/disk3/bihy/metagenomeData'
setwd(WD)

library(dplyr)
library(rlist)
library(janitor) # clean_names

## Load Data ##
metadata38 <- fread(paste0(WD,'/MGYS00005138/MGYS00005138_metadata.csv'))
metadata84 <- fread(paste0(WD,'/MGYS00005184/MGYS00005184_metadata.csv'))
metadata94 <- fread(paste0(WD,'/MGYS00005194/MGYS00005194_metadata.csv'))
metadata <- list(metadata38,metadata84,metadata94)
#summary(metadata38)

sraRunTable38 <- fread(paste0(WD,'/MGYS00005138/MGYS00005138_SraRunTable.txt'))
sraRunTable84 <- fread(paste0(WD,'/MGYS00005184/MGYS00005184_SraRunTable.txt'))
sraRunTable94 <- fread(paste0(WD,'/MGYS00005194/MGYS00005194_SraRunTable.txt'))
sraRunTable <- list(sraRunTable38, sraRunTable84, sraRunTable94)
#summary(sraRunTable)


## Cleansing - 1. Find unique key ##
#uniqKey <- intersect(metadata38$run_accession, sraRunTable38$Run)
#sraRunTable38_c1 <- sraRunTable38 %>% filter(Run %in% uniqKey)
for (i in 1:length(sraRunTable)){
  sraRunTable[[i]] <- sraRunTable[[i]] %>% 
    rename(run_accession=Run)
}
meta_merged <- list()
for (i in 1:length(sraRunTable)){
  temp <- merge(metadata[[i]],sraRunTable[[i]], all = T)
  temp <- temp %>% clean_names()
  meta_merged <- list.append(meta_merged, temp)
}

## Considering host species & host disease ##
meta_merged[[1]]$host %>% factor() %>% levels()   # all Homo sapiens
meta_merged[[1]]$host_disease %>% factor() %>% levels()  # "" exists!!!!!

meta_merged[[2]]$Organism %>% factor() %>% levels()  # all Homo sapiens
meta_merged[[2]] <- meta_merged[[2]] %>% 
  mutate(host_disease=factor(rep("Healthy",nrow(meta_merged[[2]]))))

meta_merged[[3]]$sample_species %>% factor() %>% levels()   # all Homo sapiens
meta_merged[[3]] <- meta_merged[[3]] %>% 
  mutate(host_disease=factor(rep("Chronic rhinosinusitis",nrow(meta_merged[[3]]))))

#meta_merged[[1]]$Host_disease[which(meta_merged[[1]]$Host_disease=="")] <- NA # 더 효율적으로 "" -> NA 하는 방법

intersect(meta_merged[[1]]$run_accession,meta_merged[[2]]$run_accession)   # no intersection between samples

## cleansing meta_merged ##
idx_proj <- 1:length(meta_merged)

# change blank into NA 
for (j in idx_proj){
  col_class <- matrix(lapply(meta_merged[[j]], function(x) paste(class(x), collapse = ',')))
  idx_posix <- which(str_detect(col_class,'POSIX'))
  
  for (i in 1:ncol(meta_merged[[j]])){
    if(i %in% idx_posix){
      blankCnt <- append(blankCnt,0)
      next
    }  # posix 는 제외
    idx_blank <- which(meta_merged[[j]][,..i] == "")
    meta_merged[[j]][idx_blank,i] <- NA
  }
  
}

## chechk NA ##
result <- list()
for (j in idx_proj){
  naCnt <- c()
  for (i in 1:ncol(meta_merged[[j]])){
    tmp <- sum(meta_merged[[j]][,..i] %>% is.na())    # i는 character로 읽힘
    naCnt <- append(naCnt, tmp)
  }
  result <- list.append(result,naCnt)
}

meta_merged[[1]] <- meta_merged[[1]][-which(is.na(meta_merged[[1]][,90])),] ## 90 = host diesease column
## chechk NA  one more##
idx_na <- which(result[[1]] == 1021)
print(colnames(meta_merged[[1]][,c(..idx_na)]))
meta_merged[[1]] <- meta_merged[[1]][,-..idx_na]
## chechk NA  one more##

for(j in 1:length(meta_merged)){
  print(paste0('--------------------------------------',j,'----------------------------------------'))
  idx_na <- which(result[[j]] != 0)
  print(colnames(meta_merged[[j]][,c(..idx_na)]))
}

## check same column in one project ##
same_ <- list()
for (j in idx_proj){
  fst <- c()
  snd <- c()
  col_class <- matrix(lapply(meta_merged[[j]], function(x) paste(class(x), collapse = ',')))
  idx_posix <- which(str_detect(col_class,'POSIX|Date|logical'))
  
  for (i in 1:(ncol(meta_merged[[j]])-1)){
    for(k in (i+1):ncol(meta_merged[[j]])){
      if(i %in% idx_posix | k %in% idx_posix){
        next
      } 
      if(sum(meta_merged[[j]][,..i] == meta_merged[[j]][, ..k], na.rm = T) == nrow(meta_merged[[j]])){
        print('same')
        fst <- append(fst,i)
        snd <- append(snd,k)
      }
    }
  }
  tmp <- data.frame(first=fst, second=snd)
  same_ <- list.append(same_,tmp)
}
same_

# remove same column
rm_featureName_1 <- colnames(meta_merged[[1]][,76])
rm_featureName_2 <- colnames(meta_merged[[2]][,c(20,22,23,24,32,35,39,47,49,50,57,62,
                                                 66,70,71,76,78,85,85,89,90,95,96,
                                                 97,100,101,104)])
rm_featureName_3 <- colnames(meta_merged[[3]][,c(50,53,74)])

meta_merged[[1]] <- meta_merged[[1]] %>% 
  select(-rm_featureName_1)
meta_merged[[2]] <- meta_merged[[2]] %>% 
  select(-rm_featureName_2)
meta_merged[[3]] <- meta_merged[[3]] %>% 
  select(-rm_featureName_3)
rm_featureName <- data.frame(project1=rm_featureName_1,
                             project2=rm_featureName_2,
                             project3=rm_featureName_3)
#write.csv(rm_featureName,'remove_featureName.csv')
## check column class mismatch ##
intersect_1_2 <- intersect(colnames(meta_merged[[1]]),colnames(meta_merged[[2]]))
tmp <- meta_merged[[1]] %>% 
  select(all_of(intersect_1_2))
tmp2 <- meta_merged[[2]] %>% 
  select(all_of(intersect_1_2))
tmp3 <- matrix(lapply(tmp, function(x) paste(class(x), collapse = ',')))
tmp4 <- matrix(lapply(tmp2, function(x) paste(class(x), collapse = ',')))
tmp5 <- data.frame(pjt1=tmp3, pjt2=tmp4) # 39,43,55
meta_merged[[1]] <- meta_merged[[1]] %>% 
  rename(collection_date_1 = collection_date)

sum(meta_merged[[2]]$library_name!='unspecified')
meta_merged[[2]] <- meta_merged[[2]] %>% 
  select(-c('library_name'))     # unspecified
meta_merged[[3]] <- meta_merged[[3]] %>% 
  select(-c('library_name'))     # unspecified
meta_merged[[1]] <- meta_merged[[1]] %>% 
  rename(host_age_1=host_age)

## Making huge metadata ##
temp <- merge(meta_merged[[1]], meta_merged[[2]], all = T)
meta_final <- merge(meta_merged[[3]], temp, by = intersect(colnames(meta_merged[[3]]),colnames(temp)),all=T)


write.csv(meta_final,'meta_final_v2.csv')









