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

otuTable38 <- fread(paste0(WD,'/MGYS00005138/MGYS00005138_SSU_OTUtable.csv'))
otuTable84 <- fread(paste0(WD,'/MGYS00005184/MGYS00005184_OTU_Table.tsv'))
otuTable94 <- fread(paste0(WD,'/MGYS00005194/MGYS00005194_OTU_Table.tsv'))
otuTable <- list(otuTable38, otuTable84, otuTable94)


for (i in 1:length(otuTable)){
  tmp <- otuTable[[i]] %>% 
    separate(SampleID,c('Superkingdom','Kingdom','Phylum','Class','Order','Family','Genus','Species'),';')
  otuTable[[i]] <- tmp
}

Phylum <- list()
for(i in 1:length(otuTable)){
  #tmp <- subset(otuTable[[i]]$Phylum,!is.na(Phylum),drop = FALSE)
  tmp <- subset(otuTable[[i]],!is.na(otuTable[[i]]$Phylum),drop = FALSE)
  Phylum <- list.append(Phylum,tmp)
}

Class <- list()
for(i in 1:length(otuTable)){
  tmp <- subset(otuTable[[i]],!is.na(otuTable[[i]]$Class),drop = FALSE)
  Class <- list.append(Class,tmp)
}

Order <- list()
for(i in 1:length(otuTable)){
  tmp <- subset(otuTable[[i]],!is.na(otuTable[[i]]$Order),drop = FALSE)
  Order <- list.append(Order,tmp)
}

Family <- list()
for(i in 1:length(otuTable)){
  tmp <- subset(otuTable[[i]],!is.na(otuTable[[i]]$Family),drop = FALSE)
  Family <- list.append(Family,tmp)
}

Genus <- list()
for(i in 1:length(otuTable)){
  tmp <- subset(otuTable[[i]],!is.na(otuTable[[i]]$Genus),drop = FALSE)
  Genus <- list.append(Genus,tmp)
}

#sum(is.na(Order[[3]]$Order))
#write.csv(Phylum[[1]],'test_phylum1.csv')
tmp <- merge(Phylum[[1]],Phylum[[2]],all = T)
final_phylum <- merge(tmp,Phylum[[3]],all = T)
final_phylum <- final_phylum %>% 
  mutate(across(everything(), .fns = ~replace_na(.,0))) %>% 
  select(-c(1,2,4,5,6,7,8))
final_phylum <- aggregate(.~Phylum, data=final_phylum, sum)
final_phylum <- final_phylum[-1,]

tmp <- merge(Class[[1]],Class[[2]],all = T)
final_class <- merge(tmp,Class[[3]],all = T)
final_class <- final_class %>% 
  mutate(across(everything(), .fns = ~replace_na(.,0))) %>% 
  select(-c(1,2,3,5,6,7,8))
final_class <- aggregate(.~Class, data=final_class, sum)
final_class <- final_class[-1,]

tmp <- merge(Order[[1]],Order[[2]],all = T)
final_order <- merge(tmp,Order[[3]],all = T)
final_order <- final_order %>% 
  mutate(across(everything(), .fns = ~replace_na(.,0))) %>% 
  select(-c(1,2,3,4,6,7,8))
final_order <- aggregate(.~Order, data=final_order, sum)
final_order <- final_order[-1,]

tmp <- merge(Family[[1]],Family[[2]],all = T)
final_family <- merge(tmp,Family[[3]],all = T)
final_family <- final_family %>% 
  mutate(across(everything(), .fns = ~replace_na(.,0))) %>% 
  select(-c(1,2,3,4,5,7,8))
final_family <- aggregate(.~Family, data=final_family, sum)
final_family <- final_family[-1,]

tmp <- merge(Genus[[1]],Genus[[2]],all = T)
final_genus <- merge(tmp,Genus[[3]],all = T)
final_genus <- final_genus %>% 
  mutate(across(everything(), .fns = ~replace_na(.,0))) %>% 
  select(-c(1,2,3,4,5,6,8))
final_genus <- aggregate(.~Genus, data=final_genus, sum)
final_genus <- final_genus[-1,]

## Store final OTU result ##
write.csv(final_phylum, file='final_phylum.csv')
write.csv(final_class, file='final_class.csv')
write.csv(final_order, file='final_order.csv')
write.csv(final_family, file='final_family.csv')
write.csv(final_genus, file='final_genus.csv')