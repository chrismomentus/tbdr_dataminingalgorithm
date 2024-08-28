install.packages("tbdr19prediction")
library(arules)
library(dplyr)
library(tidyverse)
library(stringr)
library("ggplot2")
library(tibble)

##############################################
#data praparations

mdr_data <- read.csv(file="data_MDR_TB_nova.csv", header = TRUE, sep =",")

mdr_data_new <- select(mdr_data,sexo2, raca3, fxet3, escol4, mumigual, hiv3, alcool, comdiab, comtabagis, comdrogas, ppl, pdres5, cavitacao, bilateral, tipores3, encerra2)
mdr_data_new <- na.omit(mdr_data_new, action = TRUE) # listwise deletion of missing

mdr_data_new$situa_ence <- as.factor(mdr_data_new$encerra2)
mdr_data_new$situa_ence2 <- as.numeric(mdr_data_new$situa_ence)
class(mdr_data_new$encerra2)

write.csv(mdr_data_new, "df_colSelecionado.csv")

df <- read.csv("df_colSelecionado.csv", header = TRUE, colClasses = c(rep('factor',8)))
mdr_data_new$X <- NULL

## função para criar samples (base Aleatoria) de tamnaho 66%
df_sample <- function(dt_frame1,seed_val){
  smp_size <- floor(0.66 * nrow(dt_frame1))
  ## set the seed to make your partition reproducible
  #set.seed(seed_val)
  t_i <- sample(seq_len(nrow(dt_frame1)), size = smp_size)
  return(t_i)
}

rule_prep <- function(rule_dt){
  h <- length(rule_dt)
  lab_rules <- labels(rule_dt)
  lhs_rules = unlist(strsplit(lab_rules,"=> "))[seq(1,2*h,by=2)]#extrai os items do lhs das regras
  rhs_rules = unlist(strsplit(lab_rules,"=> "))[seq(2,2*h,by=2)]#extrai os items do rhs das regras
  
  qual_rules <- cbind.data.frame(rule_dt@quality$count, rule_dt@quality$lift,
                                 rule_dt@quality$confidence,rule_dt@quality$support)
  
  res2 = gsub("{", "", lhs_rules, fixed = TRUE) #expressao regular pra eliminar "{LHS}"
  reg_lhs <-  as.data.frame( gsub("}", "", res2, fixed = TRUE))
  
  res3 = gsub("{", "", rhs_rules, fixed = TRUE)#expressao regular pra eliminar "{RHS}"
  reg_rhs <- as.data.frame(gsub("}", "", res3, fixed = TRUE))
  
  data_rules <- cbind.data.frame( reg_lhs, qual_rules, reg_rhs)
  data_rules <- data_rules[-c(1:109),]
  colnames(data_rules) <- c("Regras", "count", "lift","confidence", "support","Desfechos")
  data_rules$LHSandRHS = paste(data_rules$Regras, data_rules$Desfechos, sep = " ")
  return(data_rules)
}


###################
df_sample1 <- df_sample(df, seed_sample)

df_treino <- df[df_sample1, ] #66% para treinar 2/3
df_teste <- df[-df_sample1, ]#34 para testar 1/3

set.seed(101)
trans_treino <- as(df_treino, "transactions")
reg_treino <- apriori(data = trans_treino, parameter = list(supp = 0.03, conf = 0.05), 
appearance <- list(rhs=c("encerra2=0","encerra2=1", "encerra2=2","encerra2=3"), default = "lhs"), control = list(verbose=F))

df_treino_data2 <- rule_prep(reg_treino)

df_treino_data_lift <- subset(df_treino_data, lift >= 1.1, select= c(LHSandRHS,count,lift, confidence, support))
##################
#################

Apriori_func <- function(dt_frame, seed_apriori){
  
  set.seed(seed_apriori)
  trans_treino <- as(dt_frame, "transactions")
  reg_treino <- apriori(data = trans_treino, parameter = list(supp = 0.03, conf = 0.05), 
  appearance <- list(rhs=c("encerra2=0","encerra2=1", "encerra2=2","encerra2=3"), default = "lhs"), control = list(verbose=F))
  
  return(reg_treino)
}

###########################################################

###############################################

function_valida <- function(df,seed_sample, seed_apriori){
  df_sample1 <- df_sample(df, seed_sample)
  
  df_treino <- df[df_sample1, ] #66% para treinar 2/3
  df_teste <- df[-df_sample1, ]#34 para testar 1/3
  
  df_test_reg <- Apriori_func(df_teste, seed_apriori )
  df_treino_reg <- Apriori_func(df_treino, seed_apriori)
  
  df_test_data <- rule_prep(df_test_reg)
  df_treino_data <- rule_prep(df_treino_reg)
  
  #Pegar varios em cima de lift 1.1
  df_test_data <- subset(df_test_data, lift >= 1.1, select= c(LHSandRHS,count,lift, confidence, support))
  df_treino_data <- subset(df_treino_data, lift >= 1.1, select= c(LHSandRHS,count,lift, confidence, support))
  #Compara duas bases e retorna uma base valida
  return(as.data.frame( filter(df_test_data,df_test_data[,1] %in% df_treino_data[,1] )))
}


#magic_for(print, silent = TRUE,temp=FALSE) # call magic_for()
run = 30
valor <- data.frame( matrix(nrow = run, ncol = 2))

for (i in 1:1) {
  #rdnum <- sample(1:200, 2, replace=FALSE)
  seed1 = 110
  seed2 = 122
  
  coincedente = function_valida(df, seed1, seed2)
  
  valor[i,1] <- i
  valor[i,2] <- nrow(coincedente)
  print(paste(valor[i,1:2]))
  
  
  for (k in 2:run) {
    
    rdnum <- sample(1:300, 2, replace=FALSE)
    
    Regras_valida1 <- function_valida(df, seed1, seed2)
    
    #Intersect Das k frames de Regras 
    
    coincedente = merge(Regras_valida1, coincedente, by.x="LHSandRHS", by.y = "LHSandRHS")
    
    valor[k,1] <- k
    valor[k,2] <- nrow(coincedente)
    print(paste(valor[k,1:2]))
    
  }
  plot(valor,type = "o")
  
}
#NOTE:: Exporta o arquivo e separe o atributo "Desfecho" para outra coluna
write.csv(coincedente, "coincedente0.03&0.05.csv")
#-----------------------
#-------------------------

coincedentes <- read.csv(file="coincedente0.03&0.05.csv", header = TRUE, sep =",")

#newfile1 <- within(coincedente0.1, LHSandRHS<-data.frame(do.call('rbind', strsplit(as.character(LHSandRHS), ' ', fixed=TRUE))))
newfile1 <- data.frame(do.call('rbind', strsplit(as.character(coincedentes$LHSandRHS),' ',fixed=TRUE)))
newfile1$X2 = NULL
newfile2 <- select(coincedentes,-c(1,2))
newfile3 <- coincedentes$X
newfile <-  data.frame(cbind(newfile3, newfile1, newfile2)) 
names(newfile)[1] <- "ID"
names(newfile)[2] <- "LHSandRHS"
names(newfile)[3] <- "Desfechos"

write.csv(newfile, "newcoincedente0.03&0.05.csv")



desfechos0 <- subset(newfile, Desfechos == "encerra2=0", select=(ID:support.y.1))
desfechos1 <- subset(newfile, Desfechos == "encerra2=1", select=(ID:support.y.1))
desfechos2 <- subset(newfile, Desfechos == "encerra2=2", select=(ID:support.y.1))
desfechos3 <- subset(newfile, Desfechos == "encerra2=3", select=(ID:support.y.1))

write.csv(desfechos0, "desf0.csv")
write.csv(desfechos1, "desf1.csv")
write.csv(desfechos2, "desf2.csv")
write.csv(desfechos3, "desf3.csv")




