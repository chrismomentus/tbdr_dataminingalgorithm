library(arules)
library(dplyr)
library(tidyverse)
library(stringr)
library("ggplot2")
library(tibble)
library(tidyr)

##############################################
#data praparations

mdr_data <- read.csv(file="data_MDR_TB_nova.csv", header = TRUE, sep =",")

mdr_data_new <- select(mdr_data,sexo2, raca3, fxet3, escol4, mumigual, hiv3, alcool, comdiab, comtabagis, comdrogas, ppl, pdres5, cavitacao, bilateral, tipores3, encerra2)
mdr_data_new <- na.omit(mdr_data_new, action = TRUE) # listwise deletion of missing

write.csv(mdr_data_new, "Apriori_Redundancia_NaoComum/df_colSelecionado.csv")

df <- read.csv("Apriori_Redundancia_NaoComum/df_colSelecionado.csv", header = TRUE, colClasses = c(rep('factor',8)))
mdr_data_new$X <- NULL

#66 e 80
## funÃ§Ã£o para criar samples (base Aleatoria) de tamnaho 66%
df_sample <- function(dt_frame1,seed_val){
  
  smp_size <- floor(0.66 * nrow(dt_frame1))
  ## set the seed to make your partition reproducible
  set.seed(seed_val)
  t_i <- sample(seq_len(nrow(dt_frame1)), size = smp_size)
  return(t_i)
}

###################

Apriori_func <- function(dt_frame, seed_apriori){
  trans_treino <- as(dt_frame, "transactions")
  set.seed(seed_apriori)
  reg_treino <- apriori(data = trans_treino, parameter = list(supp = 0.03, conf = 0.05), 
  appearance <- list(rhs=c("encerra2=0","encerra2=1", "encerra2=2","encerra2=3"), default = "lhs"), control = list(verbose=F))
  
  return(reg_treino)
}

###########################################################
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
#######################################
freq_RegrasNaoComum <- function(df){
  require(plyr)
  qual<- c("LHSandRHS")
  countedCombos <- plyr::ddply(df, qual, nrow )
  
  names(countedCombos) <- c(qual, "Frequencia")
  
  countedCombos <- countedCombos[with(countedCombos, order(-Frequencia)),]
  write.csv(countedCombos, "Apriori_Redundancia_NaoComum/NaoComum/RegrasNaoComum.csv")
  return(countedCombos)
}

freq_RegrasComum <- function(df){
  require(plyr)
  qual<- c("LHSandRHS")
  countedCombos <- plyr::ddply(df, qual, nrow )
  
  names(countedCombos) <- c(qual, "Frequencia")
  
  countedCombos <- countedCombos[with(countedCombos, order(-Frequencia)),]
  write.csv(countedCombos, "Apriori_Redundancia_NaoComum/NaoComum/RegrasComum.csv")
  return(TRUE)
}
BuscaBinario <- function(df1, df2){
  df1 <- df1[order(df1$LHSandRHS,decreasing = TRUE),]
  df1 <- df1 %>% distinct(LHSandRHS, .keep_all = TRUE)
  
  resultingDF <- data.frame(cbind(df1, Frequencia = df2$Frequencia))
  
  return(resultingDF)
}
regrasNaoCoincedente <- function(df1, df2, nomearq){
  grau = paste("NaoCoincedentegrau", nomearq,"arq.csv", sep = "-")
  
  Nao_coincedente <- bind_rows(
    anti_join(df1, df2, by='LHSandRHS'),
    anti_join(df2, df1, by='LHSandRHS'),
  )
  #Nrow_val<- nrow(Nao_coincedente)
  write.csv(Nao_coincedente, paste("Apriori_Redundancia_NaoComum/naoCoincedente", grau, sep = "/"))
  return(Nao_coincedente)
}

###############################################

function_valida <- function(df,seed_sample, seed_apriori, nomearq){
  
  df_sample1 <- df_sample(df, seed_sample)
  
  df_bootstrap <- df[df_sample1, ] #66% da base
  
  regrasApriori <- Apriori_func(df_bootstrap, seed_apriori)
 
  #regrasElimRedun <- regrasApriori[ !is.redundant( regrasApriori, measure="confidence" ) ]
  
  regrasPreparadas<- rule_prep(regrasApriori)

  #Pegar varios em cima de lift 1.1
  regrasPreparadas <- subset(regrasPreparadas, lift >= 1.1, select= c(LHSandRHS,count,lift, confidence, support))
  #Compara duas bases e retorna uma base valida
  return(regrasPreparadas)
}

##------------Separete Desfechos to Column-------------
splitDesfechos2ColumnFunction <- function(DF){
  uncomumDFSplit <- data.frame(do.call('rbind', strsplit(as.character(DF$LHSandRHS),' ',fixed=TRUE)))
  rmColLHSandRHS <- select(DF,-c("LHSandRHS"))
  splitDF <- data.frame(cbind(setNames(uncomumDFSplit["X1"], "LHSandRHS"), setNames(uncomumDFSplit["X3"], "Desfechos"), rmColLHSandRHS)) 
  return(splitDF)
}

###############################################
run = 100
valor <- data.frame( matrix(nrow = run, ncol = 2))

for (i in 1:1) {
  
  seedValue1 = as.numeric(format(Sys.Date(), format="%m"))
  seedValue2 = seedValue1*2
  
  coincedente = function_valida(df, seedValue1, seedValue2)
  
  GeneralInfo <- data.frame(runCount=i, ValorSeed_1=seedValue1, ValorSeed_2=seedValue2, Coincedente1=nrow(coincedente), Coincedente2=0, RegraComum=nrow(coincedente), RegraNaoComum=0)
  
  valor[i,1] <- i
  valor[i,2] <- nrow(coincedente)
  print(paste(valor[i,1:2]))
  
  rbind_RegrasNaoComum = data.frame()
  #rbind_RegrasComum = data.frame()
  
  for (k in 2:run) {
    
    seedValue1 = seedValue1+3
    seedValue2 = seedValue2+3
    Regras_valida1 <- coincedente
    Regras_valida2 <- function_valida(df, seedValue1, seedValue2)
    
    #Intersect Das k frames de Regras 
    
    coincedente = merge(Regras_valida2, coincedente, by.x="LHSandRHS", by.y = "LHSandRHS")
    coincedente_df = as_tibble(coincedente, .name_repair = "minimal")
    
    df_naocomum <- regrasNaoCoincedente(Regras_valida2,coincedente_df,k)
    
    GeneralInfo[nrow(GeneralInfo)+1,] <- list(k, seedValue1, seedValue2, nrow(Regras_valida1), nrow(Regras_valida2) ,nrow(coincedente), nrow(df_naocomum))
    
    rbind_RegrasNaoComum <- rbind(rbind_RegrasNaoComum, df_naocomum[,1:5])
    
    #rbind_RegrasComum <- rbind(rbind_RegrasComum, coincedente_df[,1:5])
    
    
    valor[k,1] <- k
    valor[k,2] <- nrow(coincedente)
    print(paste(valor[k,1:2]))
    
    if(valor[k,2] == valor[k-1,2])
      return(print("Criterio de parada encontrado"))
    
  }
  #plot(valor,type = "o")
  
}
######################

write.csv(GeneralInfo, "Apriori_Redundancia_NaoComum/coincedente/GeneralInformacaoTreino.csv")
#NOTE:: Exporta o arquivo e separe o atributo "Desfecho" para outra coluna
#coincedenteDesfechos <- splitDesfechos2ColumnFunction(coincedente)
write.csv(coincedente, "Apriori_Redundancia_NaoComum/coincedente/coincedente.csv")
#-----------------------
#-------------------------

coincedentes <- read.csv(file="Apriori_Redundancia_NaoComum/coincedente/coincedente.csv", header = TRUE, sep =",")

##-------------Coincedente Nao comum <99%--------
FrequenciaRegrasNaoComum <- freq_RegrasNaoComum(rbind_RegrasNaoComum)
FrequenciaRegrasNaoComumDesfechos <- splitDesfechos2ColumnFunction(FrequenciaRegrasNaoComum)

NaoComumDistinct <- BuscaBinario(rbind_RegrasNaoComum, FrequenciaRegrasNaoComum)
write.csv(NaoComumDistinct, "Apriori_Redundancia_NaoComum/EliminaRedundancia/NaoCoincedente0.03&0.05.csv")


################################################################
################################################################
#--------Função que Elimina Redeundancia de NAO-COMUM----

leftDF <- NaoComumDistinct
count2 <- 0
leftDF$confidence <- round(leftDF$confidence, 2)
i <- 1
while(i <= nrow(leftDF)) {
  selfDF_split <- unlist(strsplit(leftDF[i,"LHSandRHS"], "[, ]+"))
  selfLength <- length(selfDF_split)
  
  for (l in 1:nrow(leftDF)) {
    count2 <- l
    downDF_Split <- unlist(strsplit(leftDF[l,"LHSandRHS"], "[, ]+"))
    downLength <- length(downDF_Split)
    
    for (k in 1:selfLength) {
      
      if(selfDF_split[k] %in% downDF_Split){
        
        if(k == selfLength && !is.na(leftDF[i,"confidence"])){
          selfDF_conf <- leftDF[i,"confidence"]
          downDF_conf <- leftDF[l,"confidence"]
          
          if(selfDF_conf >= downDF_conf && selfLength < downLength){
            leftDF <- leftDF[-l,]
            
          } 
          
        } 
      } else {
        break
      }
    }
  }
  
  print(paste("Loop ", i, " of loop ", l," and row of", nrow(leftDF) ))
  i <- i+1
}

#write.csv(leftDF, "Apriori_TesteTreino_NaoComum/EliminaRedundancia/ElimRedundanciaNaoComum.csv")

RedundanciaNaoComumSplit <- splitDesfechos2ColumnFunction(leftDF)
write.csv(RedundanciaNaoComumSplit, "Apriori_Redundancia_NaoComum/EliminaRedundancia/NaoCoincedenteElimRedundancia.csv")
