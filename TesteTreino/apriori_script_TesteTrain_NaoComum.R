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

write.csv(mdr_data_new, "Apriori_TesteTreino_NaoComum/df_colSelecionado.csv")

df <- read.csv("Apriori_TesteTreino_NaoComum/df_colSelecionado.csv", header = TRUE, colClasses = c(rep('factor',8)))
mdr_data_new$X <- NULL

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
#################################################
freq_RegrasNaoComum <- function(df){
  require(plyr)
  qual<- c("LHSandRHS")
  countedCombos <- plyr::ddply(df, qual, nrow )
  
  names(countedCombos) <- c(qual, "Frequencia")
  
  countedCombos <- countedCombos[with(countedCombos, order(-Frequencia)),]
  write.csv(countedCombos, "Apriori_TesteTreino_NaoComum/NaoComum/RegrasNaoComum.csv")
  return(countedCombos)
}

freq_RegrasComum <- function(df){
  require(plyr)
  qual<- c("LHSandRHS")
  countedCombos <- plyr::ddply(df, qual, nrow )
  
  names(countedCombos) <- c(qual, "Frequencia")
  
  countedCombos <- countedCombos[with(countedCombos, order(-Frequencia)),]
  write.csv(countedCombos, "Apriori_TesteTreino_NaoComum/NaoComum/RegrasComum.csv")
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
  write.csv(Nao_coincedente, paste("Apriori_TesteTreino_NaoComum/naoCoincedente", grau, sep = "/"))
  return(Nao_coincedente)
}

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

##------------Separete Desfechos to Column-------------
splitDesfechos2ColumnFunction <- function(DF){
  uncomumDFSplit <- data.frame(do.call('rbind', strsplit(as.character(DF$LHSandRHS),' ',fixed=TRUE)))
  rmColLHSandRHS <- select(DF,-c("LHSandRHS"))
  splitDF <- data.frame(cbind(setNames(uncomumDFSplit["X1"], "LHSandRHS"), setNames(uncomumDFSplit["X3"], "Desfechos"), rmColLHSandRHS)) 
  return(splitDF)
}

##---------------Mutate DF into Matrix---------------------------
mutateDFFunction <- function(DFToBeMutated){
  newSelected_df <- select(DFToBeMutated, LHSandRHS)
  
  newMutateDF <- newSelected_df %>% dplyr::mutate(ID =  row_number(), flg = 1) %>% tidyr::separate_rows(LHSandRHS, sep = ",") %>% spread(LHSandRHS, flg)
  
  newMutateDF$a=NULL
  
  aSelected <- select(DFToBeMutated, Desfechos:ncol(DFToBeMutated))
  matrix_coin <-  data.frame(cbind(newMutateDF, aSelected)) 
  
  return(matrix_coin)
}

######################################
######################################


#magic_for(print, silent = TRUE,temp=FALSE) # call magic_for()
run = 60
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
  plot(valor,type = "o")
  
}


#-----------------------
#-------------------------
#NOTE:: Exporta o arquivo e separe o atributo "Desfecho" para outra coluna
write.csv(GeneralInfo, "Apriori_TesteTreino_NaoComum/coincedentes/GeneralInformacaoTesteTrain.csv")
coincedenteDesfechos <- splitDesfechos2ColumnFunction(coincedente)
write.csv(coincedente, "Apriori_TesteTreino_NaoComum/coincedentes/coincedente0.03&0.05.csv")

coincedentes <- read.csv(file="Apriori_TesteTreino_NaoComum/coincedentes/coincedente0.03&0.05.csv", header = TRUE, sep =",")

##-------------Coincedente Nao comum <99%--------
FrequenciaRegrasNaoComum <- freq_RegrasNaoComum(rbind_RegrasNaoComum)
FrequenciaRegrasNaoComumDesfechos <- splitDesfechos2ColumnFunction(FrequenciaRegrasNaoComum)

NaoComumDistinct <- BuscaBinario(rbind_RegrasNaoComum, FrequenciaRegrasNaoComum)
write.csv(NaoComumDistinct, "Apriori_TesteTreino_NaoComum/EliminaRedundancia/NaoCoincedente0.03&0.05.csv")


################################################################
################################################################
#--------Função que Elimina Redeundancia----

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
        
        if(k == selfLength && is.numeric(leftDF[i,"confidence"])){
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
write.csv(RedundanciaNaoComumSplit, "Apriori_TesteTreino_NaoComum/EliminaRedundancia/NaoCoincedenteElimRedundancia.csv")

#RedundanciaNaoComumMatix <- mutateDFFunction(RedundanciaNaoComumSplit)


###############################################################
################################################################



#---------------Coincedente 100% Comum----------------
#-----------------------Mutate DF into Matrix-----------------------------------------
#newDFResult <- splitDesfechos2ColumnFunction(coincedentes)

#write.csv(newDFResult, "Apriori_TesteTreino_NaoComum/coincedentes/newcoincedente0.03&0.05.csv")

#CoincedenteMutatedMatrix <- mutateDFFunction(newDFResult)

#desfechos0 <- subset(newDFResult, Desfechos == "encerra2=0", select=(1:ncol(newDFResult)))
#desfechos1 <- subset(newDFResult, Desfechos == "encerra2=1", select=(1:ncol(newDFResult)))
#desfechos2 <- subset(newDFResult, Desfechos == "encerra2=2", select=(1:ncol(newDFResult)))
#desfechos3 <- subset(newDFResult, Desfechos == "encerra2=3", select=(1:ncol(newDFResult)))

#write.csv(desfechos0, "Apriori_TesteTreino_NaoComum/arquivo/desf0.csv")
#write.csv(desfechos1, "Apriori_TesteTreino_NaoComum/arquivo/desf1.csv")
#write.csv(desfechos2, "Apriori_TesteTreino_NaoComum/arquivo/desf2.csv")
#write.csv(desfechos3, "Apriori_TesteTreino_NaoComum/arquivo/desf3.csv")





#MatrixDesfechos00 <- subset(CoincedenteMutatedMatrix, Desfechos == "encerra2=0", select=(1:ncol(CoincedenteMutatedMatrix)))
#MatrixDesfechos11 <- subset(CoincedenteMutatedMatrix, Desfechos == "encerra2=1", select=(1:ncol(CoincedenteMutatedMatrix)))
#MatrixDesfechos22 <- subset(CoincedenteMutatedMatrix, Desfechos == "encerra2=2", select=(1:ncol(CoincedenteMutatedMatrix)))
#MatrixDesfechos33 <- subset(CoincedenteMutatedMatrix, Desfechos == "encerra2=3", select=(1:ncol(CoincedenteMutatedMatrix)))


#write.csv(MatrixDesfechos00, "Apriori_TesteTreino_NaoComum/desfechos/desf0Cura_coincedente.csv")
#write.csv(MatrixDesfechos11, "Apriori_TesteTreino_NaoComum/desfechos/desf1Abandono_coincedente.csv")
#write.csv(MatrixDesfechos22, "Apriori_TesteTreino_NaoComum/desfechos/desf2Falencia_coincedente.csv")
#write.csv(MatrixDesfechos33, "Apriori_TesteTreino_NaoComum/desfechos/desf3Obito_coincedente.csv")


###################TODO
#Count each item and save into new df
#Plot a graph showing - ctrl + shift + C
# desfechoCura <- desfechos00 %>% select( -c(Desfechos:ncol(desfechos00)))
# desfechoAbandon <- desfechos11 %>% select( -c(Desfechos:ncol(desfechos11)))
# desfechoFalencia <- desfechos22 %>% select( -c(Desfechos:ncol(desfechos22)))
# desfechoObito <- desfechos33 %>% select( -c(Desfechos:ncol(desfechos33)))
# 
# matrizToFreqGraph <- function(matrixDF){
#   countAttr <- sapply(matrixDF, FUN=table)
#   countedItemsDF <- as.data.frame(t(as.data.frame(do.call(cbind, countAttr))))
#   countedItemsDF <- tibble::rownames_to_column(countedItemsDF, "items")
#   colnames(countedItemsDF)[2] = "frequencies"
#   decreOrdfreqDataSet <- countedItemsDF[order(countedItemsDF$frequencies, decreasing=TRUE),]
#   return(barplot(decreOrdfreqDataSet$frequencies[1:10], names.arg = decreOrdfreqDataSet$items[1:10]))
# }
# 
# matrizToFreqGraph(desfechoCura)
# matrizToFreqGraph(desfechoAbandon)
# matrizToFreqGraph(desfechoFalencia)
# matrizToFreqGraph(desfechoObito)




