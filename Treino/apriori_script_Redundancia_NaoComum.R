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

write.csv(mdr_data_new, "Apriori_Redundancia_NaoComum/df_colSelecionado.csv")

df <- read.csv("Apriori_Redundancia_NaoComum/df_colSelecionado.csv", header = TRUE, colClasses = c(rep('factor',8)))
mdr_data_new$X <- NULL

#66 e 80
## função para criar samples (base Aleatoria) de tamnaho 66%
df_sample <- function(dt_frame1,seed_val){
  
  smp_size <- floor(0.66 * nrow(dt_frame1))
  ## set the seed to make your partition reproducible
  set.seed(seed_val)
  t_i <- sample(seq_len(nrow(dt_frame1)), size = smp_size)
  return(t_i)
}

###################

Apriori_func <- function(dt_frame, seed_apriori){
  
  set.seed(seed_apriori)
  trans_treino <- as(dt_frame, "transactions")
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
regrasNaoCoincedente <- function(df1, df2, nomearq){
  grau = paste("NaoCoincedentegrau", nomearq,"arq.csv", sep = "-")
  
  Nao_coincedente <- bind_rows(
    anti_join(df1, df2, by='LHSandRHS'),
    anti_join(df2, df1, by='LHSandRHS'),
  )
  Nrow_val<- nrow(Nao_coincedente)
  write.csv(Nao_coincedente, paste("Apriori_Redundancia_NaoComum/naoCoincedente", grau, sep = "/"))
  return(Nrow_val)
}

###############################################

function_valida <- function(df,seed_sample, seed_apriori, nomearq){
  
  df_sample1 <- df_sample(df, seed_sample)
  
  df_bootstrap <- df[df_sample1, ] #66% da base
  
  regrasApriori <- Apriori_func(df_bootstrap, seed_apriori)
 
  regrasElimRedun <- regrasApriori[ !is.redundant( regrasApriori, measure="confidence" ) ]
  
  regrasPreparadas<- rule_prep(regrasElimRedun)

  #Pegar varios em cima de lift 1.1
  regrasPreparadas <- subset(regrasPreparadas, lift >= 1.1, select= c(LHSandRHS,count,lift, confidence, support))
  #Compara duas bases e retorna uma base valida
  return(regrasPreparadas)
}

run = 60
valor <- data.frame( matrix(nrow = run, ncol = 2))

for (i in 1:1) {
  
  seedValue1 = as.numeric(format(Sys.Date(), format="%d"))
  seedValue2 = as.numeric(format(Sys.Date(), format="%d"))*2
  
  coincedente = function_valida(df, seedValue1, seedValue2, i)
  
  GeneralInfo <- data.frame(runCount=i, ValorSeed_1=seedValue1, ValorSeed_2=seedValue2, Coincedente1=nrow(coincedente), Coincedente2=0, RegraComum=nrow(coincedente), RegraNaoComum=0)
  
  valor[i,1] <- i
  valor[i,2] <- nrow(coincedente)
  print(paste(valor[i,1:2]))
  
  for (k in 2:run) {
    
    seedValue1 = seedValue1+3
    seedValue2 = seedValue2+3
    Regras_valida1 <- coincedente
    Regras_valida2 <- function_valida(df, seedValue1, seedValue2, k)
    
    #Intersect Das k frames de Regras 
    coincedente = merge(coincedente, Regras_valida2, by.x="LHSandRHS", by.y = "LHSandRHS")
    coincedente_df = as_tibble(coincedente, .name_repair = "minimal")
    
    nrow_naocomum <- regrasNaoCoincedente(Regras_valida2,coincedente_df,k)
    
    GeneralInfo[nrow(GeneralInfo)+1,] <- list(k, seedValue1, seedValue2, nrow(Regras_valida1), nrow(Regras_valida2) ,nrow(coincedente), nrow_naocomum)
    
    valor[k,1] <- k
    valor[k,2] <- nrow(coincedente)
    print(paste(valor[k,1:2]))
    
    if(valor[k,2] == valor[k-1,2])
      return(print("Criterio de parada encontrado"))
    
  }
  plot(valor,type = "o")
  
}
write.csv(GeneralInfo, "Apriori_Redundancia_NaoComum/coincedente/GeneralInformacaoTreino.csv")
#NOTE:: Exporta o arquivo e separe o atributo "Desfecho" para outra coluna
write.csv(coincedente, "Apriori_Redundancia_NaoComum/coincedente/coincedente.csv")
#-----------------------
#-------------------------

coincedentes <- read.csv(file="Apriori_Redundancia_NaoComum/coincedente/coincedente.csv", header = TRUE, sep =",")

#newfile1 <- within(coincedente0.1, LHSandRHS<-data.frame(do.call('rbind', strsplit(as.character(LHSandRHS), ' ', fixed=TRUE))))
newfile1 <- data.frame(do.call('rbind', strsplit(as.character(coincedentes$LHSandRHS),' ',fixed=TRUE)))
newfile1$X2 = NULL
newfile2 <- select(coincedentes,-c(1,2))
newfile3 <- coincedentes$X
newfile <-  data.frame(cbind(newfile3, newfile1, newfile2)) 
names(newfile)[1] <- "ID"
names(newfile)[2] <- "LHSandRHS"
names(newfile)[3] <- "Desfechos"

write.csv(newfile, "Apriori_Redundancia_NaoComum/coincedente/newcoincedente0.03&0.05.csv")


desfechos0 <- subset(newfile, Desfechos == "encerra2=0", select=(ID:support.y.1))
desfechos1 <- subset(newfile, Desfechos == "encerra2=1", select=(ID:support.y.1))
desfechos2 <- subset(newfile, Desfechos == "encerra2=2", select=(ID:support.y.1))
desfechos3 <- subset(newfile, Desfechos == "encerra2=3", select=(ID:support.y.1))

write.csv(desfechos0, "Apriori_Redundancia_NaoComum/arquivo/desf0FullRules.csv")
write.csv(desfechos1, "Apriori_Redundancia_NaoComum/arquivo/desf1FullRules.csv")
write.csv(desfechos2, "Apriori_Redundancia_NaoComum/arquivo/desf2FullRules.csv")
write.csv(desfechos3, "Apriori_Redundancia_NaoComum/arquivo/desf3FullRules.csv")

#----------------------------------------------------------------

df <- select(newfile,ID, LHSandRHS)

df4 <- df %>% mutate(ID =  row_number(), flg = 1) %>%
  separate_rows(LHSandRHS, sep = ",") %>%
  spread(LHSandRHS, flg)

df4$a=NULL
a <- select(newfile, Desfechos:support.y)
matrix_coin <-  data.frame(cbind(df4, a)) 

desfechos00 <- subset(matrix_coin, Desfechos == "encerra2=0", select=(ID:support.y))
desfechos11 <- subset(matrix_coin, Desfechos == "encerra2=1", select=(ID:support.y))
desfechos22 <- subset(matrix_coin, Desfechos == "encerra2=2", select=(ID:support.y))
desfechos33 <- subset(matrix_coin, Desfechos == "encerra2=3", select=(ID:support.y))


write.csv(desfechos00, "Apriori_Redundancia_NaoComum/desfecho/desf0Cura_coincedenteMatrix.csv")
write.csv(desfechos11, "Apriori_Redundancia_NaoComum/desfecho/desf1Abandono_coincedenteMatrix.csv")
write.csv(desfechos22, "Apriori_Redundancia_NaoComum/desfecho/desf2Falencia_coincedenteMatrix.csv")
write.csv(desfechos33, "Apriori_Redundancia_NaoComum/desfecho/desf3Obito_coincedenteMatrix.csv")


###################TODO
#Count each item and save into new df
#Plot a graph showing
desfechoCura <- desfechos00 %>% select( -c(Desfechos:ncol(desfechos00)))
desfechoAbandon <- desfechos11 %>% select( -c(Desfechos:ncol(desfechos11)))
desfechoFalencia <- desfechos22 %>% select( -c(Desfechos:ncol(desfechos22)))
desfechoObito <- desfechos33 %>% select( -c(Desfechos:ncol(desfechos33)))

matrizToFreqGraph <- function(matrixDF){
  countAttr <- sapply(matrixDF, FUN=table)
  countedItemsDF <- as.data.frame(t(as.data.frame(do.call(cbind, countAttr))))
  countedItemsDF <- tibble::rownames_to_column(countedItemsDF, "items")
  colnames(countedItemsDF)[2] = "frequencies"
  decreOrdfreqDataSet <- countedItemsDF[order(countedItemsDF$frequencies, decreasing=TRUE),]
  return(barplot(decreOrdfreqDataSet$frequencies[1:10], names.arg = decreOrdfreqDataSet$items[1:10]))
}

matrizToFreqGraph(desfechoCura)
matrizToFreqGraph(desfechoAbandon)
matrizToFreqGraph(desfechoFalencia)
matrizToFreqGraph(desfechoObito)




