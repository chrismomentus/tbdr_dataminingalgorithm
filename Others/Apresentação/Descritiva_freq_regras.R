# Descritiva das regras não coincidentes do algoritmo do Chris

# 06 08 2024


setwd("E:/0__ASSUNTO/Alunos/X_13_PIBIC_TCC_Chris_2018/Analises ago 2024")

library (data.table)  # carregar arquivos csv, txt
library (foreign)     # carrega arquivos dbf e outros formatos;exporta dataframe para vários formatos de pacotes estatísticos
library (epiDisplay)  # tabelas para categóricas oranizadas e gráficos
library (gtsummary)   # tabelas resumo organizadas com tipos diferentes de variáveis

naocoinc <- read.csv("E:/0__ASSUNTO/Alunos/X_13_PIBIC_TCC_Chris_2018/Analises ago 2024/NaoCoincedenteElimRedundancia.csv", sep=',')

table(naocoinc$Frequencia, useNA = 'ifany' )

naocoinc$perc <- naocoinc$freq/51*100

