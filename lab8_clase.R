#############################################################33
######### LABORATORIO ANALISIS ESTADISTICO DE ADN  ############3

#Muy importante la noralizacion, 
#significa hacer todas las muestras comparables.
#objetivo--> elimar sesgos sistematicos para quedanros solo con la variabilidad biológica


#PORBLEMA DE COMPARACIOENS MULTIPLES: 
#No puede ser nunca mayor que 1 el ajuste

BiocManager::install("edgeR")
library('edgeR')


BiocManager::install("pasilla")
library('pasilla')

#funcion que me devuelve str con esta ruta: 
datafile<- system.file("extdata/pasilla_gene_counts.tsv", packages="pasilla")


# 1.1 leer fichero---------------------------------
#el nombre de la fila es el nombre del gen con row.names=1:
rawCountTable<- read.table(datafile,header=TRUE, row.names= 1) #indica donde tengo almacenado el nombre de la fila
dim(rawCountTable)
head(rawCountTable)


# 1.2 Leer los datos fenotipicos de la condicion:-----------------------------
#como unico aegumento la ruta completa del archivo:
SampleAnno <- read.csv(file.path(system.file("extdata", package = "pasilla"),
                                 "pasilla_sample_annotation.csv"))
SampleAnno




##############3############################################
###### EXPLORACION DE DATOS Y CONTROL DE CALIDAD ##########

# 1. transofrmacion de los datos 
#EXPLORACION DE LOS DATOS DE CONTROL DE CALIDAD:
pseudoCount<- log2(rawCountTable+1)

#para construir gg plot, necesito tener toda la info en el mismo dataframe:
# data frame con los datos a representar:
df <- data.frame(rcounts = rawCountTable[ ,1], prcounts = pseudoCount[ ,1])

library(ggplot2)
library(gridExtra)

#quiero que calcule la densidad, frecuencia relativa de rcounts, no le doy
#nombre de variabel sino un calculo que debe hacer con los datos.
#si escribo desity con dos puntos espero columna que s ellame densiti

#la siguiente capa es geom histrogram

#tambien puedo modificar el color de las barras con fill y la anchura
#de las clases que forzamos a 2000 mediante binwidth. 

#el resto de capas tienen que ver con el formato opotimizado para poryector. 

p1 <- ggplot(data=df, aes(x = rcounts, y = ..density..)) +
  geom_histogram(fill = "#525252", binwidth = 2000) + 
  theme_bw() +
  ggtitle("Count distribution") +
  xlab("counts")

grid.arrange(p1, p2, ncol=2)

p2 <- ggplot(data=df, aes(x = prcounts, y = ..density..)) +
  geom_histogram(colour = "white", fill = "#525252", binwidth = 0.6) + 
  theme_bw() +
  ggtitle("Log2 Count distribution") +
  xlab(expression(log[2](counts + 1)))

grid.arrange(p1, p2, ncol = 2)

###########################################################
########### DISTRIBUCION DE LAS MEUSTRAS ##################

#Para representar todas las muestras a la vez distribuye 
#mediante reshape ::melt(), pawuete disponibel en cran que podemos instalar
#directamente con install::packages

install.packages("reshape")
library("reshape")


df<-data.frame(samp=colnames(rawCountTable), t(pseudoCount))
dim(df)
#la dimension de este dataframe son 7 filas y 14600 columnas. 
#aplicando reshape tengo 3 columnas y todos los demas valores en 
#la misma fila.
df<-reshape::melt(df)
dim(df)

# añadir una columna condition con el tipo de muestra sin el número
df$condition <- gsub("[0-9]","",df$samp)
head(df,n=15)

p <- ggplot(data = df, aes(x = samp, y = value, fill = condition)) + 
  geom_boxplot() + theme_bw() + ggtitle("count distributions") + 
  xlab("sample") + ylab("counts") + 
  theme(axis.text.x = element_text(angle = 90))
p

#########################################
#importante--> a la vista de un boxplot podemos decir si neciestamos normalizar los datos
#las cajas normalizadas son iguales. En global todas las distibuciones eran iguales para
#que un valor sea comparable con los valores en otras muestas. 

#GRAFICO DE DESNIDAD:
#Tambien es un ggplot con el mismo data-frame y la capa geom es geom_density
#alpha sirve para darle transparecia. 

#CAPA facet_wrap--> consigo colocar variso graficos en la MISMA ESCALA, de forma
#que peudo comparar uno con el otro.
#Me represnta un grafico para los tratado y otro para los no tratadoss pero la grafica la conserva,m pro lo que podemos ver los dos a la vez sueprpuestos.

#La función expression: sirve para aportar notacion matematica en lso titulos

###############################

#MEDIR ACUERDO:
#para medir el acuerdo se utilzian los fraficos de blan alman. Comparamos
#la diferencia entre las mediciones hechas (M) con (A) que coge las medidas y hace la media de las dos

#Me interesa saber en que rango me estoy moviendo, dos aparatos peuden ser itnercambiables
#en rangos bajos pero no en altos, Es importante saber el rango en el que comparamos. 
#Se espera observar que en la mayoria esten de acuerdo 

df <- data.frame("A" = (pseudoCount[ ,3] + pseudoCount[ ,4])/2,
                 "M" = pseudoCount[ ,3] - pseudoCount[ ,4])

p <- ggplot(data = df, aes(x = A, y = M)) + geom_point(alpha = 0.2) +
  theme_bw() + 
  geom_smooth(color = "red3",size=2)+
  geom_hline(yintercept=0, color = "blue3",size=2)
p

#Curva roja--> curva suavizada que ajusto empiricamnete a los datos. 
#DEBEN PARECERSE LO MAS POSIBLE LA ROJA A LA AZUL

#Si la nuve de putnso esta depslazada nos dice que una de las nube de puntos sistematicamente nos esta 
#dando uno de los valores muy diferentes a la otra. UNA DE LAS MUESTRAS SOBREEXPRESADA EN ESTE GEN CONCRETO. 


## FILTRADO DE DATOS:
#Nos quedamos con conjutno de genes que digan informacion relevante, 
#no nos quedemos quedar con lso que no den informacion relevante para ninguna muestra.

keep<- rowSums(rawCountTable)>0  #genes cuya suma es mayor que 0, los 0 se suelen quitar siempre porque no nos dicen nada
#(podemos encontrar umbrales diferentes depende del estudio)

dim(rawCountTable)

filtCount <-pseudoCount[keep,]

# representar otra vez: 


#NORMALIZACION:
dge2<.calcNormFactors(dge2, method="TTM")
dge2

#3,
pseudo.TMM<-log2(cpm(dge2)).....

#ANALISIS EXPRESION DIFERENCIAL
condition<-as.factor(gsub("[0-9]","",colnames(rawCountTable)))
dge<-DGEList(rawCountTable, grup=condition, remove.zeros=TRUE)

dge<-calcNormFactors(dge,)...

#Una vez estimada la varianza aplicamos el test:

#MATRIZ DE DISEÑO:
#CON model.matrix

design.matrix<-model.matrix(~ dge$samples$group +
                              dge$samples$replicate)
design.matrix

#podemos eliminar la constante B0 definiendola como 0


####################

fit<-glmfit(dge, design.matrix)

dgeLRTtest<-glmf

#CLUSTERING
#Clustering jerarquico 
