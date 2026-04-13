###################### LABORATORIO 7 ##########################

#utilizamos GEO (de bioconductor):
if(!("GEOquery" %in% installed.packages()))
  BiocManager::install("GEOquery")
library("GEOquery")

#Leemos el dataset: 
geo <- GEOquery::getGEO("GSE24460", AnnotGPL = TRUE) #almacenamos en objeto geo
class(geo)

#podemos cargar geo desde el escritorio: 
load("geo.RData")
geo<-geo[[1]] #OBJETO EXPRESSIONSET

#vemos todas las caracteristicas que esconden  los objetos de
#clase expressionset:
n_muestras<-sampleNames(geo)
n_genes<-featureNames(geo)
fenotipo<-pData(geo)
inf_completa<-View(geo)
expresion<-exprs(geo)


#tambien podemos extraer info con arrayexpress:
if(!("ArrayExpress" %in% installed.packages()))
  BiocManager::install("ArrayExpress")
library("ArrayExpress")
array<-ArrayExpress::ArrayExpress(geo)

###################################################################
##################### ANALISIS PCA ##################################
#geo<-geo[[1]]  #nos quedamos con el elemento de clase expressionset

#1. objeto con la matriz de expresion calculada antes: 
expresion<-exprs(geo)
# para que aparezcan las variables en las filas, trasnponemos: 

#2. componentes principales: 
componentes<-prcomp(t(expresion))
head(componentes)

#summary para obtener resumen de PCA:
summary(componentes)

#Con las dos primeras variables ya somos capaces de explicar el 98%
# de la variabilidad, nos podemos quedar con estas dos. 

# 3. REPRESENTAR PARA EVALUAR PATRON EN DATOS:
library(ggplot2)

#como las coordenadas de cada muestra se han guardado en myPca$x, sabemos
#que tenemos que representar myPca$s[:,1] y myPca$s[:,2].

#Crear data.frame con los datos a representar:
resultados<-data.frame(PC1=componentes$x[,1],PC2=componentes$x[,2],
                       Group= c(rep("Control",2), rep("Caso",2))) #con gropu luego podemos dar color 

#Representar: 
ggplot(resultados, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size=3) +
  labs(x = "PC1",
       y = "PC2")

########################################################################
################# ANALISIS DE EXPRESION DIFERENCIAL ####################

#utilizamos paquete limma: Analisis basado en modelos lineales
if(!("limma" %in% installed.packages()))
  BiocManager::install("limma")
library("limma")

if(!("antiProfilesData" %in% installed.packages()))
  BiocManager::install("antiProfilesData")
library("antiProfilesData")

#leemos data base:
data(apColonData)
apColonData

#crear tabla de frecuencias:
#las variables fenotipocas se obtienen con pdata(): !!!!!!!!!!!!!!!!!!!!!
tabla<-table(pData(apColonData)$Status)
tabla

#METODO DE LIMMA: AJUSTAR MODELO LINEAL ANALIZAR EXPRESION DIFERENCIAL.

# 1. Definir la matriz:
design <- model.matrix(~0 + factor(pData(apColonData)$Status))
#class-->matrix array

# Definir la matriz de diseño con subtipo:
subtipo<-table(pData(apColonData)$SubType)
design2 <- model.matrix(~0 + factor(pData(apColonData)$SubType))
colnames(design2) <- levels(factor(pData(apColonData)$SubType))
head(design2,n=3)
#una columna por cada nivel de factor.

# 2. Ajustar modelo en cada gen 
ajuste<-lmFit(apColonData, design)
class(ajuste) #MArrayLM limma

# 3. Estadistico empirical bayes t-test: 
estadistico<-eBayes(ajuste)
class(estadistico)  #misma clase

# 4. Ordenar los genes diferencialmente (mayor a menor signficacion):
orden_genes<-topTable(estadistico, adjust="fdr", sort.by = "B", number=Inf) 

#como resultado tenemos un data.frame corresponde con la cantidad ded epxresion diferencial
table(pData(orden_genes)$Status)

DE <- orden_genes[orden_genes$adj.P.Val<0.0001 & abs(orden_genes$AveExpr)>2,]
dim(DE)
rownames(DE)[1:10] 

###################################################################
############## PROBLMEA CON MULTIPLES ESTADOS ########################

#diseño de matriz: 
design2<-model.matrix(~0 +factor(pData(apColonData)$SubType))

#modificamos los nombres de las columas porque lo sutilizaremos como identificadores de grupo: 
colnames(design2)<-levels(factor(pData(apColonData)$SubType))

# 3. Ajusta los modelos lineales: 
ajuste_previo<-lmFit(apColonData,design2) #lmfit(datos, diseño matriz para datos)

# Definir pares de niveles a comparar. --> CONSTRUYE MATRIZ DE CONTRASTE
contrast.matrix<-makeContrasts(adenoma-normal,colorectal_cancer-normal,tumor-normal,
                               levels=c(colnames(design2)))
contrast.matrix

fit2 <- contrasts.fit(ajuste_previo,contrast.matrix) 

# estadisitico bayes sobre el modelo defintiivo: 
fitE2<-eBayes(fit2)
head(fitE2$coefficients)

#genes diferenciados expresados: 
orden2_genes<-topTable(fitE2, adjust="fdr", sort.by = "B", number=Inf, coef=1)

#selecionare genes diferencialmente expresados: 
tested2 <- orden2_genes[orden2_genes$adj.P.Val<0.0001 & abs(orden2_genes$AveExpr)>2,]
tested2

##### REPRESENTAR RESULTADOS DE ANALISIS DE EXPRESION DIFERENCIAL ##########
plot(tested2$logFC,-log10(tested2$adj.P.Val),pch=20,
     xlim=c(-10,10),ylim=c(0,15),
     xlab="log2foldchange",ylab="-log10p-value")
abline(v=c(-2,2),lwd=2,col="pink") #LINEAS VERTICALES
abline(h=-log10(0.0001),lwd=2,col="purple") #LINEA HORIZONTAL
#enrojo losgenesDE
points(tested2$logFC[is.element(rownames(tested2),rownames(DE2))],-log10(tested2$adj.P.Val[is.element(rownames(tested2),rownames(DE2))]),
       pch=20,col="purple")

#version mas sencilla con volcano plot:
volcanoplot(fitE2,pch=20,
            highlight=5, names=row.names(fitE2))

##################################################
####### ANALISIS CLUSTER CON MICROARRAYS #######

#matriz de expresion.
apColon.DE<-exprs(apColonData)[rownames(DE),]

# 1. Defitnir medida de la distancia o disimilaridad.
# Agerupamos aquellos genes con nivles de expresion mas 
# relacionado. 

matriz <- function(x) as.dist((1-cor(t(x)))/2)
 
# 2. Funcion clustering:
clust.fun <- function(x) hclust(x, method="complete")


# 3. REPRESENTAR MAPA DE COLOR
# mapa de colores
hmcol <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(256)

# tipo de muestra
labels <- pData(apColonData)$SubType
colores<-c("#dfeaf4","#f4dfdf","#f2cb98","#e0d9f4")[as.numeric(as.factor(labels))]


# mapa de calor
hp <- heatmap.2(apColon.DE,                  # datos
                scale="row",                 # estandarizamos por gen (row)
                density.info="none",         # no se muestra la densidad en la leyenda de color
                trace="none",                # no se muestran líneas dentro heat map
                col=hmcol,                   # se define el mapa de colores
                labRow=FALSE,                # no se muestran las etiquetas por filas
                dendrogram="both",           # se muestra el dendograma en ambos ejes
                labCol=labels,               # etiquetas para las muestras
                ColSideColors=colores,       # cada tipo de muestra de un color
                hclust=clust.fun,            # definicion del metodo de clustering
                distfun=matriz,            # definición de la distancia
                margins = c(8,3))

#listado de genes de cada cluster:
clusters <- cutree(as.hclust(hp$rowDendrogram),k=10) # 10 grupos
table(clusters)

# genes del cluster 5
names(clusters)[clusters==5]

