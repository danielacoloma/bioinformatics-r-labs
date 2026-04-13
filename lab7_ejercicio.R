###############################################
##########   EJERCICIOS LAB 7 #################
if(!("GEOquery" %in% installed.packages()))
  BiocManager::install("GEOquery")
library("GEOquery")


# EJERCICIO 1.---------------------------------------------------------
datos<-GEOquery::getGEO("GSE76999",AnnotGPL = TRUE)

# EJERCICIO 2. Solo 12 primeros arrays.---------------------------------
geo<-datos[[1]] #OBJETO EXPRESSIONSET
geo<-(geo[,1:12])
pData(geo)[,c("title","tissue:ch1")]   # poner titulo 


#comprobamos mediante: 
pData(geo) #--> se han seleccionado 12

#bariables con info relevante:
# en title y tissue:ch1
typeSample <- factor(pData(geo)[,"tissue:ch1"],
                     levels=c("Bone Marrow","Fetal Liver","Yolk Sac"),
                     labels=c("BM-MOs", "FL-MOs", "YS-Macs"))

# EJERCICIO 3. Estudiar genes diferenciados.----------------------------
matriz_diseño<-model.matrix(~0 +typeSample)

# cambiamos los nombres de los grupos
tmp <- as.factor(paste("G", as.numeric(typeSample), sep=""))  
design2 <- model.matrix(~0 + tmp)
colnames(design2) <- levels(tmp)

#modificamos los nombres de las columas porque los utilizaremos como identificadores de grupo: 
colnames(matriz_diseño)<-levels(factor(typleSample)

# Ajusta los modelos lineales: 
ajuste_previo<-lmFit(geo,design2) #lmfit(datos, diseño matriz para datos)

# Definir pares de niveles a comparar. --> CONSTRUYE MATRIZ DE CONTRASTE
contrast.matrix<-makeContrasts(G1-G2,G1-G3,G2-G3,levels=design2)
contrast.matrix

fit2 <- contrasts.fit(ajuste_previo,contrast.matrix) 

# estadisitico bayes sobre el modelo defintiivo: 
fitE2<-eBayes(fit2)
head(fitE2$coefficients)

#genes diferenciados expresados: 
orden2_genes<-topTable(fitE2, adjust="fdr", sort.by = "B", number=Inf)


# EJERCICIO 4. Seleccionar genes diferencialmente epxresados
testear<-topTable(fitE2, adjust="fdr", sort.by = "B", number=Inf,coef=1)
de <- testear[testear$adj.P.Val<0.0001 & abs(testear$logFC)>2,]
nrow(test)

testear2<-topTable(fitE2, adjust="fdr", sort.by = "B", number=Inf,coef=2)
de2 <- testear2[testear2$adj.P.Val<0.0001 & abs(testear2$logFC)>2,]
nrow(test)

testear3<-topTable(fitE2, adjust="fdr", sort.by = "B", number=Inf,coef=3)
de3 <- testear3[testear3$adj.P.Val<0.0001 & abs(testear3$logFC)>2,]
nrow(test)

# 5. Volcano plot para cada comparacion--------------------------------
par(mfrow=c(1,3)) #rep juntas

plot(testear$logFC,-log10(testear$adj.P.Val),pch=20,
     xlim=c(-10,10),ylim=c(0,15),
     xlab="log2foldchange",ylab="-log10p-value", main="MB frente a FL")
abline(v=c(-2,2),lwd=2,col="navy")
abline(h=-log10(0.0001),lwd=2,col="pink")
#enrojo losgenesDE
points(testear$logFC[is.element(rownames(testear),rownames(de))],-log10(testear$adj.P.Val[is.element(rownames(testear),rownames(de))]),
       pch=20,col="pink")

#--------------------------------------
plot(testear2$logFC,-log10(testear2$adj.P.Val),pch=20,
     xlim=c(-10,10),ylim=c(0,15),
     xlab="log2foldchange",ylab="-log10p-value", main="MB frente a YS")
abline(v=c(-2,2),lwd=2,col="navy")
abline(h=-log10(0.0001),lwd=2,col="pink")
#enrojo losgenesDE
points(testear2$logFC[is.element(rownames(testear2),rownames(de2))],-log10(testear2$adj.P.Val[is.element(rownames(testear2),rownames(de2))]),
       pch=20,col="pink")

#-----------------------------------------------
plot(testear3$logFC,-log10(testear3$adj.P.Val),pch=20,
     xlim=c(-10,10),ylim=c(0,15),
     xlab="log2foldchange",ylab="-log10p-value", main="FL frente a YS")
abline(v=c(-2,2),lwd=2,col="navy")
abline(h=-log10(0.0001),lwd=2,col="pink")
#enrojo losgenesDE
points(testear3$logFC[is.element(rownames(testear3),rownames(de3))],-log10(testear3$adj.P.Val[is.element(rownames(testear3),rownames(de3))]),
       pch=20,col="pink")
 


#6. ANALISIS DE COMPONENTES PRINCIPALES ----------------------

myPca <- prcomp(t(exprs(geo)))
summary(myPca)
plot(componentes$x[,1:2], col=as.numeric(typeSample), pch=19)
legend("topright", inset=.05, levels(typeSample), pch=19, col=1:3, horiz=TRUE)

# 7. ANALISIS DE CLUSTER:

