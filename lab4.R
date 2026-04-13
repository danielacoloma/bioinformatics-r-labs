library('SNPassoc')
install.package('rms')
data(asthma, package = "SNPassoc")
str(asthma, list.len=9)
data(SNPs, package="SNPassoc")
mySNP<-setupSNP(SNPs,colSNPs=6:40,sep="")
class(mySNP)
association(casco~snp10002,data=mySNP, model="all")

# (leer tabla) el riesgo de enfermear entre los que portan CC es 1.76 veces mayor 
# que TT
#ninguno es estadisticamente significativo porque los intervalos de confianza contienen el 1
# AIC Critewrio de Información çde Arcaique, 
# para ver si es significativo hay que fijarse en el intervalo de confianza y 
# en el pvalor


association(casco~snp10001+sex+blodd.pre, data=mySNP, model = "all")
# riesgo asumiendo como constante el sexo y la presión sanguínea (entre 
# indiviudos del mismo sexo y misma presión sanguínea)

#por ej queremos dominante y codominante
association(casco~snp10001, data=mySNP, model = c("dominant", "codominant"))

#casco es caso control
#las variables confusoras son presión sanguínea y sexo
WGassociation(casco~sex+blood.pre,data=mySNP, model ="all")
# donde pone monomorphic tienen valor constante en todos los individuos

plot(WGassociation(casco~sex+blood.pre,data=mySNP, model="all"))
tableHWE(mySNP,as.factor(mySNP$casco))
plotMissing(mySNP)
#MARCA EN NEGRO LOSMISSIN
association(casco~dominant(snp10001)*factor(sex), data = mySNP)
#lo q quiero evaluar e sla relacion del nip1 con caso control es diferente para hombres y mujeres, poreso factor sex
if(!("GWAStools" %in% installed.packages()))
  BiocManager::install(("GWASTools"))
library("GWAStools")

pVal <- runif(1000)

qqPlot(pVal)

set.seed(25683)
n <- 1000
pvals <- sample(-log10((1:n)/n), n, replace=TRUE)
chromosome <- c(rep(1,100), rep(2,150), rep(3,80), rep(4,90), rep(5,100))
manhattanPlot(pvals, chromosome, signif=1e-5)