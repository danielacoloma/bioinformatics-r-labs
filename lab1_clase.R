#---------------------------PRACTICA 1----------------------------
#------------------------Introduccion a R-------------------------

#Instalar librerias o paquetes de R:
#Son una funcionalidad que se ha creado implementada para resolver un problema concreto. 

#bioconductor--> reservorio de paquetes para problemas especificos de la bioinformatica
#cran--> revervorio pars apquetes en general

#PAQEUTES:
#1) INSTALAR:
#se instalan de distinta forma depdnediendo de cual sera el reservorio
#los de cran-->install.packages('nombre') tambien desde el boton de packages/install
#los de bioconductor-->
install.packages('FMM')

#2) CARGAR:
library('FMM')

getOption('defaultPackages')

#para ver cuantos hemso instalado: 
print(.packages()) #vemos cuantos estan cargados en ese momento

#La ayuda de R es muy buena, importante usarla. 
#para que un paqeute de R pueda poners een cran y biocon todas las ayudas estan documentadas. 
help("mean") 
?mean

#NOTA: los argumentos se pueden pasar bien  por posicion o por nombre, cuando escribimos el primero 
#por nombre, todos deben ir por nombre. 

#La forma mas facil de crear un vector, es la funcion c, la funcion de CONCATENAR:
x<-c(3,6,3,6,7,1,4,7,4,4)  #OPERADOR ASIGNACION <-, PERO TAMBIEN VALE =
class(x) #clase de argumentos que tiene x ( numero, caracter..)

#---------------------ACCEDER A ELEMENTO DE UN VECTOR----------------------

x[1] #valor en la posicio 1, el primer elemento del vector x
#--VECTOR CON ELEMENTO DE X--
x[1:3]
x[c(1:3,5)] #del 1 al 3 con el de la posicion 5

#(NOTA---> R opera elemento a elemento)

#----------CONDICIONES PARA SELECIONAR ELEMENTO DE VECTOR----------------------
x>=4  #devuelve verdadero en aquuellas posiciones donde se cumple 
x[x>=4] #devuelve el valor almacenado en las posiciones donde se cumple la condicion 
y<-x[x>4 & x<=6] #metemos dentro de y el resultado para guardarlo 
x==4 #estrictamente igual 
!(x==4) #complementario de x==4

#DETERMINAR LA POSICION DONDE SE ENCUENTRA UN VALOR: 
which(x==4) #nos da la sposiciones en las que se encuentra 4
x[which(x==4)] 


#-------------------------------------------------------------------------------
#------------------------OPERACIONES INTERESANTES------------------------------
#ORDENAR VECTOR:
sort(x) #menor a mayor
sort(x, decreasing=TRUE) #de mayor a menor, decreciendo verdadero
sort(x, decreasing=TRUE, index.return=TRUE) #nos devuelve el orden de donde se encontraban los valores 
#en la lista x original.


#----------------------DISTRIBUCIÓN DE PROBABILIDAD---------------------------------
n.data<-rnorm(n=100, mean=1, sd=0.1)
head(n.data) # visualizamos los primeros elementos
hist(n.data)
plot(density(n.data))


#-------------------------------------MATRICES-------------------------------------
#Simulamos 10 genes con 4 muestras de cada uno con una distribucion de poisson:
set.seed(123)
x <- matrix(rpois(10*4,lambda=5),nrow=10)  #40 distribuciones generadas por poisson en 10 filas
x[1,3] #fila y columna 
x[1:2,3:4]
x[,1]

## Funcion byrow=TRUE:
matrix(1:9,nrow=3)
matrix(1:9,nrow=3,byrow=TRUE) #observaciones por filas en vez de por columas

#hay un argumento SI QUEREMOS QUE SIGA MANTENIENDO LA CLASE de ese objeto
#si queremos que mantega matriz y no vector:
x[,1,drop=FALSE]

#CONCATENACIÓN DE STR:
paste0('sampe', 1:4) 

#poner nombre a las filas y col:
colnames(x) <- c("sample1","sample2","sample3","sample4")
rownames(x)<-paste0('gene',1:10)

#buscamos por nombre en vez de posicion: 
x["gene6","sample2"]

#----------------------------LISTAS--------------------------------
geneset1 <- c(1,4)
geneset2 <- c(2,5,7,8)
geneset3 <- c(1,5,6)

#acceder a varios elementos a la vez es complicado:
#LISTA DE LISTAS:
(genesets <- vector("list",3)) #creamos 3 lista vacia
genesets[[1]] <- geneset1
genesets[[2]] <- geneset2
genesets[[3]] <- geneset3
genesets
names(genesets)<-paste0('geneset',1:3) #damos nombre a cada genset
genesets

#-----------------------LAPPY Y SAPPLY------------------------------
sapply(genesets, FUN=length)  #combina los reusltados
lapply(genesets,FUN=length) #devuelve lista
unlist(lapply(genesets,FUN=length)) #sacar de lista el lapply

#-------------------------DATA_FRAME----------------------------
#Abrir un fichero en R no signfica que podamos utilizarlo

#ABRIMOS EL DATA FRAME IRIS (disponible en paquetes basicos):
data('iris')
class(iris)
head(iris)

#Si tenemos nuestros propios datos podemos almacenarlo en un 
#fichero .dat o .txt, de texto plano.

iris$Species #devuelve columa que se llama especies
iris[1:2,'Species']

#-----------------------------------------------------------
#LEER DATOS:
load(file="golub.dat")
#mydata <- read.table("NOMBRE FICHERO", header = TRUE, sep="\t", row.names = 1)
mydata <- read.table("golub.dat", header = TRUE, sep="\t", row.names = 1)
#---------------------------------------------------

#-------------------LEYENDO Y ESCRIBIENDO DATOS------------------------
load(file=".RData")
#Para guardar un objeto de datos:
save(DATA, file=".RData")


# Leer fichero de excel: 
install.packages("xlsx", dependencies=TRUE)
library(xlsx)
mydata <- read.xls("mydata.xls")


#--------------FUNCION APPLY---------------------
#ME HACE LA MEDIA POR FILAS: (de todas las variables)
apply(iris[,1:4], MARGIN=1, FUN=mean) # por filas
apply(iris[,1:4], 2, mean) # por columnas

#-------------CONTRASTE DE HIPOTESIS------------------------
#Mediante un CH utilizo una muetsra para repsponder una pregunta 
#que me etoy haciendo sobre la poblacion. 

#1) definimso el contraste
#2)ndemostras si nuestros datos apoyaban una posicion o la contraria
#Las funciones que se utilizan en R son bastante homogeneas

#----CHI CUADRADO----------
#PAsamos como argumento al labra de contingencia
tbl <- matrix(c(14, 33, 7, 3), ncol = 2)
tbl
#damos nombre a las columnas: 
colnames(tbl) <- c("Familiar", "Deportivo")
rownames(tbl) <- c("Hombre", "Mujer")
tbl

test <- chisq.test(as.table(tbl))
test

#-------------WILCOXON-------------
dat2 <- data.frame(
  Beginning = c(16, 5, 15, 2, 14, 15, 4, 7, 15, 6, 7, 14),
  End = c(19, 18,20, 17, 10,17, 16, 19, 20, 9, 11, 18)
)
dat2

test <- wilcox.test(dat2$End, dat2$Beginning, paired = TRUE, alternative = "greater")
test
#el p valor: 
test$p.value

#----------------VISULAICACION DE DATOS---------------
plot(x=iris[,"Sepal.Length"], y=iris[,"Petal.Width"], xlab="Petal length", ylab="Sepal length",
     col="black", main="Variation of sepal length with petal length")
model<-lm(iris[,"Petal.Width"]~iris[,"Sepal.Length"]) # ajuste del modelo lineal
abline(model,col="tomato",lwd=2) # representación de la recta en el plot
#___________________________________________________________________
# 100 observaciones N(1,0.1), 100 N(2,0.1) y 50 N(3,0.1)
set.seed(1234)
genex <- c(rnorm(100, 1, 0.1), rnorm(100, 2, 0.1), rnorm(50, 3,0.1))


#----------------------DIAGRAMA DE CAJA----------------
boxplot(Sepal.Length~Species, data=iris, ylab="sepal length",
        xlab="Species", main="Sepal length for different species")

#-------------------------HISTOGRAMA-------------------
set.seed(2345)
x <- rnorm(1000, 3, 0.02) # 1000 observaciones de una N(3,0.02)
hist(x, probability=TRUE)
lines(density(x), col="red")

