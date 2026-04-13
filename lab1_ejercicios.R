############ PROBLEMAS LAB 1#################3
x<-c(3,6,3,6,7,1,4,7,4,4)

#-------------------- PROBLEMA 1----------------------------------
# a) primer elemento y septimo: 
x[1]
x[7]
# b) Elementos entre posicion 3 y 7: 
x[3:7]

# c)valores de x que ocupan las posiciones 1, 3, y de la 5 a la 8.
x[c(1:3,5:8)]


#-------------------- PROBLEMA 2-----------------------------------
# a)Obten los valores mayores que 4 y menores o iguales que 6.
y<-x[x>4 & x<=6]
y
# b) Menores o iguales que 4 o mayores que 6.
x[x<=4 | x>6]

# c) iguales a 4:
x[x==4]

# d) Valores iguales a 4 pero con funcion which:
x[which(x==4)]


#----------------------PROBLEMA 3 ----------------------------------
# Utiliza la función sort() con el argumento decreasing = TRUE,
# ¿qué indica este argumento?

# Indica que la linea de comando ordena el vector de mayor a menor 

#-------------------PROBLEMA 4 y 5 -----------------------------------
# establecer una semilla aleatoria para la generación de 
# números aleatorios. Esto significa que al utilizar la misma semilla, 
# se obtendrá la misma secuencia 

# El objeto es una matriz. 

#------------------- PROBLEMA 6 ----------------------------------
colnames(x) <- paste0('sample',1:4)
rownames(x) <- paste0('gene',1:10)

#---------------------PROBLEMA 7 ------------------------------------
x[1:10,2] #todos los genes de la muestra 2
x[,2] #tambien 


##########################################################################

#Leer datos en dorma de tabla de un fichero .dat:
mydata <- read.table("file.dat", header = TRUE, sep="\t", row.names = 1)

#--------------------PROBLEMA 9------------------------

#cargamos los datos iris:
data('iris') 
#para conocer el tipo de estructura
str(iris) 
#Guardamos los datos relevantes de iris en un data frame:
# $--> despues de un nombre y seguido de un nombre de columna nos da los datos en esa columna.
myiris <- data.frame(Sepal.Length = iris$Sepal.Length, Sepal.Width = iris$Sepal.Width, Species = iris$Species)

#-------------------PROBLEMA 10 -------------------------
#Acceder a las columnas de longitud y de anchura de sepalo mas especie:
iris[,c(1,2,5)]

#obtener todas salvo la ultima:
iris[1:4]
iris[,-c(5)]

#-------------------PROBLEMAS 11.Añadir columas y filas
#añade una nueva variable:

Stalk.Length<-c(rnorm(30,0,0.1),rnorm(30,1,0.1),rnorm(30,1.3,0.1),rnorm(30,1.8,0.1),rnorm(30,2,0.1))
myiris<-cbind(iris, Stalk.Length)

#añadir nuevo indivduo: 
new<-data.frame(Sepal.Length=10.1, Sepal.Width=0.5, Petal.Length=2.5, Petal.Width=0.9, Species="myspecies")
myiris<-rbind(iris,new)

#-------------------PROBLEMA 12.-----------------------
#obtener un subgrupo del iris que cumpla una condicion:
subset(iris,Sepal.Length>7.5)

############### OPERACIONES ESTADISTICAS ############################

#--------------------PROBLEMA 13-------------------- 
summary(iris)

#--------------------PROBLEMA 14. Funciones sobre la variable longitud del sepalo
#ME HACE LA MEDIA POR FILAS: (de todas las variables)


#para hacerlo con todas:
apply(iris[,1:3], 2, mean) # por columnas (por eso el 2)
apply(iris[,1:4], MARGIN=1, FUN=mean) # por filas (por eso el 1)

apply(iris[,1:4], 2, mean)
apply(iris[,1:4], 2, sd)
apply(iris[,1:4], 2, min)
apply(iris[,1:4], 2, max)
apply(iris[,1:4], 2, median)
apply(iris[,1:4], 2, quantile)


#para hacerlo con una columa: 

mean(iris[,1])
sd(iris[,1])
quantile(iris[,1]) #puedo añadirle los porcentajes
median(iris[,1])
range(iris[,1]) # o con min y max

#--------------------PROBLEMA 15-------------------
#Calculamos correlaciones. (de Pearson por defecto)

cor(iris[1],iris[2])
cor(iris[1:4]) #correlacion de las cuantitativas

#--------------------PROBLEMA 16----------------------
#Calcula la matrizd e varianza-covarianza para als variabels cuantitativas:
cov(iris[1:4])
