
#--------------------------PRACTICA 3-------------------------------
##############FICHEROS FASTA-->TODOS TIENEN EL MISMO FORMATO###################

# comprobamos antes si ya los tenemos instalado
if(!("seqinr" %in% installed.packages()))
  BiocManager::install("seqinr")
library("seqinr")

#cambiamos directorio en files:
setwd("C:/Users/Propietario/OneDrive - UVa/INGENIERÍA BIOMÉDICA/3º IB/2 CUATRI/Bioinformática/Laboratorio BINF/Trabajo veterana (ALEJANDRA)/RStudio/lab3")

#Comprobamos las bases de datos disponibles:
seqinr::choosebank()

#Y seleccionamos la primera:
seqinr::choosebank("genbank", timeout = 20)

#Seleccionamos al secuencia del gen humano BRCA1:
BRCA1<-seqinr::query("BRCA1", query = "SP=Homo sapiens AND K=BRCA1")

#load("BRCA1.RData")
BRCA1<-BRCA1
#obtenemos: 
#BRCA1--> data
#myseq--Z values


#accedemos a secuencias de este objeto:
myseq <- seqinr::getSequence(BRCA1$req[[1]])

#Obtener nombre d ela funcion: 
nombre<-seqinr::getName(BRCA1)

#Para cambiar minuscula y mayuscula: 
toupper(myseq)

#--------------PROBLEMA 1:---------------------------------------
#FUNCION write.fasta(): para crear un fichero .fasta con esta informacion. 
fichero<-seqinr::write.fasta(sequences=toupper(myseq), names="sequence1",file.out="MyBRCA.fasta")
fichero2<-seqinr::write.fasta(sequences=toupper(myseq), names=c("sequence1","seq2"),file.out="MyBRCA2.fasta")

#--------------PROBLEMA 2:------------------------------------
#Para cambiar minuscula y mayuscula: 
toupper(myseq)

# -------------PROBLEMA 3:------------------------------------
fichero2<-seqinr::write.fasta(sequences=,toupper(myseq), names=c("sequence1","seq2"),file.out="MyBRCA2.fasta")


#--------------PROBLEMA 4:--------------------------------------------
#cargamos el fichero fasta que hemOS creado con dos secuencias: 
#creamos dichero:
myFastaFile<-read.fasta('MyBRCA2.fasta')  #solo escribimos nomrbe del fichero proque le tengo en mi ruta de trabajo actual



##################### ANALISIS DE CONTENIDO DE SEQ #########################################

choosebank("genbank",timeout=20)

## 1. CREAR LISTAS
#Creamos dos listas con las secuencias completas
# secuencia de la M.tuberculi
actino <- query(listname = "actino", query="SP=Mycobacterium tuberculosis AND K=rpoB")
# secuencia de la E.coli
proteo <- query(listname = "proteo", query="SP=Escherichia coli AND K=rpoB")


## 2. ELEGIR SECUENCIAS DE CADA ESPECIE
#Elegimos dos secuencias de cada especie. Vamos a averiguar la posicion de estas secuencias:
(myActino.pos<-which(getName(actino$req)=="JX303316"))
(myProteo.pos<-which(getName(proteo$req)=="AM946981.RPOB"))

#Obtenemos las secuencias para cada especie:
myActino <- getSequence(actino$req[[myActino.pos]])
myProteo <- getSequence(proteo$req[[myProteo.pos]])

#Comprobamos que son de una longitud comparable:
length(myActino)
length(myProteo)


## 3.  ESCRIBIR FICHERO FASTA CON SEQ DE INTERES
secuencia_interes<-list(myActino,myProteo)
nombres_secuencias<-c("JX303316","AM946981.RPOB")
actino_proteo<-seqinr::write.fasta(sequences=,secuencia_interes, names=nombres_secuencias,file.out="GC_Bacterias.fasta")

#---------PROBLEMA 1: LEE FICHERO Y CREA VECTORES--------
myFastaFile_GC<- seqinr::read.fasta('GC_Bacterias.fasta')
actino<-myFastaFile_GC[[1]]
proteo<-myFastaFile_GC[[2]]
#---------PROBLEMA 2: NUMERO DE APARECIONES DE CADA BASE
table(actino)
table(proteo)
#---------PROBLEMA 3: CONTENIDO DE G+C--------------------
GC_actino<-seqinr::GC(actino)
GC_proteo<-seqinr::GC(proteo)


#NO EJECUTAR---> APUNTES CLASE PARA LOS EJS:########################################
#Accedemos a solo uno de los valores de la lsita:
mynames <- c("JX303316","AM946981.RPOB")
seqinr::write.fasta(sequences = myseqs, names=mynames, file.out = "GC_Bacterias.fasta")
myFastaFile<-seqinr::read.fasta(file="GC_Bacterias.fasta")
myactino<-myFastaFile[[1]]
myproteo<-myFastaFile[[2]]
table(myactino)
myseqs<-list(myactino, myproteo);
table(myactino)/length(myactino)

#me da la frecuencia relatina del contenido de GC en esta secuencia:
seqinr::GC(myactino)  
seqinr::GC(myproteo)

#Calculo de percentil: 
#el 0.975 lo sacamos del intervalode confianza del 95%, 1-alfa
qnorm(0.975, mean=0, sd=1)  #-->obtenemos el 1.96
##################################################################3


#---------PROBLEMA 4: INTERVALO DE CONFIANZA---------------
myseqs<-list(actino, proteo);

#INTERVALO DE CONFIANZA PARA MYACTINO: 
#Porcentaje de C+G en la secuencia:
phat<-GC(actino) #contenido del actino
phat

#lo utilizamos en el calculo del intervalo:
se<-sqrt(phat*(1-phat)/length(myActino))
phat-qnorm(0.975)*se #sería el limite inferior
phat+qnorm(0.975)*se #sería el limite superior
phat+c(-1,1)*qnorm(0.975)*se #me da un vector con el limite superior e inferior

#INTERVALO DE CONFIANZA PARA MYPROTEO:
phat2<-GC(proteo) #contenido de GC en actino
se2<-sqrt(phat2*(1-phat2)/length(myProteo))
phat2-qnorm(0.975)*se2 #sreia el limite inferior
phat2+qnorm(0.975)*se2 #sería el limite superiorç
phat2+c(-1,1)*qnorm(0.975)*se2 #me da un vector con el limite superior e inferior


#CONSTRUIMOS UN DATA FRAME con los valores a representar: 
data <- data.frame(GC=c(GC(myActino), GC(myProteo)),
                   SD=c(sqrt((GC(myActino)*(1-GC(myActino)))/length(myActino)),
                        sqrt((GC(myProteo)*(1-GC(myProteo)))/length(myProteo))),
                   name =c("Actinobacteria","Proteobacteria"))

#importante, calcular el intervalo de confianza con funcion sapply: 
(IC<-sapply(1:2,function(i) data[i,1]+c(-1,1)*qnorm(0.975,mean=0,sd=1)*data[i,2]))

##################### 5. GRAFICO DE BARRAS ####################################
library(ggplot2)
#como son barras de error necesito un minimo y un maximo: 
data <- data.frame(GC=c(GC(myActino), GC(myProteo)),
                   SD=c(sqrt((GC(myActino)*(1-GC(myActino)))/length(myActino)),
                        sqrt((GC(myProteo)*(1-GC(myProteo)))/length(myProteo))),
                   lim.inf=IC[1,],lim.sup=IC[2,],
                   name =c("Actinobacteria","Proteobacteria"))

#añadimos una capa mas a la rperesentacion: 
#ggplot(data)+ geom_bar(aes(x=name, y=GC))
ggplot(data) +
  geom_bar( aes(x=name, y=GC), stat="identity", fill="pink", alpha=0.5) +
  geom_errorbar( aes(x=name, ymin=lim.inf, ymax=lim.sup), width=0.4, colour="orange", alpha=0.9, size=1.3) +
  coord_flip()+
  ggtitle("Contenido de GC en bacteria actino y proteo")


#coord_flip() gira la grafica

##---------------PROBLEMA 6. Frecuencia de palabras---------------
#-----------------------------------------------------------------
#FRECUENCIA ABSOLUTA:
#Calcula la frecuencia absoluta de codones: 
abs<-(seqinr::count(proteo,wordsize=3)) #ponemos 3 porque contamos codones

#FRECUENCIA RELATIVA:
#Si quiero la frecuencia relativa: 
rel<-(seqinr::count(proteo,wordsize=3) /length(proteo))

#------------- PROBLEMA 7. SOBRE O INFRA REPRESENTADA---------

#hay funcion para esto: 
seqinr::rho(myseq,wordsize=2)

#--------------------------------------------------------------------
#-------------------ALINEACIÓN DE SECUENCIAS------------------------

#cargamos Biostrings
library(Biostrings)

#Necesitamos dos secuencias. Definimos las secuencias: 
sequence1 <- "GAATTCGGCTA"
sequence2 <- "GATTACCTA"   #no tienen porque sere igual de largas

#Definimos sintema de puntuacion (creamos unamatriz de sustitucion):
(myScoringMat <- pwalign::nucleotideSubstitutionMatrix(match = 1, mismatch =-1, baseOnly = TRUE))

#Penalización por huecos. Para indicar al algoritmo si queremos que favorezca
#la aparicion de gap en el alineamiento o no. 
#Determinamos: 
gapOpen<-2
gapExtend<-1


#Utilizamos la función para llevar a cabo un alineamiento global 
#Especificamos si quiere que sea global o local en type
(myAlignment <- pairwiseAlignment(sequence1, sequence2, substitutionMatrix = myScoringMat,
                                  gapOpening = gapOpen, gapExtension = gapExtend,
                                  type="global", scoreOnly = FALSE)) #TRUE  solo devuelve puntuacion de alineamiento

#ALINEAMIENTO DE PROTEINAS -->secuencias de aminoácidos: 
data(package="Biostrings") #--> conjuntos de datos dentro de ese paquete que se pueden usar para ilustrar ejemplos
sequence1 <- "PAWHEAE"
sequence2 <- "HEAGAWGHE"

#En este caso vamos a utilizar unamatriz de sustituacion del paquete biostring:
#el cambio que tenemos que hacer en al funcion es en substititionMatriz definir el nombre del tipo de blosum que vamos a utilizar. 

#GLOBAL: 
(myAlignment <- pairwiseAlignment(sequence1, sequence2, substitutionMatrix = "BLOSUM62",
                                  gapOpening = gapOpen, gapExtension = gapExtend,
                                  type="global", scoreOnly = FALSE))
#LOCAL:
(myAlignment <- pairwiseAlignment(sequence1, sequence2, substitutionMatrix = "BLOSUM62",
                                  gapOpening = gapOpen, gapExtension = gapExtend,
                                  type="local", scoreOnly = FALSE))


#------------------------------------------------------------------
#---SIGNIFCACION ESTADISTICA DE UN ALINEAMIENTO GLOBAL-------------

#-Contraste de tipo bootstrap que nos ayudena evlauar como de lejos esta el alinemiwntro de lo que obtendrimaso con secuenci aleatoria. 
#-Simularemos muchas secuencias aleatorias y comprobar si los alinemaientos de la que considero de referencia se aleja mucho o no. 
#-Alineamos ambas secuencias que peuden tener cierta similitud y veremos si es esta alineacion mejor que la que esperariamos por azar. 

#CREAMOS UNA FUNCION QUE GENERA UN NUMERO X DE VECES LA SECUENCIA:
generateSeqsWithMultinomialModel<-function(inputsequence,X){
  require("seqinr")#Necesitamos maquete seqinr
  
  #Cambia la secuencia en un vector de letras:
  inputsequencevector<-s2c(inputsequence)
  
  #Encuentra la frecuencia de las letras:
  mylength<-length(inputsequencevector)
  mytable <-table(inputsequencevector)
  
  #Encuentra el nombre de las letras:
  letters<-rownames(mytable)
  numletters<-length(letters)
  
  #Crea un vector para alamcenar las probabilidades
  probabilities<-numeric()
  for(i in 1:numletters) {
    letter<-letters[i]
    count<-mytable[[i]]
    probabilities[i]<-count/mylength
  }
  
  #MakeXrandomsequencesusingthemultinomialmodelwithprobabilities"probabilities"
  seqs<-numeric(X)
  for(j in 1:X) {
    seq<-sample(letters,mylength, rep=TRUE,prob=probabilities) #Sample con remplazo
    seq<-c2s(seq)           #argumento prob con la probabilidas calculada anres en el vector
    seqs[j]<-seq
  }
  
  #Devuelve vector:
  return(seqs)
}

#PARA EJECUTAR LA FUNCION: 
generateSeqsWithMultinomialModel(sequence1, 3) #EL NUMERO DETERMINA EL Nº DE SECUENCIAS UE TE DA


#---------- PROBLEMA 1. Generar 1000 secuencias aleatorias-----------------------------------------

#Creamos lista con las secuencias generadas pro al funcion: 

set.seed(569)
random_seq<-generateSeqsWithMultinomialModel('PAWHEAE',1000) #EL NUMERO DETERMINA EL Nº DE SECUENCIAS UE TE DA

#----------PROBLEMA 2. aLINEAMIENTO DE SEQ CON LAS RANDOM-------------
seq4 <- "HEAGAWGHEE"
randomscores <- double(1000) # Create a numeric vector with 1000 elements
for (i in 1:1000) {
  score <- pairwiseAlignment(seq4, random_seq[i], substitutionMatrix = "BLOSUM50",
                             gapOpening = -2, gapExtension = -8, scoreOnly = TRUE)
  randomscores[i] <- score
}


#-------PROBLEMA 3. Alineamiento de las random con las originales----------
pairwiseAlignment(seq4, "PAWHEAE", substitutionMatrix = "BLOSUM50", 
                  gapOpening = -2,gapExtension = -8, scoreOnly = FALSE)

pairwiseAlignment(seq4, sequence2, substitutionMatrix = "BLOSUM50", 
                  gapOpening = -2,gapExtension = -8, scoreOnly = FALSE)
#-------PROBLEMA 4. ESTUDIAR PROBABILIDAD DE OBTENER PUNTUACION MAYOR QUE LA OBTENIDA------------
hist(randomscores, col="pink") 

#Para ver cuantas secuencias tienen un alineamiento mayor o igual que -5:
sum(randomscores>=-5)

#LA probabilidad de obtener un valor superior a el -5 es:
prob_masque5<-(sum(randomscores>=-5)/length(randomscores))

