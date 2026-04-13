#########################################################################3
##########################  EJERCICIOS LAB 3  ###############################################

# ------------ANALISIS DE CONTENIDO DE SECUENCIA------------------------------------------

#---------------------Problema 1.----------------------------------
#Leer fichero fasta.
SeqFastadna <- seqinr::read.fasta(file = "den1.fasta")

#lognitud de la secuencia:
longitud<-length(SeqFastadna[[1]])

#---------------------Problema 2.------------------------------------
dengueseq <- SeqFastadna[[1]]
dengueseq[(length(dengueseq)-20+1):length(dengueseq)]

#---------------------Problema 3.-------------------------------------
#frecuencia absoluta:
abs<-table(SeqFastadna)

#frecuencia relativa:
rel<-(abs/longitud)

#contenido en GC:
GC(SeqFastadna[[1]])

#---------------------Problema 4.------------------------------------
complementaria<-comp(SeqFastadna[[1]])
#frec en cadena complementaria.
abs_comp<-table(complementaria)

#---------------------Problema 5.------------------------------------
GC(complementaria)
#en total los pares PARA LA COMPLEMENTARIA:
pares<-seqinr::count(complementaria, wordsize=2) #CC=787 y CG=500

#en la sprimeras 100 posiciones:
pares_100primeros<-seqinr::count(complementaria[0:100], wordsize=2)

#---------------------Problema 6.------------------------------------
#Frecuencia de codones observada:
codones<-seqinr::count(dengueseq, wordsize=3)
count(dengueseq,3)["gac"]

#Frecuencia de codones esperada de GAC:
frec_esperada<-seqinr::rho(dengueseq,wordsize=3)["gac"] #-->1.06 esta proximo a 1

###################################################################33
# ----------------ALINEAMIENTO DE PARES DE SECUENCIAS------------------------------------------
library(seqinr)
library("Biostrings") # para DNAString

#----------PROBLEMA 1. ALINEAMIENTO DE PARES E DOAGRAMA DE PUNTOS--------

#Constuir dos secuencias aleatorias de 100 nucleotidos:
secuencia1 <- sample(c("a","c","g","t"),100,rep=TRUE,prob=c(0.1,0.4,0.4,0.1))
secuencia2 <- sample(c("a","c","g","t"),100,rep=TRUE,prob=c(0.1,0.4,0.4,0.1))
secuencia1_inversa<-reverse(secuencia1)
#DIAGRAMA DE PUNTOS:
par(mfrow=c(1,2))
dotPlot(secuencia1, secuencia2, wsize = 1,wstep = 1) #primera contra segunda
dotPlot(secuencia1, secuencia1, wsize = 1,wstep = 1) # y primera contra primera

#Alineamiento se pares:
(myAlignProt.local <- pairwiseAlignment(secuencia1, secuencia1, substitutionMatrix = "BLOSUM62", 
                                        gapOpening = gapOpen, gapExtension =gapExtend, 
                                        type="local", scoreOnly = FALSE))

#DIBUJAMOS EN UN GRAFICO SOLO: 
par(mfrow=c(1,1)) # volvemos a representar un unico gráfico


#APARTADO b)
#Convertimos secuencia a string:
(seq1.DNAString <- DNAString(paste(secuencia1,collapse="")))

#invertimos:
(seq1.inv <- reverse(seq1.DNAString))

#metemos en una lista:
(seq1.inv <- unlist(strsplit(paste(seq1.inv,collapse=""), split="")))

par(mfrow=c(1,2))
dotPlot(seq1, seq1, wsize = 1,wstep = 1)
dotPlot(seq1, seq1.inv, wsize = 1,wstep = 1)


#--------------PROBLEMA 2. SECUENCIA FICHERO---------------------------
#APARTADO A)
seq_fasta <- seqinr::read.fasta(file = "Alig_Ej2.fasta")
brugiaseq<-seq_fasta[[1]]
loaseq<-seq_fasta[[2]]

#APARTADO B) su alineamiento optimo global:
data(BLOSUM50)                              
brugiaseqstring <- c2s(brugiaseq)           
loaseqstring <- c2s(loaseq)                 
brugiaseqstring <- toupper(brugiaseqstring) 
loaseqstring <- toupper(loaseqstring) 
(myglobalAlign <- pairwiseAlignment(brugiaseqstring, loaseqstring, substitutionMatrix = "BLOSUM50",
                                    gapOpening = -9.5, gapExtension = -0.5, scoreOnly = FALSE))


#APARTADO C)
(alineamiento_brugloa <- pairwiseAlignment(brugiaseqstring, loaseqstring, substitutionMatrix = "BLOSUM62", 
                                           gapOpening = -9.5, gapExtension =-0.5, 
                                           type="global", scoreOnly = TRUE))



#APARTADO D)
randomseqs <- generateSeqsWithMultinomialModel(brugiaseqstring,1000)
randomscores <- double(1000) 
for (i in 1:1000) {
  score <- pairwiseAlignment(loaseqstring, randomseqs[i], substitutionMatrix = "BLOSUM50",
                             gapOpening = -9.5, gapExtension = -0.5, scoreOnly = TRUE)
  randomscores[i] <- score
}
sum(randomscores >= 752.5) #vemos si se obtienen muchos valore ssuperiores al obtenido en el estudio 

#NO HAY NINGUNO POR LO QU EEL P VALOR ASOCIADO ES MUY PEQUEÑO Y EN CONSECUENCIA EL ALINEAMIENTO ES SIFDNCIATIVO