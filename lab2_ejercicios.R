#---------------------LABORATORIO 2-------------------------

#PROBLEMA 1:


#PROBLEMA 2:
myEIDs<-c("1","10","100","1000","37690")
(mySymbols<- AnnotationDbi::mget(myEIDs, org.Hs.egGO,ifnotfound=NA ))
names(myGO_all)

#PORBLEMA 3: 
#LAs tres ontologias son GOID, ontology y evidence. 
#La ontología a la que pertenece, que puede ser: 
#proceso biológico (BP: Biological Process), función molecular 
#(MF: Molecular Function) or componente celular (CC: Cellular 
#Components).

#PROBLEMA 4:
#Los nombres de cada elemento de 
#la lista se corresponden con los identificadores GO. Así, 
#para el primer gen.
names(myGO_All[[1]])

#para obtener info de todos los genes:
lapply(myGO_All,names)

#si queremos construir un vector:
unlist(lapply(myGO_All,names))

#para eliminar los repetidos:
unique(unlist(lapply(myGO_All,names)))

#PROBLEMA 5:

########################## SECUENCIAS ##################################3

#Notaciones de los acidos nucleicos:
DNA_ALPHABET

#------------PROBLEMA 1:---------------------------------------------
ds<-paste(sample(DNA_ALPHABET[c(1:4,16)],30,replace=TRUE),collapse = "")

#------------PROBLEMA 2-------------------------------------------
(ds <- DNAString(ds)) #Ppara convertir  wl string en un objeto de clase dnastring

#apartado a)
alphabetFrequency(ds) #vemos cuento se repiten cada uno de la cadena DNA ALPHABET


#apartado b):
complement(ds)
reverse(ds)
reverseComplement(ds)

#apartado c):
subseq(ds, start=c(3),end=c(11))

#apartado d): 
(ds.v <- Views(ds,start=c(3,6,10),end=c(8,15,20)))# podemos utilizar el argumento width en lugar de end

#apartado e):
ds.v[[1]] #accedemso al primer elemento guardado en el xstring 
ds[1] #obtenemos en primer valor de al lista de aminacidos

#vector ds.set con 5 secuencias aleatorias:
set.seed(223)
ds.set <- NULL
for(i in 1:5)
  ds.set <- c(ds.set,
              paste(sample(DNA_ALPHABET[c(1:4,16)],sample(10:15,1),replace=TRUE),collapse=""))
ds.set
ds.set<-DNAStringSet(ds.set)

###################### CONTENIDO GC EN GENOMA ###########################
if(!("BSgenome.Hsapiens.UCSC.hg19" %in% installed.packages()))
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
library("BSgenome.Hsapiens.UCSC.hg19")

# solo nos interesa el cromosoma 8
seqChr8 <- unmasked(Hsapiens$chr8) #secuencia sin enmaskarar 
cpglocs <- read.table("model-based-cpg-islands-hg19.txt",header=TRUE)

# solo nos interesa el cromosoma 8 (columna 1)
# las columnas 2 y 3 son las posiciones iniciales y finales de las islas CpG
cpglocs8 <- cpglocs[which(cpglocs[,1]=="chr8"),2:3]
head(cpglocs8, n=5)

#dimension, ver las islas con las que contamos;
dim(cpglocs8)

#Con la primera funcion hacemos un string con una coleccion de secuencias
#a partir de la secuencia del cromosoma 8 uytilizando las posiciones de 
#inicio y fin. 
seqChr8Islands <- DNAStringSet(seqChr8, start=cpglocs8[,1], end=cpglocs8[,2])

#con la segunda calculamos la frecuencia,la funcion vocuntpattern(x,y) 
#cuenta las veces que aparece x en la lista y. y lo dividimos entre la longitud de
#y para obtener la frecuencia. 
freqIslands <- vcountPattern("CG", seqChr8Islands) / width(seqChr8Islands)


# --------PROBLEMA 3--------------------------------------------------
matriz<-matrix(0,2856,2,byrow=FALSE)
for valor in matriz:
  
