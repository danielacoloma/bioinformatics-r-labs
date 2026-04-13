#-----------------LABORATORIO 2----------------

#cargar paquetes:
library(AnnotationDbi) 
library(org.Hs.eg.db)

#generar listado de identificadores:
myEIDs<-c("1","10","100","1000","37690")

#varios tipos de anotaciones disponibles en Entrez:
ls("package:org.Hs.eg.db") 

(mySymbols<- AnnotationDbi::mget(myEIDs, org.Hs.egSYMBOL, ifnotfound=NA)) #NA para que, aunque haya missing data, sigan saliendo los demás valores por pantalla

unlist(mySymbols) #nos devuelve el contenido en forma de vector

#volvemos a asignar a mysymbols el resultado en vector: 
mySymbols<-unlist(mySymbols)

#funciones is. devuelve verdadero cuando encuentra un NA.
is.na(mySymbols)

#el complementario de is.na---> !is.na.

mySymbols[!is.na(mySymbols)] #elimina los que son NA

#Asignamos en y ese vector de caracteres si en NA: 
y<-mySymbols[!is.na(mySymbols)]
(myEIDs <- unlist(mget(y, org.Hs.egSYMBOL2EG, ifnotfound=NA)))

AnnotationDbi::mget(y, org.Hs.egSYMBOL2EG, ifnotfound=NA)

#si lo quiero en forma de vector: 
unlist(AnnotationDbi::mget(y, org.Hs.egSYMBOL2EG, ifnotfound=NA)) 

#------------------------------------------------------------------
#-------------------BASE DE DATOS GO------------------------------------------

#instalado y cargado go y paquetes que correspondan
y<-mySymbols

# comprobamos si lo tenemos instalado antes de hacer nada
if(!("GO.db" %in% installed.packages()))
  BiocManager::install("GO.db") 
library(GO.db) #cargamos el GO

ls("package:org.Hs.eg.db") 

#creamos vector: 
myEIDs<-c("1","10","100","1000","37690") #identificadores de los genes

#Para obtener en que categorias de go estan dando estos genes utilizo: 
(mySymbols<- AnnotationDbi::mget(myEIDs, org.Hs.egGO,ifnotfound=NA ))  #importante el ifnotfound
#mySymbol<-AnnotationDbi::mget(nombredelvector , loquequeremosbuscarenlalistade ls )


myGO_all<-unlist(mySymbols)
myGO_all[!is.na(mySymbols)] #Solo los que no tengan NA
names(myGO_all) #me da los nombres de los elementos de la lista
names(myGO_all[[1]])


#obtengo los 132 genes que estan anotados con el gen 1:
AnnotationDbi::mget(names(myGO_all[[1]]), org.Hs.egGO2EG)  

#-------------BUSQUEDAS EN PUBMED-----------------------
#CARGAMOS LIBRERIAS
#instalamos rismed, todas ellas estan EN CRAN, NO EN BIOCONDUCTOR:

if (!("RISmed"%in% installed.packages()))
  install.packages("RISmed")
#hay otra manera mas simple de instalarlo
library("RISmed")  #IMPORTANTE CARGAR CON LIBRARY
if(!("rentrez" %in% installed.packages()))
  install.packages("rentrez")
library("rentrez")
if(!("dplyr" %in% installed.packages()))
  install.packages("dplyr")
library("dplyr")

#de la misma manera instalamos las demas de la pagina 3 que nos piden 

search_topic<-'(BRCA1[TI] AND Breast Cancer [TI] AND Humans) AND"2023"[EDAT]'
#EN STR (STRING) VALEN TANTO COMILLAS SIMPLES COMO DOBLES

#BUSQUEDA: utilizamos paquete rismed, hay dos funciones que utilizamos par ahacerlo: 

search_query<-RISmed::EUtilsSummary(search_topic, type='esearch', db='pubmed') #funcion que busca el resumen del topico que buscamos 
records<-RISmed::EUtilsGet(search_query, type='efetch', db="pubmed")

#OBTENEMOS EL NUMERO DE RESULTADOS CON ESOS CRITERIOS DE BUSQUEDA--> Records:  56 

class(records)

#----------------ANALISIS ESTADISTICO DEL CONJUNTO EXTRAIDO-------------------
#variables: 
pmid<-PMID(records) #argumento el objeto medline y devuelve identificador de pubmed
titulo<-ArticleTitle(records) #titulo de los articulos 
mes<-MonthEntrez(records) #mes de publicacion 
dia<-DayEntrez(records) #dia en el que se han publicado
tipo<-PublicationType(records)
doi<-DOI(records)
idioma<-Language(records)
pais<-Country(records)

#vamos a juntar las variables en el dataframe: 
datos_busqueda<-data.frame(pmid=PMID(records),
                           mes=MonthEntrez(records),
                           dia=DayEntrez(records),
                           titulo=ArticleTitle(records),
                           doi=DOI(records),
                           idioma=language(records),
                           pais=Country(records),
                           tipo=PublicationType(records))

head(datos_busqueda,n=2) #para buscar dos publicaciones

# Añadir la información del número de citas de cada uno de los artículos para obtener el NUMERO DE CITAS:
pmid_summary<-entrez_summary(db='pubmed', id=pmid)

#fuerza a que sea un numero lo que me devuelve:
as.numeric(extract_from_esummary(pmid_summary, 'pmcrefcount')) 
### [1] 41
#MEDIANTE PMREFCOUNT--> OBTENEMOS EL NUEMRO DE CITAS 

(num_citas <- as.numeric(extract_from_esummary(pmid_summary, "pmcrefcount")))

is.na(num_citas) #verdadero donde hay NA
num_citas[is.na(num_citas)]<-0 #sustituye los NA por un 0, porque tiene 0 referencias 

datos_busqueda$num_citas<-num_citas #ella lo llama pubmed_data

head(datos_busqueda,n=2)


#--------------------PAQUETE DT: FACILITA SINTAXIS BASICA--------------------
#PAQUETE DT PARA VER TABLA: 
# para presentarlo en forma de tabla necesitamos el paquete DT
if(!("DT" %in% installed.packages()))
  install.packages("DT")
library("DT")

datos_busqueda %>% mutate(title=substr(titulo,1,20)) %>% # y palabras clave a los 20 primeros
dplyr::distinct(pmid,.keep_all = TRUE) %>% DT::datatable()    #ME DA ERROR RARO

#------------------------PAQUERE GGPLOT2----------------------------
#PARA GRAFICAS MUY VISUALES
#tabla:
datos_busqueda %>% mutate(title=substr(titulo,1,20)) %>% # y palabras clave a los 20 primeros
  dplyr::distinct(pmid,.keep_all = TRUE) %>% DT::datatable()    #ME DA ERROR RARO


#grafica:
#no concuerda si condigo con mis nombres, mirarlo!!!!!!
library(dplyr)
library(ggplot2)

datos_busqueda %>% distinct(pmid,.keep_all = TRUE) %>% group_by(MES=factor(mes)) %>%
  summarise(n=n()) %>%
  ggplot(aes(x=MES,y=n,fill=MES)) +  #esta linea es obligatoria
  geom_bar(stat = "identity") +  #voy a representar barras
  geom_text(aes(label=n),nudge_y = 1) +  #para que me añada el número en cada barra
  theme_bw() + theme(legend.position = "none") +
  ggtitle("Número de publicaciones por mes")

#numero de publicacion: 
datos_busqueda %>% distinct(pmid,.keep_all = TRUE) %>% group_by(MES=factor(mes)) %>%
  summarise(n=n())

#ABSTRACT DEL ARTICULO 1:
getAbstracts(datos_busqueda$pmid[1]) #primer articulo

#FUNCION CLEAN ABSTRACT: nos permite construir esta distribucion de frecuencias: 
abstract<-getAbstracts(datos_busqueda$pmid)
#obtengo todos los datos de esos articulos (los abstract) y los he metido dentro de el abstrac

length(abstract)

#RECUENTO PALABRAS EN ABSTRACT:
head(cleanAbstracts(abstract))
#con esto hace recuento de palabras 

#hacemos objeto con distrubucion de frecuencias:
cleanAbs<-cleanAbstracts(abstract) 
class(cleanAbs) #es un data frame que se lo pasamos a la funcion siguiente: 


# NUMERO DE CITACIONES: 
datos_busqueda %>% distinct(pmid,.keep_all = TRUE) %>% group_by(tipo) %>% 
  summarise(Citation=sum(num_citas)) %>% na.omit() %>%
  ggplot(aes(x="",y=Citation,fill=tipo)) +
  geom_col() + coord_polar(theta = "y") +
  theme_bw()+
  geom_text(aes(label = Citation),
            position = position_stack(vjust = 0.5)) +
  theme(legend.position="bottom",legend.text = element_text(size=8)) +
  ggtitle("Número de citaciones por tipo de publicación")

#REPRESENCIACION NUBE DE PALABRAS CLAVE:
#(instalamos paquete wordcloud2)
wordcloud2(cleanAbs)
