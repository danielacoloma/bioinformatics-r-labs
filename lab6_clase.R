#----------LAB 6 ESTRUCTURA DE PORTEINA CON R-------------------------------
install.packages("bio3d",dependencies=TRUE) #dependencies =TRUE para instalar todas las depenendias 
library(bio3d)

#1. leemos fichero pbd: 
pdb<-read.pdb("4q21")

#2. Atributos de un objeto
attributes(pdb)

#clase de pdb:
class(pdb)

#si veo que pdb me sale la informacion, si me meto en atributes me da la info 
#atom es un data frame que contiene atomo a atomo , tengo el numero de residuo 
plot.bio3d(pdb$atom)
(ca.inds<-atom.select(pdb,'calpha')) #extraigo los carbonos alfa
attributes(ca.inds)


# 3. Seleccionar todos los atomos excepto moleculas de agua:
nowat.ind<-atom.select(pdb,string="water",inverse=TRUE)
nowat.ind

#atom tiene por fila y columna info de ese atomo, podemos 
#acceder a la info de cada uno con los corchertes:

#buscamos los atomos que son de agua:
head(pdb$atom[nowat.ind$atom,])
help(atom.select)

# 4. Combinar proteina y GDP:
atom.select(pdb, string="protein", resid="GDP",operator="OR")
#cuando utilizo varios residuos con varios argumentos utilizo operator qaue peude ser "and" o "or"

# 5. Combine select:
seleccion<-combine.select(string="protein", resid="GDP", operator="OR")

###################################################################
#-----------------TRABAJAR CON MULTIPLES FICHEROS-----------------
# leer los dos pdb
a <- read.pdb("4q21")
b <- read.pdb("4lhy")

## trim.pdb() crea un nuevo objeto PDB basado en una seleccion
# cadena A de 4q21:
a1 <- trim.pdb(a, chain="A")

#Cadenas A,E,F:
b1 <- trim.pdb(b, chain="A")
b2 <- trim.pdb(b, chain="E")
b3 <- trim.pdb(b, chain="F")

#creamos nueva proteina juntando varias estructuras que
#podriamos estudiar:
new <- cat.pdb(a1, b1, b2, b3, rechain=TRUE)
# rechain para asignar nuevos identificadores de cadena
unique(new$atom$chain)

#para guardar:
write.pdb(new, file="4Q21-4LHY.pdb")

# 6. Determinar sitio de union: 
sitio<-binding.site()

# 7. Leer archivo con multiples modelos. 
modelos<-read.pdb("1d1d",multi=TRUE)

# 8. Regiones con menor flexibilidad--> rigidas. 
rigidos<-geostas(modelos)

# 9. Identificar la region más invariante. 
invariante<-core.find(modelos)

###########################################################################3
#------------------VISUALIZACION 3D DE PROTEINAS-------------------
#paquete disponible en cran: 
if(!("NGLVieweR" %in% installed.packages()))
  install.packages("NGLVieweR")
library("NGLVieweR")

#para ver representacion tridimensional utilizamos: 
#------>  con %>% concatenamos instrucciones
NGLVieweR("1BG2") %>%
  addRepresentation("cartoon")


# 2. Representar proteina superponiendo representacioens: 
NGLVieweR("1BG2") %>%
  addRepresentation("cartoon") %>%
  addRepresentation("ball+stick",
                    param = list(sele = "233-248",
                                 colorValue = "yellow",
                                 colorScheme = "element")) %>%
  addRepresentation("label",param = list(sele = "20",
                                         labelType = "format",
                                         labelFormat = "[%(resname)s]%(resno)s",
                                         labelGrouping = "residue",
                                         color = "white", fontFamiliy = "sans-serif",
                                         xOffset = 1, yOffset = 0, zOffset = 0,
                                         fixedSize = TRUE, radiusType = 1, radiusSize = 1.5)) %>%
  setSpin()

# 3. Añadir animacion:
#setSpin() #-->giro constante
#setRock() --> hace va y ven 


# 4. Cambiar el aspecto de algunos residuos.
NGLVieweR("1BG2") %>%
  addRepresentation("cartoon",
                    param= list(colorScheme="residueindex")) #modifica el color de la estructura

#con param podemos repsresntar solamente un conjunto de residuos con el argumento sele e indicando las posiciones:
#hemos cambiado esquema de repserensacion y selecionado los residuos a represnetar
NGLVieweR("1BG2") %>%
  addRepresentation("ball+stick",  
                    param = list(sele = "233-248", #indicamos solo ciertos residuos
                                 colorValue = "pink",
                                 colorScheme = "element"))


#para superponer a la proteina la region de interes




###########################################################3
#-------ESQUEMATIZAR LA INFO DISPONIBLE EN UNIPLOT------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("drawProteins")

# 1. Obtener características de una proteina:
rela<-prot.data<-drawProteins::get_features("Q04206")

# 2. Transofrmamos a data.frame:
rela.data<-drawProteins::feature_to_dataframe(rela)

# 3. Esquematizar la info.

#Vamos a ir añadiendo capas. 
#Lo primero que tenemos que ejecutar es al funcion
#draw_canvas que nos representa el background pero no dibuja nada:
#DEBEMOS EJECUTARLA SIEMPRE PARA QUE FUNCIONE UNA REPRESENTACION:
area<-drawProteins::draw_canvas(rela.data)

# 3.2 draw_chains() añade el rectangulo que es la cadena,
#TAMBIEN MNECESARIO PARA QUE FUNCIONE:
p<-drawProteins::draw_chains(area,rela.data)
p

#3.3 dibuja dominio de la proteina: 
#Esto es un ggplot a este objeto p podriamos modificarle cosas
#al igual que hariamos con cualquier ggplot. 
dom<-drawProteins::draw_domains(p, rela.data)
dom

#Podemos modificarle personalizando el tamaño, fondo, ejes...
library("ggplot2")
p <- p + theme_bw(base_size = 12) + # white background
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank()) +
  theme(axis.ticks = element_blank(),
        axis.text.y = element_blank()) +
  theme(panel.border = element_blank())
p

########## AÑADIR OTRAS CARACTERISTICAS#############3

# 1. Podemos añadir regiones de interes:
a<-drawProteins::draw_regions(p, rela.data)
a

#2. podemos añadir motivos, la señal de localizacion nuclear y el 
#dominio de transactivacion de 9 aminoacidos:
p2<-drawProteins::draw_motif(p,rela.data)
p2

#3. Añadimos los sitios de fosforilacion--> Los pone como puntos los puntos de fosfo
#selecionamos un tamaño no muy grande de punto.
p3<-drawProteins::draw_phospho(p, rela.data, size=8)
p3

#todo esto es añadible podmeos dibujar todo en la misma figura

#4. Elementos comunes conf. tridimiensional
tri<-drawProteins::draw_folding(p,rela.data)
tri

# AÑADIMOS TODO EN UN AMISMA FIGURA:
rela.subtitle <- paste0("circles = phosphorylation sites\n",
                        "RHD = Rel Homology Domain\n
                        source:Uniprot")
p <- drawProteins::draw_regions(p, rela.data)
p <- drawProteins::draw_phospho(p, rela.data, size = 8)
p <- p + labs(title = "Rel A/p65",
              subtitle = rela.subtitle)
p

####################################################################
#-------------------COMPARAR VARIAS PORTEINAS----------------------
prot.data<-drawProteins::get_features("Q04206 Q01201 Q04864 P19838 Q00653")
prot.data<-drawProteins::feature_to_dataframe(prot.data) 
p<-drawProteins::draw_canvas(prot.data)
p

#legend.position--> podemos cambiar la ubicacion de la leyenda
#añaduimos titulos tambien 


prot.subtitle <- paste0("circles = phosphorylation sites\n",
                        "RHD = Rel Homology Domain\n
                        source:Uniprot")
p <- p + labs(title = "Schematic of human NF−kappaB proteins",
              subtitle = prot.subtitle)



#################################################################
#---------------------EJERCICIOS PROPUESTOS.--------------------
#necesitamos instalar drawprotein cosa que no me deja hacer
kv<-drawProteins::get_features("P22001")
kv.data<-drawProteins::feature_to_dataframe(kv)

#primero dibujar el area:
area<-drawProteins::draw_canvas(kv.data)

#tambien rectangulo, cadena para que funcione 
p<-drawProteins::draw_chains(area,kv.data)
p

# esquema de la organizacion:
org<-drawProteins::draw_recept_dom(p,kv.data)
org

# COMPARAR ESPECIES:
especies<-drawProteins::get_features("P22001 P15384 P16390")
especies.data<-drawProteins::feature_to_dataframe(especies)
area1<-drawProteins::draw_canvas(especies.data)
p<-drawProteins::draw_chains(area1,especies.data)
esq_especie<-drawProteins::draw_recept_dom(p,especies.data)
esq_especie

#-------------------VISUALIZACION 3D-----------------------
NGLVieweR("73J1") %>%
  addRepresentation("cartoon")

