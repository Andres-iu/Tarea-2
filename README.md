# Tarea-2
tarea archivo NEXUS
> setwd("C:/Users/Alumno24-B/Desktop")
> install.packages("seqinr")
> library(seqinr)
> install.packages("ape")
>  library("ape")
>  data<-read.nexus.data("sequence_nex.nexus")
> data
> num_sequencies<-length(data)
> num_sequencies
> longitud_secuencia <- unique(sapply(data,length))
> longitud_secuencia
> library(msa)
> dna_sequences <- as.DNAbin(data)
> dna_sequences
>  library(Biostrings)
> sequences <- sapply(data, function(x) paste(x, collapse = ""))
>  dna_set <- DNAStringSet(sequences)
> dna_set
> aliniamiento <- msa(dna_set, method = "ClustalW")
> aliniamiento
> aliniamiento_convertido <- as(aliniamiento, "DNAMultipleAlignment")
> dna_bin <- as.DNAbin(aliniamiento_convertido)
> dist_matrix <- dist.dna(dna_bin, model = "raw")
> tree <- nj(dist_matrix)
> plot(tree, main = "Arbol filogenético")
![54873a87-f634-4448-bbad-236e2df54023](https://github.com/user-attachments/assets/1b34dc0d-dd1d-41b1-8ae5-bb1b0919c7c7)
>  motivo <- "atgcg"
> ids_con_motivo <- names(sequences)[grepl(motivo, sequences)]
> ids_con_motivo
> motivo <- "agcat"
> ids_con_motivo <- names(sequences)[grepl(motivo, sequences)]
> ids_con_motivo
> install.packages("pegas")
>  diversidad_nucleotidica<-nuc.div(dna_bin)
> library(pegas)
> diversidad_nucleotidica<-nuc.div(dna_bin)
> diversidad_nucleotidica<- nuc.div(dna_bin)
> cat("Diversidad nucleotidica (π):", diversidad_nucleotidica,"/n")
> alignment_matrix <- as.matrix(dna_bin)
> variable_sites <- which(apply(alignment_matrix, 2, function(col) length(unique(col[col != "-"])) > 1))
> variable_sites
>  haplotipos <- haplotype(dna_bin)
> haplotipos
>  red_haplotipica <- haploNet(haplotipos)
> plot(red_haplotipica, main = "Red haplotípica", show.mutations = TRUE)
![da1b6d87-6825-4433-ab9c-ed110fa5da58](https://github.com/user-attachments/assets/f90e0e81-02dd-494f-8399-adbc057c3b47)

> frecuencias_nucleotidos <- alphabetFrequency(dna_set, as.prob = TRUE)
> frecuencias_nucleotidos
> frecuencia_total <- colSums(frecuencias_nucleotidos)
> frecuencia_total
> tajima_result <- tajima.test(dna_bin)
> tajima_result
