# integrating-functional-genomics-into-phylogenomics

Integrating Functional Genomics into Phylogenomics

Here we will outline the processs to identify where a UCE loci is found in the genome, next to chatagorize it (intronic, exotic, intergenic) . Once the genetic class has been identified we then go through how to merge loci by gene and conduct species tree analyses in ASTRAL-III. Please cite:

*Van Dam M. H., Henderson J. B., Esposito L., Trautwein M. Incorporating functional genomics into UCE  phylogenomics improves gene and species tree reconstructions.  (2018).*

Corresponding author: Matthew H. Van Dam, email: ([matthewhvandam@gmail.com](mailto:matthewhvandam@gmail.com))

[TOC]



### 1. **Identify where a UCE loci is found in the genome**

First we need to extract the base genome's probes from the *.fasta* probe file generated by PHYLUCE. We will use the *get_individual_probes_from_list.sh* script below. You will need to run this in the terminal. 

```shell
########################## get_individual_probes_from_list.sh ############################

## open terminal tab change directory to where your bait files are, type pwd to see full copy and change the below path in script to match yours
 
 nano get_individual_probes_from_list.sh

## paste in below, modify your path for your bait file & the taxon e.g. "athros1"


critter=$1
[ -z $1 ] && critter=athros1

awk -v critter=$critter 'BEGIN{re="^>.*" critter "[0-9]*"}
     $0 ~ re {in_probe=1; print; next}
     /^>/ {in_probe=0}
     in_probe {print}' /home/hymv2baitsetfilescomplete/hymenoptera-v2-PRINCIPAL-bait-set.fasta.txt


######### exit nano Ctrl + o, "return", Ctrl + x

chmod 775 get_individual_probes_from_list.sh
##########################################################
###

```



You will need to modify what "critter"  you want to by its name in the fasta header in the probe file e.g. "*critter=athros1*".

```sh
>uce-125_p10 |design:hymenoptera-v2,designer:faircloth,probes-locus:uce-125,probes-probe:10,probes-source:athros1,probes-global-chromo:KB467657.1,probes-global-start:99052,probes-global-end:99172,probes-local-start:40,probes-local-end:160
```

 In addtion change the last line of the "*get_individual_probes_from_list.sh*" to the path of your probe set, typically I would use the full path. To run the script  and save output simply run by (this assumes it is in your current working directory).

```shell
./get_individual_probes_from_list.sh > athros1_individual_probes.fasta
```

Now you should have a new file called "*apimel4_individual_probes.fasta*" . If you want to list the files by date modified in your directory type 

```sh
ls -lhtr
```

This will be helpfull when you have alot of files generated and you want to find the most recent one in the shell. 



### 1a. **baltq.sh the probes** to find where they are hiding in the "base genome"

Get blatq.sh and install first https://github.com/calacademy-research/BLATq

Then run as follows blatq.sh "genome.fasta" "individual-probes.fasta" "output file.m8".

```shell
blatq.sh GCA_000344095.2_Aros_2.0_genomic.fna athros1_individual_probes.fasta athros_to_athros-uce-probes-matches.m8
```



### 2. Now we are going to **insert introns into the general feature format (GFF) file** of the "base genome"

You can typically find these annotation files from genbank. Here you will want to locate the **RefSeq** assembly for your taxon, and dowload the "**genomic.gff.gz**" file and unzip it into your directory. Once you have it, the introns can then be inserted into the file. This in not nessessary to do for matching the UCEs to genes in the GFF, but will come in handy when you want to get fancy with partitioning or categorizing the UCEs by feature (intronic, exonic, intergenic).

Below is **R** the code we used in the paper (link) and we have an updated version of the code here (https://github.com/matthewhvandam/integrating-functional-genomics-into-phylogenomics/blob/master/add_introns_to_gff.py). It has been tested on linux. 



```R
# interminal type R

####################################### tetrapod UCEs CHICKEN

#set your working direcoty and read in some R libraries 

setwd("~/data/synteny/bird_genomes/annotations")
library(gdata)
library(ape)
library(ips)
library(dplyr)
library(parallel)


##########################################################
################################################   
################ read in .gff and insert introns by gene
################################################   
##########################################################


gff = read.gff('~/annotations/GCF_000002315.4_Gallus_gallus-5.0_genomic.gff', na.strings = c(".", "?"))  ## make sure it is refseq version of gff

gff_mod = gff

geneID =list()
for (i in 1:nrow(gff_mod)){
geneID[[i]] = unlist(strsplit(as.character(gff_mod[i,9]), "[;]", perl=TRUE ))[1]
}

gff_mod = gff_mod[,-9]


gff_mods = cbind(gff_mod, unlist(geneID))
 
generange = cbind(as.numeric(as.character(gff_mods$start)) , as.numeric(as.character(gff_mods$end))) 
generange_max = apply(generange, 1, max)
generange_min = apply(generange, 1, min)
 
gff_mods =  cbind(gff_mods, generange_min, generange_max)

gff_mods$newcolumn = "n"

curgeneid = "n"

class(curgeneid)
gff_mods = data.frame(gff_mods, stringsAsFactors = FALSE)
gff_mods[,9] = as.character(gff_mods[,9])

##### this can take a few hours
for (i in 1:nrow(gff_mods)){ 

	if(gff_mods[i,3 ]=="gene"){
		curgeneid = gff_mods[i,9 ]
	}
	gff_mods[i, 12]= curgeneid 
	print(row.names(gff_mods[i,]))	
}
##### 


##### this can take a few hours

 make_geneid = function(dataframes) {
 dataframes =  gff_mods
 for (i in 1:nrow(dataframes )){ 
	if(dataframes[i,3 ]=="gene"){
		curgeneid = dataframes[i,9 ]
	}
	dataframes[i, 12]= curgeneid 	
	print(row.names(dataframes[i,]))
}}


mclapply(gff_mods, make_geneid, mc.cores=getOption("mc.cores", 93))

write.csv(gff_mods, file="gff_mods_complete.csv")

##### 


exons = gff_mods[gff_mods[,3]=="exon",]

by_gene_intron_exon = data.frame(stringsAsFactors =F)
by_gene_intron_exonL =  data.frame(stringsAsFactors =F)
by_gene_intron_exonLsigs = data.frame(stringsAsFactors =F)
intron = data.frame(stringsAsFactors =F)

###have to deal with trna and rna exons, then read back in 

gff_mods = read.csv("test_gff_mods.csv", header=TRUE, stringsAsFactors = FALSE)
row.names(gff_mods)= gff_mods$X
gff_mods = gff_mods[,-1]

##########################################################
################################################   
########### if much more than 20,000 genes in .gff you may want to break this up into parallel runs, if it gets bogged down 
################################################   
################################################## insert introns into genes 
		
gene_exons = gff_mods[gff_mods[,3]=="exon", ]

for (i in 1:1440027){
	        #1440104 number of rows to go through 
		
## find a gene's start & stop range and extract the exons it contains
	if(gff_mods[i,3 ]=="gene"){
	gene_exons_start = as.numeric(as.character(gff_mods[i, 10]))
	gene_exons_end = as.numeric(as.character(gff_mods[i, 11])) 
	gene_exons_startID = substr(gff_mods[i, 12], 1, 13)
	gene_exone_startID = substr(gff_mods[i, 12], 1, 13)
          
   gene_exons_alls = gene_exons[gene_exons[, 10]==gene_exons_start , ]
    if(nrow(gene_exons_alls) < 1  ) next # bail if no exon in gene for "pseudogenes" 
    	
  gene_exons_alls = gene_exons_alls[substr(gene_exons_alls[, 12],1,13)==gene_exons_startID, ]
  
  
               if(length(row.names(gene_exons[gene_exons[, 11]==gene_exons_end , ])) < 1 ){
               	gene_exons_end = as.numeric(as.character(gff_mods[i+1, 11]))
               
               }
  
  gene_exons_alle = gene_exons[gene_exons[, 11]==gene_exons_end , ] 
   gene_exons_alle = gene_exons_alle[substr(gene_exons_alle[, 12],1,13)==gene_exone_startID, ]
			
		if(nrow(gene_exons_alle) > 1  ){ ### deal with  duplicates & alternative splicings
			gene_exons_alle = gene_exons_alle[nrow(gene_exons_alle),]
			gene_exons_alles = 	gene_exons[as.character(as.numeric(rownames(gene_exons_alls):as.numeric(rownames(gene_exons_alle)))), ]
			# remove any NA rows generated 
  	gene_exons_alles = gene_exons_alles[complete.cases(gene_exons_alles[ , 10]), ]
  			#get all unique pairs, removes duplicate pairs
  	gene_exons_alles = gene_exons_alles %>% distinct(gene_exons_alles, generange_min, generange_max, .keep_all = TRUE)

		}else {
	gene_exons_alles = gene_exons[as.character(as.numeric(rownames(gene_exons_alls):as.numeric(rownames(gene_exons_alle)))), ] 
}
 
    gene_exons_alles = gene_exons_alles[order(gene_exons_alles[,10]), ]
    colnames(gene_exons_alles)[9:12] = c("geneID", "start1", "end1", "geneID1")
    ### now deal with alternative overlapping splicings
			if(sum(duplicated(gene_exons_alles[,10]), na.rm=TRUE) | sum(duplicated(gene_exons_alles[,11]), na.rm=TRUE) >= 1  ){
				gene_exons_alles = as.data.table(gene_exons_alles)
				if(sum(duplicated(gene_exons_alles[,10]), na.rm=TRUE) >= 1){

				gene_exons_alles = gene_exons_alles[ , .SD[end1 == max(end1)], by = start1]
## for what the the .SD means see https://stackoverflow.com/questions/8508482/what-does-sd-stand-for-in-data-table-in-r
				
					}
				if(sum(duplicated(gene_exons_alles[,11]), na.rm=TRUE) >= 1){
				gene_exons_alles = 	gene_exons_alles[ , .SD[start1 == min(start1)], by = end1]
					}
					setcolorder(gene_exons_alles, c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "geneID", "start1", "end1", "geneID1" ))
					gene_exons_alles = as.data.frame(gene_exons_alles, stringsAsFactors = FALSE)
				}				
            
intron = gene_exons_alles
intron = intron[-nrow(intron), ]
if(nrow(intron) >=1 ){ 
intron[,1] = as.character(intron[,1])
}else {
	}

intron = data.frame(lapply(intron, as.character), stringsAsFactors=FALSE)
intron[, 4] = as.numeric(intron[, 4])
intron[, 5] = as.numeric(intron[, 5])
intron[, 10] = as.numeric(intron[, 10])
intron[, 11] = as.numeric(intron[, 11])


## now working only with rows that contain exons of a gene, if # of rows longer than 1, insert intron if exons are non-consecutive
 if(nrow(gene_exons_alles)-1 >=1){
 	
 	#gene_exons_alles[with(ranges, startgene <= startUCE & endgene >= endUCE), ]
 	
       for (i in 1:(nrow(gene_exons_alles)-1))  { ##  for number of exons -1, as after end of last exon is intergenic
    			
    			print(gene_exons_alles[i,])
				intron_start = as.numeric(as.character(gene_exons_alles[i, 11]))+1
				intron_end = as.numeric(as.character(gene_exons_alles[i+1, 10]))-1
				#intron[i,1]=gene_exons_alles$V1 
				intron[i,3]="intron"
				intron[i,4]=intron_start
				intron[i,10]=intron_start
				intron[i,5]=intron_end
				intron[i,11]=intron_end
	      	}
          by_gene_intron_exon =  rbind(data.frame(gene_exons_alles), intron)
  
     }
      
              else {
	     	by_gene_intron_exon_sigs = rbind(gene_exons_alles)  ## no rb 
	     	print(by_gene_intron_exon)
	     	print("IIIIINNNNTTTTTTRRRROOONNNN!!IS!!NNNNAAAAA")
	     	by_gene_intron_exonLsigs = bind_rows(by_gene_intron_exonLsigs, by_gene_intron_exon_sigs)
	     #	## notes from testing trouble shooting ##  #
                  #by_gene_intron_exonLsigs= cbind(by_gene_intron_exonLsigs, as.numeric(as.character(by_gene_intron_exonLsigs[,4])))
	     	#by_gene_intron_exonLsigs= cbind(by_gene_intron_exonLsigs, as.numeric(as.character(by_gene_intron_exonLsigs[,5])))
	      	}
	   by_gene_intron_exonL = bind_rows(by_gene_intron_exonL, by_gene_intron_exon)
 			# by_gene_intron_exonL[[i]] =  by_gene_intron_exon
 			
      } }
      


by_gene_intron_exonLDFs = rbind(by_gene_intron_exonL, by_gene_intron_exonLsigs)  



write.csv(by_gene_intron_exonLDFs, file="by_gene_intron_exonL_tetrapods.csv")



```





### 3. **Match the results from blastq with genomic GFF data**

Here we will use **R** to match the two datasets together, use the csv result from above. The .csv should look something like the example below.

```shell
head by_gene_intron_exonL_tetrapods.csv

"","seqid","source","type","start","end","score","strand","phase","geneID","start1","end1","geneID1"
"1","NC_006088.4","Gnomon","exon",117,317,NA,"+",NA,"ID=id1",117,317,"ID=gene0"
"2","NC_006088.4","Gnomon","exon",1677,1769,NA,"+",NA,"ID=id2",1677,1769,"ID=gene0"
"3","NC_006088.4","Gnomon","exon",1977,2027,NA,"+",NA,"ID=id3",1977,2027,"ID=gene0"
"4","NC_006088.4","Gnomon","exon",2030,2157,NA,"+",NA,"ID=id4",2030,2157,"ID=gene0"
"5","NC_006088.4","Gnomon","exon",2352,2458,NA,"+",NA,"ID=id5",2352,2458,"ID=gene0"
"6","NC_006088.4","Gnomon","exon",2548,2861,NA,"+",NA,"ID=id6",2548,2861,"ID=gene0"
"7","NC_006088.4","Gnomon","intron",318,1676,NA,"+",NA,"ID=id1",318,1676,"ID=gene0"
"8","NC_006088.4","Gnomon","intron",1770,1976,NA,"+",NA,"ID=id2",1770,1976,"ID=gene0"
"9","NC_006088.4","Gnomon","intron",2028,2029,NA,"+",NA,"ID=id3",2028,2029,"ID=gene0"

```



```R
##########################################################
################################################   
################ now match UCEs to GFF exons-introns by gene-ID
################################################   
##########################################################
#set your working direcoty and read in some R libraries 

setwd("~/data/synteny/bird_genomes/annotations")
library(gdata)
library(ape)
library(ips)
library(dplyr)
library(parallel)

## get a filtered m8 file of UCE blastq hist to the genome
gff_mods = read.csv("test_gff_mods.csv", header =T)
gff_mods = gff_mods[,-1]

### find replace " with nothing and same for spaces in file
fh = read.table('Tetrapods_5K_individual_probes_chicken.m8')

m8table_filt = fh[(fh[,3]>99.0 & fh[,4]>119.0), ]
write.table(m8table_filt, file="filt2_noNT.m8", sep="\t")
#nrow(m8table_filt) #should be same as number of probes for the base taxon, assuming this is on the base taxon

m8table_filt = read.table("~/annotations/filt2_noNT.m8", sep="\t")

## get range by scaffold
UCE_range = m8table_filt[,c(2,1, 9:10)]
gene_ranges = gff_mods[gff_mods[,3 ]=="gene", c(1,9, 10:11)]
colnames(UCE_range) = c("scaffold", "feature", "start", "end")
colnames(gene_ranges) = c("scaffold", "feature","start", "end")

#check make sure all start cols less than end cols
gene_ranges[gene_ranges$"start" > gene_ranges$"end", c("scaffold", "feature","start", "end")] <- gene_ranges[gene_ranges$"start" > gene_ranges$"end",  c("scaffold", "feature","end", "start")] 
#gene_ranges

UCE_range = UCE_range[order(UCE_range[,1], UCE_range[,3]), ]
gene_ranges = gene_ranges[order(gene_ranges[,1], gene_ranges[,3]), ]

##match UCE interval with Gene interval

ranges = merge(UCE_range, gene_ranges, by="scaffold", suffixes=c("UCE", "gene"))
overlap = ranges[with(ranges, startgene <= startUCE & endgene >= endUCE), ]

nrow(overlap)
#3296
write.csv(overlap, file="gene_UCE_overlap_birds.csv")

gene_ranges = gene_ranges[order(gene_ranges[,1], gene_ranges[,3]), ]

order_by_geneID = overlap[order(overlap[,5], overlap[,2]), ]

sapply(order_by_geneID, function(x) length(unique(x)))
 
unique_genes = unique(order_by_geneID$featuregene)

get_unique_UCE_order_by_geneID = data.frame(lapply(order_by_geneID, as.character), stringsAsFactors=FALSE)
get_unique_UCE_order_by_geneID[] = lapply(get_unique_UCE_order_by_geneID, as.character)
get_unique_UCE_order_by_geneID = apply(get_unique_UCE_order_by_geneID, 2 , function(y) gsub("_p10", "", y))
get_unique_UCE_order_by_geneID = apply(get_unique_UCE_order_by_geneID, 2 , function(y) gsub("_p9", "", y))
get_unique_UCE_order_by_geneID = apply(get_unique_UCE_order_by_geneID, 2 , function(y) gsub("_p8", "", y))
get_unique_UCE_order_by_geneID = apply(get_unique_UCE_order_by_geneID, 2 , function(y) gsub("_p7", "", y))
get_unique_UCE_order_by_geneID = apply(get_unique_UCE_order_by_geneID, 2 , function(y) gsub("_p6", "", y))
get_unique_UCE_order_by_geneID = apply(get_unique_UCE_order_by_geneID, 2 , function(y) gsub("_p5", "", y))
get_unique_UCE_order_by_geneID = apply(get_unique_UCE_order_by_geneID, 2 , function(y) gsub("_p[0-9]", "", y))
get_unique_UCE_order_by_geneID = apply(get_unique_UCE_order_by_geneID, 2 , function(y) gsub("_p[0-9]", "", y))

get_unique_UCE_order_by_geneID = data.frame(get_unique_UCE_order_by_geneID, stringsAsFactors=FALSE)
get_unique_UCE_order_by_geneID = get_unique_UCE_order_by_geneID[!duplicated(get_unique_UCE_order_by_geneID$featureUCE), ]
                              
write.csv(get_unique_UCE_order_by_geneID, file="ONLY_unique_UCE_order_by_geneID.csv")

get_unique_UCE_order_by_geneID = read.csv("ONLY_unique_UCE_order_by_geneID.csv")
get_unique_UCE_order_by_geneID = get_unique_UCE_order_by_geneID[,-1]

get_unique_UCE_order_by_geneIDList = cbind(as.character(get_unique_UCE_order_by_geneID[,2]), as.character(get_unique_UCE_order_by_geneID[,5]))
get_unique_UCE_order_by_geneIDList = data.frame(get_unique_UCE_order_by_geneIDList)
#unique(get_unique_UCE_order_by_geneIDList[,2])
#2958 


genes_to_concat = list()
unique_genes_list = unique(get_unique_UCE_order_by_geneIDList[,2])
for(i in 1:length(unique(get_unique_UCE_order_by_geneIDList[,2]))){
	genes_to_concat[[i]] = get_unique_UCE_order_by_geneIDList[get_unique_UCE_order_by_geneIDList[,2]==unique_genes_list[i], ]
}


genes_to_concat_keepers = list()
for(i in 1:length(genes_to_concat)){
if(nrow(genes_to_concat[[i]])>1){
	genes_to_concat_keepers[[i]] = genes_to_concat[[i]]
	} else {
	   	  print("do NA")
}}

genes_to_concat_keepers[sapply(genes_to_concat_keepers, is.null)] = NULL
length(genes_to_concat_keepers)

gene_to_be_concatted = genes_to_concat_keepers
gene_to_be_concatted = lapply(gene_to_be_concatted, function(x) { x[[1]] = as.character(x[[1]]) })
#
############ now write each set into its own dir and concat

gene_to_be_concatted_file_set = lapply(gene_to_be_concatted, function(x) { x = paste(x, ".phylip", sep="") })

dir.create("nexus_for_concat")

## don't close R down just yet, we will use these objecs down below

```



### 3a. **Writing the files you have identified to concatenate by gene**

At this point you should have identified the files (UCEs) that are within a gene and now actually need to concatenate those files and put them into directories, so things will not get too crazy all in one directory, as you may have thousands of loci and will run thousands or analyses it will be very difficult to easily find where things are. So, I like to put each UCE or concatenated set of UCEs (those found in a gene like those found above) into its own directory. It also makes moving those directories around easier.

```R

processFilelength <- function(f) {
    #df <- read.phy(f)
    df = read.nex(f)
}
        
                                       
# !!!!!
 # !!!!! read in nexus indidividula alignemtns made by PHYLUCE these should be in one #directory !!!!!
# !!!!!
                                       
                                       
files = dir(".", recursive=TRUE, full.names=TRUE, pattern="\\.nexus$")
#if not too many nexus alignments
result = sapply(files, processFilelength) 
#if thousands of nexus alignments
results =  mclapply(files, read.nex, mc.cores=getOption("mc.cores", 16))                                  

#make names for rewritng nexus to phylip
namess = files
namess = gsub(".nexus",".phylip",namess)
namess = gsub("./","",namess)

#now write .nexus files to .phylip into a directory
for (i in 1:length(result)) {
    namess[[i]] = paste("/synteny/nexus_for_concat/", namess[[i]], sep="")
    write.phy(result[[i]], file=namess[[i]], interleave=F)
}

                                       

matches=list()
sumlist = list()
phylip_file_list = list.files(path = ".", pattern = ".phylip$")
dirlist=list()
for(i in 1:length(gene_to_be_concatted_file_set)){
matches[[i]] = gene_to_be_concatted_file_set[[i]] %in% phylip_file_list
 }
 for(i in 1:length(gene_to_be_concatted_file_set)){
if(length(which(matches[[i]]==TRUE)) > 1){
dirlist[[i]] = paste("concat_set","_",i, sep="")
dir.create(dirlist[[i]])

} else {
	print("file_list_too_short")
	print(gene_to_be_concatted_file_set[[i]])
} }

FINAL_gene_to_be_concatted_file_set = list()
for(i in 1:length(gene_to_be_concatted_file_set)){
if(length(which(matches[[i]]==TRUE)) > 1){
	FINAL_gene_to_be_concatted_file_set[[i]] = gene_to_be_concatted_file_set[[i]]
}
}

FINAL_gene_to_be_concatted_file_set[sapply(FINAL_gene_to_be_concatted_file_set, is.null)] = NULL

dir_list = list.dirs(path = ".", recursive = FALSE)
dir_list =dir_list[-c(406:407)] ### ! this is to remove any other directories from the list that may be in your dir if you don't have any other spurious dirs in your lint no need to run this line

for(i in 1:length(FINAL_gene_to_be_concatted_file_set)){
file.copy(FINAL_gene_to_be_concatted_file_set[[i]], dir_list[i])
}

                                       
filePATH_list=list()
for(i in 1:length(dir_list)){
	filePATH_list[[i]] = list.files(path = dir_list[i], pattern="*.phylip", full.names = TRUE, recursive = TRUE)
	
	}
	
	filePATH_list[lapply(filePATH_list,length)>0] 

gene_to_be_concatted_file_set_fin_NEXUS = lapply(filePATH_list, function(x) { x = paste(x, ".nexus", sep="") })




test = unlist(gene_to_be_concatted_file_set_fin_NEXUS)
to_rm = c(".nexus")
#
test = setdiff(test, to_rm)

## write the phylips into nexus for concatenation in PHYLUCE
test_file_to_write = unlist(filePATH_list)
for(i in 1:length(test)){
	onetowrite = read.phy(test_file_to_write[i])
	write.nex(onetowrite , test[i])
}


                                       
### concatinate each nexus in the dirs created with PHYLUCE magic


#R
#setwd("~/annotations/nexus_for_concat/")
library(ips)
nex_files_list = list.files(pattern = "phylip.nexus$", recursive = TRUE)
dir_list = list.dirs(path=".", recursive = FALSE)
dir_list = dir_list[-c(406:407)] ### ! this is to remove any other directories you dont want

cmd =list()
dir_list = gsub("./", "", dir_list)
## change of directoty below as this was run on remote cluster
for(i in 1:length(dir_list)){
	nex = paste("/home/ubuntu/bird_phy_fin/phylip_files/", dir_list[i], sep="") ## you will have to modify to your working directory
    
	outp = paste(nex,"concat_partitionedbygene.phy", sep="")
    subvcmd = paste("phyluce_align_format_nexus_files_for_raxml --alignments", nex, "--output", outp, "--charsets --log-path log", sep=" ")
   subvcmd = gsub("[:.:]", "", subvcmd)

cmd[[i]] = subvcmd
}


## !!!!!!!!!
 ## !!!!!!!!! in terminal shell change ? to N , this is not done in R !!!!!!!!!!!!!!
## !!!!!!!!!
                                       
mkdir log

#### find replace all ? with N for files in dir

for dir in concat_set_*; do
cd "${dir}"

for file in "${dir}"/*.phylip.nexus; do

find *.phylip.nexus -type f -exec sed -i 's/?/N/g' {} \;

cd ..
done

done


#### !!!!!!! ################## back to R ########################## !!!!!!!!!!!!!!

mclapply(cmd, system, mc.cores=getOption("mc.cores", 3))
                                       





```

### 4. **Organizing files for phylogentic analyses**

Now you should have two different sets of directories one with the original phylips and nexus files e.g. **concat_set_99** and another **concat_set_99concat_partitionedbygenephy** with the *concat_set_99.charsets* *concat_set_99.phylip*, files. Now to build the species trees we need to take or make a list of the files that got concatenated and then subtract it from the total list of phylip files. Then take these lists and feed it to your tree building method of choice, here we will use RaxML. 

Below is the shell script used to execute RAxML over files in a directory, it will be executed in a loop, in R

```shell
# !! open a terminal tab
# change path to match your working directory in R you want to have this script in the wd
############ run_RAxML.sh ##############################################
nano run_RAxML.sh
######### paste the below, modify how you like RAxML, especially your path 

#!/bin/bash
cd $(dirname $1)
id=$(basename $(dirname $1))
phy=$(basename $1)
#below are the RAxML commands, when in doubt give full path to where it lives
/home/ubuntu/standard-RAxML-master/raxmlHPC-PTHREADS-AVX2 -f a -m GTRGAMMA -N 100 -x 12345 -p 25258 -n ${id}.best.tre -s $phy -T 6


######### exit nano Ctrl + o, "return", Ctrl + x

chmod 775 run_RAxML.sh
##########################################################
###

```

Now back to R in terminal 

```R
###
phy_files_list = list.files(pattern = "concat_set_.*.phylip$", recursive = TRUE, path = getwd(), full.names = TRUE)

cmd =list()

for(i in 1:length(phy_files_list)){

cmd[[i]] = paste(getwd(), "/run_RAxML.sh ", phy_files_list[i], sep="")

}

# make sure your system has 48 cores or change the 48 to an appropriate number of cores. Some of the concatented sets are long alignments, it may go a little slower for the short alignments using more cores than needed but really slow for the long ones using far too few, it goes quickly enough

final_raxml = mclapply(cmd, system, mc.cores=getOption("mc.cores", 8))  ### 6 X mc.cores 8 = 48 cores




```



Now that we have made the gene trees for the concatenated UCEs in a gene, we need to make the gene trees for the rest plus the originals "those individual UCEs that got concatenated ". 

We will move all the original phylips that got concatenated into their own directories. Using  unix script that can be simply pasted into in the terminal, this is not strictly necessary, I did this as it helped later in comparisons between subsets of these data. 

```shell
## unix terminal
for dir in concat_set_*; do
cd "${dir}"
find . -name "uce-*.phylip" -exec sh -c 'NEWDIR=`basename "$1new" ` ; mkdir "$NEWDIR" ; mv "$1" "$NEWDIR" ' _ {} \;
cd ..
done


```



You should have something like the following in your main directory 

```shell
uce-805.phylip  ##thousands/hundreds of uce phylip files

concat_set_96 ## hundreds of concat_set_# directories
			  ## concat_set_# with multiple dirs inside e.g. uce-#.phylips 
​	   		concat_set_96/uce-1094.phylipnew/ uce-1094.phylip  
​           concat_set_96/uce-1933.phylipnew/ uce-1933.phylip
              ## the nexus files used in PHYLUCE to concatenate
​           concat_set_96/uce-1094.phylip.nexus uce-1933.phylip.nexus 

concat_set_96concat_partitionedbygenephy  # hundereds of directory with concatenated file and character set from PHYLUCE + RAxML results

​		concat_set_96concat_partitionedbygenephy/ concat_set_96.charsets
​		concat_set_96concat_partitionedbygenephy/ concat_set_96.phylip
​		concat_set_96concat_partitionedbygenephy/ RAxML_bestTree.concat_set_96concat_partitionedbygenephy.best.tre
​		concat_set_96concat_partitionedbygenephy/ RAxML_bipartitionsBranchLabels.concat_set_96concat_partitionedbygenephy.best.tre
​		concat_set_96concat_partitionedbygenephy/ RAxML_bipartitions.concat_set_96concat_partitionedbygenephy.best.tre
​		concat_set_96concat_partitionedbygenephy/ RAxML_bootstrap.concat_set_96concat_partitionedbygenephy.best.tre
​		concat_set_96concat_partitionedbygenephy/ RAxML_info.concat_set_96concat_partitionedbygenephy.best.tre


```



Now back to **R**, we will list the phylips in the top directory (the one with the original complete set of phylips), and then subtract those from the ones used in concatenation for example  "**concat_set_96/uce-1094.phylipnew**/*uce-1094.phylip*". 

```R
setwd("/PHYLIPS/birds/test_dir_oct/") ## your wd
new_dir_list = list.dirs(path=".",  recursive = TRUE,  full.names = TRUE)
new_dir_list = new_dir_list[grepl("*.phylipnew", new_dir_list) ]

library(ips)
phy_files_list= list()
for(i in 1:length(new_dir_list)){
phy_files_list[[i]] = list.files(path=new_dir_list[[i]], pattern = "*.phylip$", recursive = TRUE,  full.names = TRUE)
}
phy_files_list = unlist(phy_files_list)
cmd =list()

phy_files_list = gsub(".*phylipnew/", "", phy_files_list)

all_phy_files = list.files(".", pattern = "*.phylip$", recursive = FALSE,  full.names = FALSE)

those_not_concated_by_gene = setdiff(all_phy_files, phy_files_list)
wd = getwd()

those_not_concated_by_gene = paste(wd,"/" ,those_not_concated_by_gene, sep="")

dir.create("those_not_concated_by_gene")

file.copy(those_not_concated_by_gene, "those_not_concated_by_gene/",  )

setwd("/PHYLIPS/birds/test_dir_oct/those_not_concated_by_gene/")

## !!!!!!!
##   in unix
## !!!!!!!!
cd /PHYLIPS/birds/test_dir_oct/those_not_concated_by_gene/
find . -name "uce-*.phylip" -exec sh -c 'NEWDIR=`basename "$1new" ` ; mkdir "$NEWDIR" ; mv "$1" "$NEWDIR" ' _ {} \;


## !!!!!!!
# !!!!!! Copy/paste "run_RAxML.sh" from directoy above, change number of cores and RAxML flavor raxmlHPC-SSE3 or raxmlHPC
## !!!!!!!!


phy_files = list.files(".", pattern = "*.phylip$", recursive = TRUE,  full.names = TRUE)
wd = getwd()
for(i in 1:length(phy_files)){
phy_files[[i]] = gsub("./", "/", phy_files[[i]])
phy_files[[i]] = paste(wd,phy_files[[i]], sep="")
}

for(i in 1:length(phy_files)){
cmd[[i]] = paste(getwd(), "/run_RAxML.sh ", phy_files[[i]], sep="")
}


# make sure your  system has 48 cores or change the 48 to an appropriate number of cores.

library(parallel)
final_raxml = mclapply(cmd, system, mc.cores=getOption("mc.cores", 48))  ### 48 cores






```



### 5. **Round up the bootstrap files paths and best tree files for ASTRAL**

```shell
# in unix treminal 

mkdir tree_results
mkdir boot_trees

#find . -type f -iname "RAxML_bestTree*" -exec cp -t tree_results/ {} \;

#mkdir tree_results_mrB_msct

find . -name "RAxML_bootstrap.*.best.tre" -exec cp {} boot_trees/ \;
cd boot_trees
find `pwd` -name "RAxML_bootstrap.*.best.tre" > all.merge.boot_strapped_actualbstrees.txt
cd ..

find . -name "RAxML_bipartitions.uce-*.best.tre" -exec cp {} tree_results/ \;

cd tree_results
sed -n wall_best_tree.bird.merge.tre *.tre


# copy files to ASTRAL dir
cp all_best_tree.bird.merge.tre ~/ASTRAL-master/ #or where ever you have ASTRAL
cd ..
cd boot_trees

cp all.merge.boot_strapped_actualbstrees.txt ~/ASTRAL-master/

#run ASTRAL how you wish
java -jar astral.5.6.1.jar -i all_best_tree.bird.merge.tre -b all.merge.boot_strapped_actualbstrees.txt -r 100 -o all.bird.reg.astral.speciestr_out.tre


```







  

