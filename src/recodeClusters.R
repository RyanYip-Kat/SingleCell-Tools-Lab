library(monocle3)
library(stringr)
library(dplyr)

cds=readRDS("ips-result/monocle.rds")
rm_cluster=c(13,29,32,15,22)
cds=cds[,!colData(cds)$gp_cluster%in%rm_cluster]
cell_type=recode(colData(cds)$gp_cluster,
		 "1"="ips stem cells",
		 "14"="ips stem cells",
		 "25"="ips stem cells",
		 "26"="ips stem cells",
		 "4"="neuro-ectodermal-like cells",
		 "6"="neuro-ectodermal-like cells",
		 "7"="neuro-epithelium",
		 "8"="neuro-epithelium",
		 "18"="neuro-epithelium",
		 "3"="neuroectoderm-Neuroepithelium",
		 "10"="neuroectoderm-Neuroepithelium",
		 "28"="Krt8+ isl1+ transitional stem cell",
		 "2"="retinal progenitor cells",
		 "9"="retinal progenitor cells",
		 "11"="middle retinal progenitor",
		 "16"="middle retinal progenitor",
		 "17"="mesenchyme",
		 "19"="RGC progenitor",
		 "24"="mitotic cells",
		 "30"="mesenchyme",
		 "33"="RGC progenitor",
		 "34"="RGC progenitor",
		 "5"="late retinal progenitor",
		 "12"="late retinal progenitor",
		 "20"="photoreceptor progenitor",
		 "21"="RGC progenitor",
		 "23"="RGC progenitor",
		 "27"="photoreceptor progenitor",
		 "31"="photoreceptor progenitor",
		 "3"="mesenchyme and RPE progenitor")

seurat=readRDS("ips-result/seurat.rds")
Idents(seurat)=seurat$gp_cluster
seurat=subset(seurat,idents=rm_cluster,invert=TRUE)
seurat=RenameIdents(seurat,
		    "1"="ips stem cells",
                 "14"="ips stem cells",
		 "25"="ips stem cells",
                 "26"="ips stem cells",
                 "4"="neuro-ectodermal-like cells",
                 "6"="neuro-ectodermal-like cells",
                 "7"="neuro-epithelium",
                 "8"="neuro-epithelium",
                 "18"="neuro-epithelium",
                 "3"="neuroectoderm-Neuroepithelium",
                 "10"="neuroectoderm-Neuroepithelium",
                 "28"="Krt8+ isl1+ transitional stem cell",
                 "2"="retinal progenitor cells",
                 "9"="retinal progenitor cells",
                 "11"="middle retinal progenitor",
                 "16"="middle retinal progenitor",
                 "17"="mesenchyme",
                 "19"="RGC progenitor",
                 "24"="mitotic cells",
                 "30"="mesenchyme",
                 "33"="RGC progenitor",
                 "34"="RGC progenitor",
                 "5"="late retinal progenitor",
                 "12"="late retinal progenitor",
                 "20"="photoreceptor progenitor",
                 "21"="RGC progenitor",
                 "23"="RGC progenitor",
                 "27"="photoreceptor progenitor",
                 "31"="photoreceptor progenitor",
                 "3"="mesenchyme and RPE progenitor")


Idents(seurat)=seurat$celltype
seurat=RenameIdents(seurat,"neuro-ectodermal-like cells"="neuroectodermal-like cells",
		    "neuroectoderm-Neuroepithelium"="neuroectodermal-like cells",
		    "middle retinal progenitor"="retinal progenitor cells",
		    "mitotic cells"="retinal progenitor cells",
		    "mesenchyme and RPE progenitor"="mesenchymal cells",
		    "mesenchyme"="mesenchymal cells",
		    "Krt8+ isl1+ transitional stem cell"="mesenchymal cells",
		    "late retinal progenitor"="retinal progenitor cells")

