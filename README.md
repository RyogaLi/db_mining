#### ENSEMBL
* Convert ensg to entrez gene id using ENSEMBL api
	* Not all ensg can be found in entrez
	* Some ensg are linked to the same entrez gene ID

---

#### genome_crispr
* Parse the genome CRISPR db and grep all the genes with "hit" in cell lines
	* Use ensg as primary key

#### genome_rnai
* Parse the genome CRISPR db and grep all the genes with "hit" in cell lines
	* Use entrez gene ID as primary key
	* ensg can also be used here

---

#### omim
* Finding all the dominant disease related variants from OMIM

#### clinvar
* Finding all the dominant disease related variants from clinVar