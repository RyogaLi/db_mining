from __future__ import division
import pandas as pd
import requests, sys
import json
import ast
import time
import mechanize
from bs4 import BeautifulSoup

class genomeRNA(object):

	def __init__(self, f=None):
		self._name = "genomeRNAi"
		if f: self._f = pd.read_table(f, sep="\t") # load the data file
		self._url_genes = "http://www.genomernai.org/rest/genes/entrezId/"
		self._url_expr = "http://www.genomernai.org/rest/experiments/"
		self._url_phenotype = "http://www.genomernai.org/ajax/phenotypes/fulltable/entrezId/?"
	
	def _parse_raw_file(self):
		# drop those without any phenotype
		droped = self._f[self._f["Phenotype"]!="none"]
		droped.to_csv("dropped_none.txt", sep="\t", index=False)

	def _get_gene(self, gene):
		"""
		gene: gene with entrez gene id
		return: gene info page
		"""
		# get the gene info page and compare ensmbl gene id
		info = self._url_genes + gene
		r = requests.get(info, "html.parser")
		data = str(r.text)
		data = json.loads(data)
		return data

	def _get_experiment(self, expr):
		"""
		expr: expr id 
		return: info for this experiment (id, gene, biomodel, pubmed)
		"""
		expriment = self._url_expr + expr
		r = requests.get(expriment, "html.parser")
		print r.text
		data = str(r.text)
		# data = json.loads(data)
		output = {"KD_bioModel":str(data[u"bioModel"]), "KD_experimentId": data[u"stableId"], "KD_pubmed":data[u"pubmedId"]}
		return output

	def _get_phenotype(self, gene):
		"""
		get full table of phenotype based on gene id (entrez)
		gene: extrez id
		return: pd.DataFrame with cols: "StableId" || "PubmedId" || "Assay" || "Biomodel" || "Phenotype"
		"""
		colnames = ["StableId", "Screen title", "PubmedId", "geneId", "Assay", "Biomodel",
			"ReagentId", "Scoretype", "Score_cutoff", "Score", "Phenotype", "Comment"]
		param = "{%22findBy%22:{%22entrezId%22:"+gene+"},%22findIsRegex%22:{},%22findRegexParams%22:{},%22findOptions%22:{%22skip%22:0,%22limit%22:5000,%22sort%22:{%22screenTitle%22:1}},%22draw%22:2}"
		r = requests.get(self._url_phenotype+param)
		if r.status_code != 200:
			print "Exception:" + str(r.status_code)
			time.sleep(60)
			r = requests.get(self._url_phenotype+param)
		
		data = r.text.encode('ascii', "xmlcharrefreplace")
		data = json.loads(data)
		col = pd.DataFrame(data).data

		# convert output to pd.DataFrame and select useful info
		summary_df = pd.DataFrame(col.tolist(), index=col.index, columns= colnames)
		# drop none
		# drop no effect
		summary_not_none = summary_df[summary_df.Phenotype != "none"]
		# print summary_not_none
		summary_not_none = summary_not_none[summary_not_none.Phenotype != "no effect"]
		summary_not_none['geneId'] = summary_not_none.geneId.apply(''.join)
		# summary_w_ID = summary_not_none[summary_not_none.geneId == gene]
		
		self._sum_filtered = summary_not_none[["StableId", "PubmedId", "Assay", "Biomodel", "Phenotype"]]
		self._sum_filtered.PubmedId = self._sum_filtered.PubmedId.apply(str).fillna("N/A")
		# self._sum_filtered.PubmedId = self._sum_filtered.PubmedId.str.replace(" ", "N/A")
		# print self._sum_filtered
		return self._sum_filtered

	def _convert_df(self):
		# convert self._sum_filtered to dictionary 
		# with gene name as key
		d = {}
		stableId = " | ".join(self._sum_filtered.StableId.tolist())

		pubmed = " | ".join(self._sum_filtered.PubmedId.tolist())
		assay = " | ".join(self._sum_filtered.Assay.tolist())
		biomodel = " | ".join(self._sum_filtered.Biomodel.tolist())
		phenotype = " | ".join(self._sum_filtered.Phenotype.tolist())
		d = {"KD_StableId": str(stableId), "KD_pubmed": pubmed, "KD_assay": str(assay), "KD_bioModel": str(biomodel), "Phenotype": str(phenotype)}
		
		return d


def main():
	table = pd.read_table(summary_file, sep="\t")
	entrez_id = table.entrez.tolist()
	KD_d = {}
	for gene in entrez_id:
		# init
		rnai = genomeRNA()

		tmp_sum = rnai._get_phenotype(str(gene))
		KD_d[gene] = rnai._convert_df()
		# print gene
	df = pd.DataFrame.from_dict(KD_d, orient="index")
	df.to_csv("rnai_different_id.csv")

def test_one_gene(gene):
	rnai = genomeRNA()
	tmp_sum = rnai._get_phenotype(str(gene))
	print rnai._convert_df()


if __name__ == '__main__':
	# main()
	f = "./GenomeRNAi_v17_Homo_sapiens.txt"
	summary_file = "genome_crispr.txt"
	main()

	# test
	# test_one_gene("19")


