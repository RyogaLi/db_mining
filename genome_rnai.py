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
		self._url = "http://www.genomernai.org/rest/genes/entrezId/"

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
		info = self._url + gene
		r = requests.get(info, "html.parser")
		data = str(r.text)
		data = json.loads(data)
		return data

	def _get_experiment(self, expr):
		"""
		expr: expr id 
		return: info for this experiment (id, gene, biomodel, pubmed)
		"""

def main():
	# init
	rnai = genomeRNA()
	rnai._get_gene("25925")


if __name__ == '__main__':
	# main()
	f = "./GenomeRNAi_v17_Homo_sapiens.txt"
	main()