# Created 2018-May-17
# Author: Roujia Li
###
# Using REST api to grep information from Ensmbl db (GRCh38 build) - current build
# Detailed description of REST api: http://grch37.rest.ensembl.org/
# Main purpose of this api:
# 1. convert gene ids (from ensmbl to entrez)

from __future__ import division
import pandas as pd
import requests, sys
import json
import ast
import time
import mechanize
from bs4 import BeautifulSoup


class Ensmbl(object):

	def __init__(self):
		self._name = "Ensmbl, GRCh37"
		self._grch37 = "http://grch37.rest.ensembl.org/"
		self._url = "https://rest.ensembl.org/"

	def _retrieve_entrez_id(self, gene):
		"""
		Using ensmbl current build
		gene: ensembl gene id
		return: entrez gene id
		doc: http://grch37.rest.ensembl.org/documentation/info/xref_id
		"""
		ext = "/xrefs/id/"+gene+"?"+"external_db=EntrezGene;all_levels=1"
		try:
			r = requests.get(self._url+ext, headers={ "Content-Type" : "application/json"})
		except Exception:
			print "Exception:" + str(r.status_code)
			time.sleep(60)
			r = requests.get(self._url+ext, headers={ "Content-Type" : "application/json"})
		
		if r.status_code != 200:
			print "Exception:" + str(r.status_code)
			time.sleep(60)
			r = requests.get(self._url+ext, headers={ "Content-Type" : "application/json"})
		
		data = str(r.text)
		data = json.loads(data)
		if len(data) == 1:
			data = data[0]
			entrez_id = data[u"primary_id"]
		elif len(data) > 1:
			entrez_id  = "more than one entrez"
		else:
			entrez_id = "N/A"

		return str(entrez_id)


def main():
	# init 
	f = "./test_all_genes.csv"
	table = pd.read_table(f, sep=",")
	ensg = table.ensg.tolist()
	entrez_id = {}
	start = time.time()
	for i in ensg:
		en = Ensmbl()
		try:
			entrez = en._retrieve_entrez_id(i)
		except Exception:
			entrez = "failed"

		entrez_id[i] = entrez

	end = time.time()
	print "processing time:" + str(end - start)

	id_convert = pd.DataFrame(entrez_id.items())
	id_convert.to_csv("id_convert.csv", index=False)

if __name__ == '__main__':
	main()