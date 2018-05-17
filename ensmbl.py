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
		r = requests.get(self._url+ext, headers={ "Content-Type" : "application/json"})
		data = str(r.text)
		data = json.loads(data)
		# print data
		entrez_id = data[u"primary_id"]
		return str(entrez_id)


def main():
	# init 
	en = Ensmbl()
	entrez_id = en._retrieve_entrez_id("ENSG00000198795")


if __name__ == '__main__':
	main()