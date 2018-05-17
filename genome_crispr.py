from __future__ import division
import pandas as pd
import requests, sys
import json
import ast
import time
from BeautifulSoup import BeautifulSoup

class genomeCrispr(object):

	def __init__(self):
		self._name = "genomeCrispr"
		self._url = "http://genomecrispr.dkfz.de/api/sgrnas/ensembl"

	def _grep_gene(self, ensg):
		"""
		send request to server - grep 
		"""
		try:
			r = requests.post(self._url, data = {'query':ensg})
		except Exception:
			print "requests error"
			print r.status_code
			time.sleep(240)
			r = requests.post(self._url, data = {'query':ensg})
			# proc = 0
			# while proc == 0:
			# 	r = requests.post(self._url, data = {'query':ensg})
			# 	if r.ok:
			# 		proc = 1
			# 	else:
			# 		print "Connection lost"
			# 		time.sleep(240)
		if r.status_code == 502 or not r.ok:
			print r.status_code
			time.sleep(240)
			r = requests.post(self._url, data = {'query':ensg})

		data = str(r.text)
		data = json.loads(data)
		new_list = []
		for i in data: # filter = only take hit == true
			if i[u'hit']==u"true":
				new_list.append(i)

		df = pd.DataFrame.from_dict(new_list, dtype=str)

		if not df.empty:
			df = df[['symbol','ensg','cellline','pubmed', 'screentype', 'condition']]
			df = df.drop_duplicates(subset =["cellline", "screentype"], keep='last')

		return df


def read_file(db):
	table = pd.read_table(db, sep=",")
	# table = table.head(300)
	genes = set(table["ensg"].tolist())
	
	return list(genes)


def read_names(ensg):

	names = pd.read_table(ensg, sep=",")
	ensg = names["genes"].tolist()

	return ensg


def convert_crispr(df, ensg, combined):
	"""
	conver the dataframe to dictionary with ensg as key and merge all the entries
	"""
	# combined = {ensg:{}}
	# df = df.drop_duplicates(subset["cellline", "screentype"], take_last=True)
	gene_name = df.symbol.drop_duplicates().tolist()[0]
	cellline = " | ".join(df.cellline.tolist())
	pubmed = " | ".join(df.pubmed.tolist())
	screentype = " | ".join(df.screentype.tolist())
	condition = " | ".join(df.condition.tolist())
	combined[ensg] = {"gene_name": str(gene_name), "cellline": str(cellline), 
			"pubmed":str(pubmed), "screentype": str(screentype), "condition": str(condition)}
	return combined


def test_gene_name(csv_f, gene_name):
	"""
	see if two different entry with the same gene name leads to the same results
	"""

	table = pd.read_table(csv_f, sep=",")
	entries = table[table.symbol == gene_name]
	l = entries.cell_line.tolist()
	l = [i.split(" | ") for i in l]
	print len(l[0])
	print len(l[1])


def main():
	db = "./GenomeCRISPR_full05112017.csv"
	ensg = "./en_names.csv"
	# en_names = pd.DataFrame({"genes": gene_list})
	# en_names.to_csv("en_names.csv", index=False)

	# gene_list = read_file(db)
	l = read_names(ensg)

	print "total genes: "+str(len(l))

	gC = genomeCrispr()

	n = 60 
	gene_list = [l[i:i + n] for i in xrange(0, len(l), n)]
	print "total batches:"+str(len(gene_list))
	bn = 1
	for batch in gene_list:
		start = time.time()
		combined = {}
		for gene in batch:
			gene_sum = gC._grep_gene(gene)
			# gene_sum.to_csv('test_all_genes.csv', mode='a', header=False, index=False)
			if not gene_sum.empty:
				combined = convert_crispr(gene_sum, gene, combined)

		df = pd.DataFrame.from_dict(combined, orient='index')
		df["ensg"] = df.index
		cols = ["ensg", "gene_name", "cellline", "pubmed", "screentype", "condition"]
		df = df[cols]
		df.to_csv('test_all_genes.csv', mode='a', header=False, index=False)
		end = time.time()
		print "processing time:" + str(end - start) + " for " + str(bn)
		print "====="
		time.sleep(60)
		bn +=1

def test_main():
	csv_f = "./test_all_genes.csv"
	gene_name = "ZNHIT2"
	test_gene_name(csv_f, gene_name)


if __name__ == '__main__':
	test_main()



