## get data from multiple db

from genome_crispr import *
from genome_rnai import *
from ensmbl import *

def main_genomeCRISPR():
	"""
	load data from genomeCrispr file
	based on the ensg in the file
	find genes from their online db
	return: none
	file: creates a file that contains all the genes in genomeCrispr, one gene per row
	"""

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
				combined = convert(gene_sum, gene, combined)

		df = pd.DataFrame.from_dict(combined, orient='index')
		df["ensg"] = df.index
		cols = ["ensg", "gene_name", "cellline", "pubmed", "screentype", "condition"]
		df = df[cols]
		df.to_csv('all_genes.csv', mode='a', header=False, index=False)
		end = time.time()
		print "processing time:" + str(end - start) + " for " + str(bn)
		print "====="
		time.sleep(60)
		bn +=1


if __name__ == '__main__':
	main()