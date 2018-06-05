import pandas as pd

class ClinvarParse(object):

	def __init__(self, variants_summary):
		self._summary = variants_summary

	def _parse(self):
		self._var = pd.read_table(self._summary, delimiter="\t")

	
	def _grep_by_gene(self, gene_symbol, type=None):
		variants_in_gene = self._var.loc[self._var['GeneSymbol'] == gene_symbol]
		return variants_in_gene

	def _grep_by_disease_type(self, disease_type):
		
		variants_in_gene = self._var.loc[self._var['PhenotypeList'].str.contains(disease_type)]
		# grep only pathogenic
		variants_in_gene = variants_in_gene[variants_in_gene["ClinicalSignificance"] == "Pathogenic"]
		# grep Ghrc 37
		variants_in_gene = variants_in_gene[variants_in_gene["Assembly"] =="GRCh37"]
		# only snp and deletion
		variants_in_gene = variants_in_gene[(variants_in_gene["Type"] == "single nucleotide variant") | (variants_in_gene["Type"] == "deletion")]
		# ignor +/-
		#variants_in_gene = variants_in_gene[variants_in_gene["Name"].str.contains("+")]
		return variants_in_gene

if __name__ == '__main__':
	summary = "./variant_summary.txt"
	Clinvar = ClinvarParse(summary)
	Clinvar._parse()
	filtered = Clinvar._grep_by_disease_type("dominant")
	# select cols
	filtered = filtered[["#AlleleID", "Type", "Name", "GeneID", "GeneSymbol", "ClinicalSignificance", "PhenotypeList", "ReviewStatus", "OtherIDs"]]
	filtered.to_csv("filtered_ClinVar.csv", index = False)
