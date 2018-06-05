# parse OMIM db files
# 1.4 How are mutations cataloged in OMIM?
# Mutations are cataloged in OMIM in the Allelic Variants section of 
# gene entries (see 1.2). For most genes, only selected mutations are 
# included. Criteria for inclusion include the first mutation to be discovered, 
# high population frequency, distinctive phenotype, historic significance, 
# unusual mechanism of mutation, unusual pathogenetic mechanism, 
# and distinctive inheritance (e.g., dominant with some mutations, 
# recessive with other mutations in the same gene). 
# Most of the allelic variants represent disease-causing mutations. A few polymorphisms are included, many of which show a positive correlation with particular common disorders.
# These criteria are spelled out in more detail in: Amberger, J. et al. Nucleic Acids Research 43:D789, 2015 (PMID 25428349).
from __future__ import division
import pandas as pd
import requests
import json
from bs4 import BeautifulSoup


class ParseOMIM(object):

	def __init__(self, omim_gene_map2, ):
		self._genemap2 = pd.read_table(omim_gene_map2)
		# requests by mim number
		self._api_server = "https://api.omim.org/api/entry?mimNumber="

	def _parse_genemap2(self):
		# get only dominant diseases
		self._genemap2 = self._genemap2.dropna(subset=["Phenotypes"])
		dom = self._genemap2[self._genemap2["Phenotypes"].str.contains("dominant")]
		return dom

	def _api_call(self, gene, entry):

		# what to include: all the variants and gene map 
		include = "&include=allelicVariantList&include=geneMap"
		apikey = "&apiKey=dbRVd26jQAaFhCSncTi0iA"
		mimnumber = str(int(gene))
		url = self._api_server+ mimnumber + include + apikey
		r = requests.get(url)

		if r.status_code != 200:
			print "Exception:" + str(r.status_code)
			time.sleep(60)
			r = requests.get(url)

		soup = BeautifulSoup(r.text, 'html.parser')
		# check prefix 
		try:
			prefix = str(soup.find_all("prefix")[0])
			prefix = prefix.split("</")[0][-1]
		except Exception:
			prefix = "N/A"
		if "^" in prefix: 
			# removed entries
			exit(0)

		# find all phenotype associated with this entry
		phenotypemimnumber = soup.find("phenotypemaplist")
		inheritance = soup.find_all("phenotypeinheritance")
		# if "dominant" in inheritance[0]
		# print len(soup.find_all("allelicvariant"))
		# print phenotypemimnumber
		# print soup.find_all("allelicvariant")[0]

		va_list = soup.find_all("allelicvariant")
		if len(va_list) != 0:
			for v in va_list:
				info = v.text.split("\n")
				phenotype = v.find_all("name")[0].text
				if "RECESSIVE" in phenotype: # only select those with dominant phenotype
					continue

				if prefix not in "*#+":
					continue
				# fixed index entries
				entry["phenotype"].append(phenotype)
				entry["variantmimnumber"].append(mimnumber+str(int(info[1])/10000)[1:])
				entry["mimnumber"].append(mimnumber)
				entry["prefix"].append(prefix)

				# mutation info contains gene_name and mutations
				mutation_info = v.find_all("mutations")
				if mutation_info != []:
					m = mutation_info[0].text
					mutation = ", ".join(m.split(", ")[1:])
					entry["mutations"].append(mutation)
				
					gene_name = m.split(", ")[0]
					entry["gene_name"].append(gene_name)
				else:
					entry["mutations"].append("N/A")
					entry["gene_name"].append("N/A")
				
				# not all entry has external ids
				# find these ids and N/A if not found
				if v.find_all("dbsnps") != []:
					dbsnps_id = str(v.find_all("dbsnps")[0]).split("<dbsnps>")[1].split("</dbsnps>")[0]
					entry["dbsnp"].append(dbsnps_id)
				else:
					entry["dbsnp"].append("N/A")

				if v.find_all("clinvaraccessions") != []:
					dbsnps_id = str(v.find_all("clinvaraccessions")[0]).split("<clinvaraccessions>")[1].split("</clinvaraccessions>")[0]
					entry["clinvaraccessions"].append(dbsnps_id)
				else:
					entry["clinvaraccessions"].append("N/A")

				if v.find_all("exacdbsnps") != []:
					dbsnps_id = str(v.find_all("exacdbsnps")[0]).split("<exacdbsnps>")[1].split("</exacdbsnps>")[0]
					entry["exacdbsnps"].append(dbsnps_id)
				else:
					entry["exacdbsnps"].append("N/A")

		# print pd.DataFrame(entry)

		# output should contain:
		# Mimnumber (with prefix), preferredtitle, phenotypemimnumber, phenotype, phenotypeinheritance, variantmimnumber, mutations, name, dbsnps, clinvaraccessions


def main():
	# init with genemap2 file
	# genemap2 = "./OMIM/genemap2.txt"
	# omim = ParseOMIM(genemap2)

	# # get all the genes for dominant disease
	# dom_genemap = omim._parse_genemap2()
	# # save this to file
	# dom_genemap.to_csv("filtered_omim_genemap2.txt", sep="\t", index=False)
	
	# entry = {"gene_name":[],"mimnumber":[], "prefix":[], 
	# 	"phenotype":[], "variantmimnumber":[], "mutations":[], 
	# 	"dbsnp":[], "clinvaraccessions":[], "exacdbsnps":[]}

	# # get all the mutations in these genes
	# gene_list = dom_genemap["Mim Number"].tolist()
	# for gene in gene_list:
	# 	omim._api_call(gene, entry)

	# df = pd.DataFrame(entry)
	# df.to_csv("test_entry_0524.csv", sep=",", index=False)

	# organize the df so it can be used to cross reference RCVaccession in clinvar db
	omim_entries = "./test_entry_0524.csv"
	clinvar = "./variant_summary.txt"
	df = pd.read_csv(omim_entries)
	df["RCVaccession"] = df["clinvaraccessions"]
	df.drop("clinvaraccessions",axis=1)
	df["RCVaccession"] = df['RCVaccession'].str.replace(';;;',';')
	# print df

	clin = pd.read_table(clinvar)

	# merge based on ID
	s1 = pd.merge(clin, df, how='inner', on=["RCVaccession"])
	
	combined = {"ClinVar_AlleleID": s1["#AlleleID"], "ClinVarName": s1.Name, "ClinVar_gene_symbol": s1.GeneSymbol, "ClinVar_ClinicalSignificance":s1.ClinicalSignificance,
			"RCVaccession": s1.RCVaccession, "dbsnp": s1.dbsnp, "OMIM_gene_name":s1.gene_name, "MimNumber":s1.mimnumber, "OMIM_mutations":s1.mutations,
			"OMIM_phenotype": s1.phenotype, "OMIM_prefix":s1.prefix, "OMIM_variantmimnumber": s1.variantmimnumber}

	df = pd.DataFrame(combined)
	df.to_csv("merged_OMIM_clinVar_0524.csv", index=False)

if __name__ == '__main__':
	main()
