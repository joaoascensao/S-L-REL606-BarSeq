'''
7.1.19
Joao Ascensao

functions for gene ontology enrichment
'''
import pandas as pd
import scipy as sp 
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import convert_Ecoli_geneIDs
from matplotlib.backends.backend_pdf import PdfPages
from goatools.obo_parser import GODag
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.base import download_ncbi_associations
from goatools.test_data.genes_NCBI_10090_ProteinCoding import GENEID2NT as GeneID2nt_mus

sns.set(style="whitegrid")


def run_GOEA(id_list,strain,cat):
	fin_gene2go = download_ncbi_associations()



	obodag = GODag("go-basic.obo")

	# Read NCBI's gene2go. Store annotations in a list of namedtuples
	objanno = Gene2GoReader(fin_gene2go, taxids=[511145])


	ns2assoc = objanno.get_ns2assc()

	for nspc, id2gos in ns2assoc.items():
	    print("{NS} {N:,} annotated E. coli genes".format(NS=nspc, N=len(id2gos)))


	ncbi_ecoli=pd.read_csv('gene_result_NCBI.txt',sep='\t')


	from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
	goeaobj = GOEnrichmentStudyNS(
	        ncbi_ecoli['GeneID'], # List of e coli protein-coding genes
	        ns2assoc, # geneid/GO associations
	        obodag, # Ontologies
	        propagate_counts = True,
	        alpha = 0.05, # default significance cut-off
	        methods = ['fdr_bh']) # defult multipletest correction method


	goea_results_all=goeaobj.run_study(id_list)
	goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]
	print(goea_results_sig)

	goeaobj.wr_xlsx('res/GOEnrichment_{}_{}.xlsx'.format(strain,cat), goea_results_sig)
	#from goatools.godag_plot import plot_gos, plot_results, plot_goid2goobj
	#plot_results('GOEnrichment_%s_%s_{NS}.png'%(strain,cat), goea_results_sig)
	



