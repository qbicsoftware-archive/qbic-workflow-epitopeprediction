#! /usr/bin/python

import sys
sys.path.append('/lustre_cfc/qbic/chris/Fred2')

import IO
from IO.MartsAdapter import MartsAdapter
from IO.RefSeqAdapter import RefSeqAdapter
from Core.Transcript import Transcript
from Prediction.PSSM import Syfpeithi
from Prediction.NetMHC import NetMHC
from IO.UniProtAdapter import UniProtDB
from Core.Base import AASequence

import time
import timeit
import pickle
import argparse
import datetime


VERSION = 'Version 1.01'

parser = argparse.ArgumentParser(description='Epitope Prediction Pipeline', prog='IRMA')

parser.add_argument(
	'--basedir', '-b',
	required=True,
	help = 'Base directory of files and for output.'
	)
parser.add_argument(
	'--somatic', '-s',
	required=True,
	help = 'File with somatic mutations.'
	)
parser.add_argument(
	'--germline', '-g',
	required=False,
	help = 'File with germline mutations.'
	)
parser.add_argument(
	'--alleles', '-a',
	required=True,
	help = 'File with HLA alleles in official nomenclature (X*xy:vz).'
	)
parser.add_argument(
	'--MHCclass', '-m',
	required=True,
	help = 'I or II for MHC-class-I or -II.'
	)
parser.add_argument(
	'--sampleID', '-d',
	required=True,
	help = 'The ID of the processed sample.'
	)
parser.add_argument(
	'--unmutated_peptides', '-u',
	action="store_true",
	help = 'Include unmutated peptides.'
	)

args = parser.parse_args()

# Starting the actual IRMA functions
mart_db = MartsAdapter()
vl = IO.read_GSvar(args.somatic, args.sampleID + '_s')

if args.germline is None:
	number_of_variants = len(vl)

	for v in vl:
		tmp = v.find_zygosity()
		tmp = v.find_coding()
		tmp = v.find_gene(mart_db) # this could be repl. by get_all_variant_genes

	genes = set(filter(None, [v.gene for v in vl]))

	#~ filter only those with a gene and mutation syntax (coding)
	soon_transcripts = list()
	for g in genes:
		soon_transcripts.append([var for var in vl if var.gene == g and var.coding and 'variant_details' in var.metadata and 'nonsynonymous SNV' in var.metadata['variant_details']])	
else:
	vl_normal = IO.read_GSvar(args.germline, args.sampleID + '_g')

	number_of_variants = len(vl)

	for v in vl+vl_normal:
		tmp = v.find_zygosity()
		tmp = v.find_coding()
		tmp = v.find_gene(mart_db) # this could be repl. by get_all_variant_genes

	# What for ?
	for v in vl_normal:
		v.specific = True

	genes = set(filter(None, [v.gene for v in vl+vl_normal]))

	#~ filter only those with a gene and mutation syntax (coding)
	soon_transcripts = list()
	for g in genes:
		soon_transcripts.append([var for var in vl+vl_normal if var.gene == g and var.coding and 'variant_details' in var.metadata and 'nonsynonymous SNV' in var.metadata['variant_details']])

soon_transcripts = filter(None,soon_transcripts)
refseq_db = RefSeqAdapter('/lustre_cfc/qbic/reference_genomes/IRMA_references/RefSeq_human.protein.faa', 66, '/lustre_cfc/qbic/reference_genomes/IRMA_references/RefSeq_human.rna.fna', 66)

ids = mart_db.get_all_variant_ids(genes=list(genes))

gwt = set() # genes without found transcripts 
transcripts = list()
for st in soon_transcripts:
	found_id = False
	for id in ids:
		pi, ps, ti, ts = [None]*4
		if 'RefSeq mRNA [e.g. NM_001195597]' in id and 'RefSeq Protein ID [e.g. NP_001005353]' in id and st[0].gene == id['UniProt Gene Name']:
			ti = id['RefSeq mRNA [e.g. NM_001195597]']
			pi = id['RefSeq Protein ID [e.g. NP_001005353]']
			found_id = True
		else:
			continue
		#TODO find & document which coding NM_s were not found
		for var in st:
			if ti not in var.coding.keys() and pi not in var.coding.keys():
				continue
		try:
			ts = str(refseq_db.get_transcript_sequence(ti).seq)
			ps = str(refseq_db.get_product_sequence(pi).seq)
			if ts and ps:
				tr = Transcript(ti, ts, pi, ps)
			for var in st:
                        	tr.add_variant(var)
                        	transcripts.append(tr)
		except:
			print 'unsuccessful retrieval',ti,refseq_db.get_transcript_sequence(ti),pi,refseq_db.get_product_sequence(pi)
	if not found_id:
		gwt.add(st[0].gene)
for u in gwt:
	print 'Gene had no associated NM/NP pair: ',u

c1 = 0
c2 = 0
peptides = list()
protein_list = list()
variantMapping = {}
for t in transcripts:
	if t.mrna and t.protein:
		try:
			protein_list = t.to_proteins()
			c1 +=1
		except:
			c2+= 1
			print t,sys.exc_info()[0]
		for p in protein_list:
			if args.MHCclass == 'I':
				peptides += p.unfold(8)
				peptides += p.unfold(9)
				peptides += p.unfold(10)
				peptides += p.unfold(11)
			elif args.MHCclass == 'II':
				peptides += p.unfold(15)
			else:
				print 'Invalid MHC class provided: %s' % args.MHCclass
				sys.exit()

for variantMap in variantMapping:
	variantMapping[variantMap] = sorted(set(variantMapping[variantMap]))
print 'Variantproteins:', c1, 'protein fail:', c2

# Filter self-peptides
up_db  = UniProtDB('sp')
up_db.read_seqs('/lustre_cfc/qbic/reference_genomes/IRMA_references/UniProt_HUMAN_v2014_07.fasta')

selfies = [str(p.seq) for p in peptides if up_db.exists(str(p.seq))]
filtered_peptides = [p for p in peptides if str(p.seq) not in selfies]

# Get unmutated peptides from peptides with included variants. TODO make it nice, include in FRED2 ?
if args.unmutated_peptides:
	unmutated_peptides = list()
	for p in filtered_peptides:
		mutated_string = list(str(p.seq))
		c = p.variants[p.variants.keys()[0]].coding.values()
		for pk in p.variants.keys():
			mutated_string[pk] = p.variants[pk].coding.values()[0].aa_mutation_syntax[2]


		unmutated_string = ''.join(mutated_string)
		unmutated_peptides.append(AASequence(unmutated_string))

al = IO.import_allele_list('%s' % (args.alleles))
alle = [y for x,y in al.iteritems()] 


if args.MHCclass == 'I':

	start = timeit.default_timer()

	ss = Syfpeithi('/lustre_cfc/qbic/reference_genomes/IRMA_references/Syfpeithi-Matrices')
	ss.make_predictions(filtered_peptides,alle)
	if args.unmutated_peptides:
		ss.make_predictions(unmutated_peptides, alle)

	stop = timeit.default_timer()
	print 'Syfpeithi took', stop - start , 'sec'

	start = timeit.default_timer()

	nn = NetMHC('netMHC-3.0','netMHCpan')
	nn.make_predictions(filtered_peptides,alle,'netMHC-3.0')
	if args.unmutated_peptides:
		nn.make_predictions(unmutated_peptides, alle, 'netMHC-3.0')

	stop = timeit.default_timer()
	print 'NetMHC took', stop - start , 'sec'

	start = timeit.default_timer()

	nn = NetMHC('netMHC-3.0','netMHCpan')
	nn.make_predictions(filtered_peptides,alle,'netMHCpan-2.4')
	if args.unmutated_peptides:
		nn.make_predictions(unmutated_peptides,alle,'netMHCpan-2.4')

	stop = timeit.default_timer()
	print 'NetPAN took', stop - start , 'sec'

elif args.MHCclass == 'II':

	start = timeit.default_timer()

	ss = Syfpeithi('/lustre_cfc/qbic/reference_genomes/IRMA_references/Syfpeithi-Matrices')
	ss.make_predictions(filtered_peptides,alle)
	if args.unmutated_peptides:
		ss.make_predictions(unmutated_peptides,alle)

	stop = timeit.default_timer()
	print 'Syfpeithi took', stop - start , 'sec'

	start = timeit.default_timer()

	from Prediction.NetMHC import NetMHC
	nn = NetMHC('netMHC-3.0','netMHCpan','netMHCII','netMHCIIpan')
	nn.make_predictions(filtered_peptides,alle,'netMHCII-2.2')
	if args.unmutated_peptides:
		nn.make_predictions(unmutated_peptides,alle,'netMHCII-2.2')

	stop = timeit.default_timer()
	print 'NetMHCII took', stop - start , 'sec'

	start = timeit.default_timer()

	nn = NetMHC('netMHC-3.0','netMHCpan', 'netMHCII','netMHCIIpan')
	nn.make_predictions(filtered_peptides,alle,'netMHCIIpan-2.0')
	if args.unmutated_peptides:
		nn.make_predictions(unmutated_peptides,alle,'netMHCIIpan-2.0')

	stop = timeit.default_timer()
	print 'NetPANII took', stop - start , 'sec'

report_template = ''' 
###################################################################

		EPITOPE PREDICTION REPORT 

###################################################################

Persons in Charge: Christopher Mohr, Mathias Walzer
Date: %s
Pipeline Version: 1.01
Workflow Version: 1.0

Sample ID: %s

Alleles:
-------------
%s

Used Prediction Methods:
-------------
%s

Binding Assessment Criteria:
-------------
Syfpeithi predictions: prediction score > half max score of corresponding allele
netMHC/netMHCpan predictions: affinity (as IC50 value in nM) <= 500

Additional Steps:
-------------
Filtering of self-peptides (UniProtKB/Swiss-Prot HUMAN.fasta.gz - 9/7/14)

Stats:
-------------
Number of Variants: %s
Number of Peptides: %s
Number of Peptides after Filtering: %s
Number of Predictions: %s
Number of Predicted Binders: %s 
Number of Predicted Non-Binders: %s
Number of Binding Peptides: %s
Number of Non-Binding Peptides: %s 

Contacts:
-------------
mohr@informatik.uni-tuebingen.de
walzer@informatik.uni-tuebingen.de
University of Tuebingen, Applied Bioinformatics,
Center for Bioinformatics, Quantitative Biology Center,
and Dept. of Computer Science,
Sand 14, 72076 Tuebingen, Germany
'''

with open('%s%s_peptides.tsv' % (args.basedir, args.sampleID), 'w') as out_file:
	if args.unmutated_peptides:
		header = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'
                new_line = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.2f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.2f\t%s\t%s\n'
		out_file.write(header % ('GENE', 'POS', 'TRANSCRIPT', 'PROTEIN', 'ALLELE', 'METHOD','SCORE', 'AFFINITY', 'PEPTIDE', 'LENGTH', 'BINDER', 'CODING', 'VARIANT_DETAILS', 'GENOTYPE','UNMUTATED_SCORE', 'UNMUTATED_AFFINITY', 'UNMUATED_PEPTIDE','UNMUTATED_BINDER' ))
	else:
		header = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'
		new_line = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.2f\t%s\t%s\t%s\t%s\t%s\t%s\n'
		out_file.write(header % ('GENE', 'POS', 'TRANSCRIPT', 'PROTEIN', 'ALLELE', 'METHOD', 'SCORE', 'AFFINITY', 'PEPTIDE', 'LENGTH', 'BINDER', 'CODING', 'VARIANT_DETAILS', 'GENOTYPE', ))
	
	num_binder = 0
	num_non_binder = 0
	num_predictions = 0
	num_peptides = 0

	unique_binders = set()
	unique_non_binders = set()
	methods = set()

	for k in range(0, len(filtered_peptides)):
		binder_unique = False
		num_peptides += 1
		gene_name = filtered_peptides[k].variants.itervalues().next().gene
		protein_id = filtered_peptides[k].id.split('|')[1]

		if args.unmutated_peptides:
			unmutated_peptide = str(unmutated_peptides[k].seq)

		transcript_ids = []
		pos = '%s:%s_%s' % (filtered_peptides[k].variants.itervalues().next().chromosome,filtered_peptides[k].variants.itervalues().next().start,
			filtered_peptides[k].variants.itervalues().next().stop)
		zygosity = ''
		variant_details = ''
		
		coding = '%s:' % gene_name
		for variant in filtered_peptides[k].variants.values():
			variant_details = '%s,' % variant.metadata['variant_details'][0]
			transcript_ids.extend(variant.coding.keys())
			zygosity += '%s,' % variant.metadata['genotype'][0]
			for w,v in variant.coding.items():
				coding += '%s:%s:%s,' % (w,v.cds_mutation_syntax,v.aa_mutation_syntax)

		transcript_ids = ','.join(transcript_ids)

		for l in range(0, len(filtered_peptides[k].scores)):

			methods.add(filtered_peptides[k].scores[l].method)
			if args.unmutated_peptides:
				unmutated_score = unmutated_peptides[k].scores[l].score
				unmutated_affinity = unmutated_peptides[k].scores[l].affinity
			num_predictions += 1

			if filtered_peptides[k].scores[l].method == 'Syfpeithi':
				binder = (float(filtered_peptides[k].scores[l].affinity) > 50.0)
			else:
				binder = (float(filtered_peptides[k].scores[l].affinity) <= 500.000)

			if args.unmutated_peptides:
				if unmutated_peptides[k].scores[l].method == 'Syfpeithi':
					unmutated_binder = (float(unmutated_peptides[k].scores[l].affinity) > 50.0)
				else:
					unmutated_binder = (float(unmutated_peptides[k].scores[l].affinity) <= 500.000)

			if binder:
				binder_unique = True
				num_binder +=1
			else:
				num_non_binder += 1

			if args.unmutated_peptides:
				out_file.write(new_line % (gene_name, pos, transcript_ids, protein_id, filtered_peptides[k].scores[l].allele, filtered_peptides[k].scores[l].method, filtered_peptides[k].scores[l].score, float(filtered_peptides[k].scores[l].affinity), filtered_peptides[k].seq, len(filtered_peptides[k].seq), binder, coding[:-1], variant_details, zygosity[:-1],unmutated_score, float(unmutated_affinity), unmutated_peptide, unmutated_binder))
			else:
				out_file.write(new_line % (gene_name, pos, transcript_ids, protein_id, filtered_peptides[k].scores[l].allele, filtered_peptides[k].scores[l].method, filtered_peptides[k].scores[l].score, float(filtered_peptides[k].scores[l].affinity), filtered_peptides[k].seq, len(filtered_peptides[k].seq), binder, coding[:-1], variant_details, zygosity[:-1]))

		if binder_unique:
			unique_binders.add(str(filtered_peptides[k].seq))
		else:
			unique_non_binders.add(str(filtered_peptides[k].seq))


with open('%s%s_prediction_info.info' % (args.basedir, args.sampleID), 'w') as report:
	report.write(report_template % (datetime.datetime.now().strftime("%d/%m/%Y"), args.sampleID, 
		'\n'.join([str(y) for x,y in al.iteritems()] ), '\n'.join(methods), number_of_variants, len(set(peptides)), len(set(filtered_peptides)), num_predictions, 
		num_binder, num_non_binder, len(unique_binders), len(unique_non_binders)))


print 'Number of Variants: %s' % number_of_variants
print 'Number of Peptides: %s' %  len(set(peptides))
print 'Number of Peptides after Filtering: %s' % len(set(filtered_peptides))
print 'Number of Predictions: %s' % num_predictions
print 'Number of Predicted Binders: %s' % num_binder
print 'Number of Predicted Non-Binders: %s' % num_non_binder
print 'Number of Binding Peptides: %s' % len(unique_binders)
print 'Number of Non-Binding Peptides: %s' % len(unique_non_binders)

