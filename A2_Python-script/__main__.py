from Bio import SeqIO
from datetime import datetime
import os
import re

from Bio.Blast.Applications import NcbiblastnCommandline
from io import StringIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from bird import *
from constants import *

def min_max_finder(ls):
    max_element = max(ls, default="None")
    min_element = min(ls, default="None")
    if max_element == min_element:
        return str(max_element)
    else:
        return str(min_element) + " - " + str(max_element)

def rename_cr(gene_name):
    """Translating the one-letter-code for an.YCR/an.CR/YCR/CR to the full name.

    :param gene_name: Name of the gene as string
    :return: Translated version as string
    """
    if gene_name == name_CR:
        return "CR"
    elif gene_name == name_YCR:
        return "YCR"
    elif gene_name == name_annotated_CR:
        return "an.CR"
    elif gene_name == name_annotated_YCR:
        return "an.YCR"
    else:
        return gene_name

if __name__ == "__main__":
    print("+++ +++ +++ Y-CR SEARCHER +++ +++ +++")

    genbank_input = input("Please input the Genbank input file name: ")
    name_output = input("Please input the output file name prefix: ")

    gb_file = SeqIO.parse('./genbank-input/' + genbank_input, 'gb')

    COMPLETE_BIRDS = []  # List of all full bird genomes objects
    #BIRDS_GENEORDER = set()  # Set of all different Bird gene orders
    current_order = None #For sorting each order in an individual fasta file
    current_genus = None

    # Informations for result paths
    date = datetime.today().strftime('%Y-%m-%d-%H:%M:%S')
    parent_directory = "./results/"
    path = os.path.join(parent_directory, date)
    os.mkdir(path)

    error_file = open('./results/' + date + '/' + name_output + '_ERROR.txt', 'w')
    names_file = open('./results/' + date + '/' + name_output + '_NAMES.txt', 'w')
    analyse_file = open('./results/' + date + '/' + name_output + '_ANALYSE.tsv', 'w')

    analyse_file.write(
        "Accession" + t + "Name" + t + "Taxonomy" + t + "Genome Size" + t + "Genome complete?" + t + "Novel gene order?" + t + "Duplication remnants" + t + "CR length" + t + "YCR length" + t + "CR-YCR Coverage" + t + "CR-YCR E-Value" + t + "CR-YCR Identity" + t + "Number of BLAST results" + t + "Genes \n")

    genus_set = set()

    for gb_record in gb_file:

        gb_record.id = Bird(gb_record)

        # Filter all GenBank entry where the genome size is too small
        if gb_record.id.length < min_len_genome:
            gb_record.id.set_completeness(False)
            error_file.write(str(gb_record.id.accession) + ": Incomplete Genome! Length! \n")
            #continue

        # Setting author, title etc.
        gb_record.id.set_author_co(gb_record.annotations["references"])

        # Search all desired genes
        gb_record.id.gene_searcher(gb_record.features[1:])

        # Filter all bird genomes where not all tRNAs of interest are present
        if gb_record.id.is_present(name_tRNAThr) is False or gb_record.id.is_present(
                name_tRNAPro) is False or gb_record.id.is_present(name_tRNAGlu) is False or gb_record.id.is_present(
            name_tRNAPhe) is False:
            gb_record.id.set_completeness(False)
            error_file.write(str(gb_record.id.accession) + ": Incomplete Genome! tRNAs! \n")
            continue

        # Filter all bird genomes where tRNA-Phe occurs twice
        if gb_record.id.is_doublet(name_tRNAPhe):
            error_file.write(str(gb_record.id.accession) + ": tRNA-Phe occurred twice! \n")
            gb_record.id.set_duplication_remnants(True)
            continue

        # Calibration of all present genes:
        calibration_length = gb_record.id.genome_calibration(gb_record.id.gene(name_tRNAPhe).start)
        if calibration_length != 0:
            phe_start = gb_record.id.gene(name_tRNAPhe).start
            for gene in gb_record.id.genes:
                if gene.name != name_tRNAPhe:
                    gene.update_location(phe_start, gb_record.id.length, calibration_length)
            gb_record.id.gene(name_tRNAPhe).update_location(phe_start, gb_record.id.length, calibration_length)
            # tRNAPhe, as calibration point for all other genes, has to be calibrated last

        # PSEUDO/CONTROL REGION SEARCH:
        # We're searching for all CR candidates and YCR candidates
        before_gene = None
        for gene in gb_record.id.sorted_genes():
            if gene.name == name_annotated_CR or gene.name == name_annotated_YCR or gene.name == name_CR or gene.name == name_YCR or gene.name == name_rRNA12S or gene.name == name_tRNAPhe:
                continue
            elif before_gene is None:
                before_gene = gene
                continue
            else:
                if gene.start - before_gene.end > min_len_cr:
                    gb_record.id.append_gene(name_CR, gb_record.id.full_genome, before_gene.end, gene.start - 1)
                before_gene = gene
        if gb_record.id.length - before_gene.end > min_len_ycr:
            gb_record.id.append_gene(name_YCR, gb_record.id.full_genome, before_gene.end, gb_record.id.length)

        # Storage of all birds of interest
        COMPLETE_BIRDS.append(gb_record.id)

    #All complete entries sorted by taxonomy
    COMPLETE_SORTED_BIRDS = sorted(COMPLETE_BIRDS, key=lambda seq_obj: seq_obj.taxonomy[0])

    #Write to file!
    for gb_record.id in COMPLETE_SORTED_BIRDS:

        # Write to File
        cr_len = 0
        ycr_len = 0

        # Check for duplication remnants:
        gb_record.id.search_duplication_remnants()

        # Check for completeness of GenBank entry
        if "complete genome" in gb_record.id.ncbi_title or "full genome" in gb_record.id.ncbi_title:
            gb_record.id.set_completeness(True)
        elif "partial genome" in gb_record.id.ncbi_title:
            gb_record.id.set_completeness(False)

        if str(current_order) != str(gb_record.id.taxonomy[0]):
            print(str(gb_record.id.taxonomy[0]))
            if current_order is not None:
                cr_file.close()
                ycr_file.close()
            current_order = gb_record.id.taxonomy[0]
            genus_set = set()
            genus_set.add(str(gb_record.id.taxonomy[-1]))
            cr_file = open('./results/' + date + '/' + name_output + '_'  + str(gb_record.id.taxonomy[0]) + '_CR.fasta', 'w')
            ycr_file = open('./results/' + date + '/' + name_output + '_' + str(gb_record.id.taxonomy[0]) + '_YCR.fasta', 'w')
        else:
            genus_set.add(str(gb_record.id.taxonomy[-1]))
            print(str(gb_record.id.taxonomy[-1]))

        # Fasta to File
        if gb_record.id.is_present(name_CR) and gb_record.id.gene(name_CR).length > min_len_cr:
            gb_record.id.set_gene_order(False)
            cr_len = gb_record.id.gene(name_CR).length
            if gb_record.id.is_present(name_YCR):
                gb_record.id.set_gene_order(True)
                ycr_len = gb_record.id.gene(name_YCR).length
                cr_file.write(">" + str(gb_record.id.accession) + ":" + str(gb_record.id.organism) + " "+ str(
                    gb_record.id.taxonomy[0]) +" "+ "Novel gene order CR: \n")  # name
                cr_file.write(str(gb_record.id.gene(name_CR).sequence) + '\n')
                ycr_file.write(">" + str(gb_record.id.accession) + ":" + str(gb_record.id.organism) + " " + str(
                    gb_record.id.taxonomy[0]) +" "+ "Novel gene order YCR: \n")  # name
                ycr_file.write(str(gb_record.id.gene(name_YCR).sequence) + '\n')
        elif gb_record.id.is_present(name_YCR) and gb_record.id.gene(name_YCR).length:
            # Hier werden alle known gene order CR (welche als YCR sortiert wurden) umsortiert
            # Das hier ist unser Bereich für Known gene order CRs die als YCR falsch annotiert wurde
            gb_record.id.set_gene_order(False)
            cr_len = gb_record.id.gene(name_YCR).length
            #cr_file.write(">" + str(gb_record.id.accession) + ";" + str(gb_record.id.organism) + t + str(
                #gb_record.id.taxonomy[0]) + t + "Known gene order CR: \n")  # name
            #cr_file.write(str(gb_record.id.gene(name_YCR).sequence) + '\n')
            if gb_record.id.is_present(name_CR):
                gb_record.id.gene(name_CR).rename("!CR")
            gb_record.id.gene(name_YCR).rename(name_CR)
        else:
            if gb_record.id.is_present(name_YCR):
                ycr_len = gb_record.id.gene(name_YCR).length
            if gb_record.id.is_present(name_CR):
                cr_len = gb_record.id.gene(name_CR).length

        if gb_record.id.novel_order is True:
            seq1 = SeqRecord(Seq(gb_record.id.gene(name_CR).sequence), id="cr")
            seq2 = SeqRecord(Seq(gb_record.id.gene(name_YCR).sequence), id="ycr")
            SeqIO.write(seq1, "seq1.fasta", "fasta")
            SeqIO.write(seq2, "seq2.fasta", "fasta")

            # Run BLAST and parse the output as XML
            output = NcbiblastnCommandline(query="seq2.fasta", subject="seq1.fasta", outfmt=5)()[0]
            blast_result_record = NCBIXML.read(StringIO(output))

            for alignment in blast_result_record.alignments:
                #print(str(alignment))
                for hsp in alignment.hsps:
                    gb_record.id.append_blast([str(hsp.align_length / blast_result_record.query_length), str(hsp.expect), str(hsp.identities / hsp.align_length), 1])
                    #print("#########################")
                    #print(str(gb_record.id.accession))
                    #print(str(hsp.align_length / blast_result_record.query_length), str(hsp.expect), str(hsp.identities / hsp.align_length), 1)
                    #print(hsp.query)
                    #print(hsp.match)
                    #print(hsp.sbjct)

     
        # Gene order to file
        analyse_file.write(
            gb_record.id.accession + t + str(gb_record.id.organism) + t + str(gb_record.id.taxonomy) + t + str(
                gb_record.id.length) + t + str(gb_record.id.completeness) + t + str(gb_record.id.novel_order) + t + str(gb_record.id.duplication_remnants) + t + str(cr_len) + t + str(ycr_len) + t + gb_record.id.print_blast() + t)
        for g in gb_record.id.sorted_genes():
            if g.name != "!CR" and g.name != "!YCR":
                analyse_file.write(str(rename_cr(g.name)) + "(" + str(g.start) + "-" + str(g.end) + ")" + t)
        analyse_file.write("\n")

        # Names to file
        names_file.write(
            str(gb_record.id.accession) + " " + str(gb_record.id.organism) + "\t" + str(
                gb_record.id.organism) + " " + str(gb_record.id.accession) + "\n")



    analyse_file.close()
    names_file.close()
    cr_file.close()
    ycr_file.close()
    error_file.close()

    print("+++ +++ +++ Done! Thanks for using Y-CR searcher! +++ +++ +++")
