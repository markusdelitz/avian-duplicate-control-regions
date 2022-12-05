from constants import *


class Sequence:

    def __init__(self, name, og_full_genome, start, end):
        """Initiation of object Sequence.

        :param name: Name of the gene
        :param og_full_genome: Full sequence of the GenBank entry
        :param start: Start position of the gene
        :param end: End position of the gene
        """
        self.name = name
        self.start = start
        self.end = end
        self.sequence = og_full_genome[self.start:self.end]
        self.length = len(self.sequence)

    def __eq__(self, other):
        """Compares Sequence objects by name, start and end.

        :param other: Sequence object
        :return:
        """
        if isinstance(other, Sequence):
            return self.name == other.name and self.start == other.start and self.end == other.end

    def update_location(self, calibration_point_start, genome_len, calibration_len):
        """Setting a new location concerning the current calibration point (tRNA-Phe) and the calibration length.

        :param calibration_point_start: Start point of the current gene, which is the calibration point
        :param genome_len: Length of the whole genome
        :param calibration_len: Length of the calibration
        :return: Method without return = None
        """
        if self.start >= calibration_point_start:
            self.start = self.start - calibration_len
            self.end = self.end - calibration_len
        else:
            self.start = self.start + (genome_len - calibration_len)
            self.end = self.end + (genome_len - calibration_len)

    def rename(self, name):
        """Renaming the gene

        :param name: Desired new name as string
        :return: Method without return = None
        """
        self.name = name


class Bird:

    def __init__(self, genbank_entry):
        """Initiation of object bird.

        :param genbank_entry: GenBank entry as object
        """
        self.genbank_entry = genbank_entry
        self.ncbi_title = str(genbank_entry.description).lower()
        self.author = None
        self.title = None
        self.journal = None
        self.novel_order = None
        self.duplication_remnants = False
        self.completeness = None
        self.accession = str(genbank_entry.id)
        self.organism = genbank_entry.annotations["organism"]
        self.taxonomy = genbank_entry.annotations["taxonomy"][14:]
        self.og_full_genome = genbank_entry.seq
        self.length = len(self.og_full_genome)
        self.full_genome = ""
        self.blast_results = []
        self.genes = []

    def set_gene_order(self, novel_order):
        self.novel_order = novel_order

    def set_duplication_remnants(self, duplication_remnants):
        self.completeness = duplication_remnants

    def search_duplication_remnants(self):
        if self.is_doublet(name_tRNAGlu) or self.is_doublet(name_tRNAPro) or self.is_doublet(name_tRNAThr) or self.is_doublet(name_ND6) or self.is_doublet(name_ND5) or self.is_doublet(name_CytB) or self.is_doublet(name_rRNA12S):
            self.duplication_remnants = True

    def set_completeness(self, completeness):
        self.completeness = completeness

    def append_blast(self, blast_result):
        if self.blast_results:
            self.blast_results[3] = int(self.blast_results[3]) + 1
            #raise Exception("More than one blast result for a CR YCR pair")
        else:
            self.blast_results = blast_result
        #Vorsicht, manchmal gibts mehr hits!
        #raise Exception("More than one blast result for a CR YCR pair")

    def print_blast(self):
        if self.blast_results:
            return str(self.blast_results[0]) +t+ str(self.blast_results[1]) +t+ str(self.blast_results[2]) +t+ str(self.blast_results[3])
        else:
            return str(None) +t+ str(None) +t+ str(None) +t+ str(0)

    def set_author_co(self, references):
        """Sets information about Author, Title and Journal to the current object.

        :param references: List of references from annotations of GenBank file
        :return: Method without return = None
        """
        for ref in references:
            self.author = ref.authors
            self.title = ref.title
            self.journal = ref.journal

    def genome_calibration(self, calibration_point):
        """We calibrate all our genomes, which don't start with a certain Gene.

        :param calibration_point: Starting position (int) of the gene with which every genome should start
        :return: We return the length of our calibration.
        """
        if calibration_point != 0:
            self.full_genome = str(self.og_full_genome[calibration_point:self.length]) + str(
                self.og_full_genome[0:calibration_point])
            calibration_len = len(self.og_full_genome[0:calibration_point])
            return calibration_len
        else:
            self.full_genome = str(self.og_full_genome)
            return 0

    def append_gene(self, name, sequence, start, end):
        """We append every found gene to the list of genes (self.genes).

        :param name: Name of the gene
        :param sequence: Full genome sequence in which the gene is present
        :param start: Start position of the gene
        :param end: End position of the gene
        :return: Method without return = None
        """
        new_gene = Sequence(name, sequence, start, end)
        if new_gene not in self.genes:
            self.genes.append(new_gene)

    def is_present(self, name_str):
        """We check if a certain gene is present in the list of self.genes.

        :param name_str: Name of the desired gene
        :return: True or False
        """
        if self.genes:
            for gene_obj in self.genes:
                if gene_obj.name == name_str:
                    return True
            return False
        else:
            return False

    def is_doublet(self, name_str):
        """We check if a certain gene occurs more than one in the list of self.genes.

        :param name_str: Name of the desired gene
        :return: True or False
        """
        gene_list = []
        if self.genes:
            for gene_obj in self.genes:
                if gene_obj.name == name_str:
                    gene_list.append(gene_obj)
            if len(gene_list) > 1:
                return True
            else:
                return False
        else:
            return False

    def gene(self, name_str, copy=False):
        """We search for every gene with the desired name in our list of self.genes, add it to a list
         and return the object with the smaller starting position if copy = False or return the object with
         the bigger starting position if copy = True.

        :param name_str: Name of the desired gene
        :param copy: True or False, if we want to obtain the object with the higher or smaller starting position
        :return: Returns the gene as object
        """
        gene_list = []
        for gene_obj in self.genes:
            if gene_obj.name == name_str:
                gene_list.append(gene_obj)
        if not gene_list:
            raise Exception("Gene was not found!")
        gene_list.sort(key=lambda seq_obj: seq_obj.start)
        if copy:
            if len(gene_list) <= 2:
                return gene_list[1]
            else:
                print(str(gene) + "=> present more than twice!")
                return gene_list[1]
                # raise Exception("Gene is present more than twice!")
        else:
            return gene_list[0]

    def sorted_genes(self):
        """Return a sorted list of self.genes, sorted by gene starting positions.

        :return: Sorted List of self.genes by starting position
        """
        return sorted(self.genes, key=lambda seq_obj: seq_obj.start)

    def gene_searcher(self, features):
        """Here we search for our desired genes: tRNAs T, P, E , F, ND5, ND6, 12S-rRNA, CytB and annotated CR/YCR.
        We set the found genes right away in self.object.

        :param features: List of features from the current record of the GenBank file
        :return: Method without return = None
        """
        for feature in features:
            for key, value in feature.qualifiers.items():
                if "trna-thr" in str(value).lower() or "trnt" in str(value).lower():
                    self.append_gene(name_tRNAThr, self.og_full_genome, feature.location.start, feature.location.end)
                if "trna-pro" in str(value).lower() or "trnp" in str(value).lower():
                    self.append_gene(name_tRNAPro, self.og_full_genome, feature.location.start, feature.location.end)
                if "trna-glu" in str(value).lower() or "trne" in str(value).lower():
                    self.append_gene(name_tRNAGlu, self.og_full_genome, feature.location.start, feature.location.end)
                if "trna-phe" in str(value).lower() or "trnf" in str(value).lower():
                    self.append_gene(name_tRNAPhe, self.og_full_genome, feature.location.start, feature.location.end)
                if "nd5" in str(value).lower() or "nadh dehydrogenase subunit 5" in str(
                        value).lower() or "nadh 5" in str(value).lower():
                    self.append_gene(name_ND5, self.og_full_genome, feature.location.start, feature.location.end)
                if "nd6" in str(value).lower() or "nadh dehydrogenase subunit 6" in str(
                        value).lower() or "nadh 6" in str(value).lower():
                    self.append_gene(name_ND6, self.og_full_genome, feature.location.start, feature.location.end)
                if "12s" in str(value).lower():
                    self.append_gene(name_rRNA12S, self.og_full_genome, feature.location.start, feature.location.end)
                if "cytochrome b" in str(value).lower() or "cytb" in str(value).lower() or "cyt b" in str(
                        value).lower():
                    self.append_gene(name_CytB, self.og_full_genome, feature.location.start, feature.location.end)
                if "pseudo" in str(value).lower() or "cr2" == str(value).lower() or "cr 2" == str(
                        value).lower() or "control region 2" in str(value).lower() or "control region II" in str(
                    value).lower():
                    self.append_gene(name_annotated_YCR, self.og_full_genome, feature.location.start,
                                     feature.location.end)
                elif "control" in str(value).lower() or "cr1" == str(value).lower() or "cr 1" == str(value).lower():
                    self.append_gene(name_annotated_CR, self.og_full_genome, feature.location.start,
                                     feature.location.end)
