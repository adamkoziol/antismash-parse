#!/usr/bin/env python 3
from glob import glob
from Bio import SeqIO
from accessoryFunctions.accessoryFunctions import *

__author__ = 'adamkoziol'


class ParseResults(object):

    def parse(self):
        """
        Parse the .xml file and extract protein sequence for each hit name provided
        """
        import xml.etree.ElementTree as ElementTree
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        printtime('Parsing biosynML.xml file', self.start)
        # Load the XML file using ElementTree
        xml = ElementTree.parse(self.biosyn)
        # Get the root of the XML tree for parsing the file
        root_element = xml.getroot()
        # Iterate through the tree
        for child in root_element:
            # Only interested in the 'genelist' section of the file
            if child.tag == 'genelist':
                for domain in child:
                    # Initialise a variable to store the gene name
                    genename = str()
                    # Iterate through all the children of this domain
                    for tags in domain:
                        # Set the name of the gene based on the 'gene_name' tag's.text
                        if tags.tag == 'gene_name':
                            genename = tags.text
                        # Only proceed if the contig name is present in the list of contigs
                        if genename in self.contigs:
                            # Search for 'gene_qualifiers'
                            if tags.tag == 'gene_qualifiers':
                                # Iterate through all the qualifiers present in the 'gene_qualifier' section
                                for qual in tags:
                                    # The sequence information is stored within the qualifiers with 'translation' 'name'
                                    # and 'genbank' 'style' attributes
                                    if qual.attrib['name'] == 'translation' and qual.attrib['style'] == 'genbank':
                                        # Create a FASTA-formatted sequence output of the query sequence
                                        # The qualifier.text is the sequence of the hit
                                        record = SeqRecord(Seq(qual.text),
                                                           id=genename,
                                                           description='')
                                        # Do not create multiple identical records - use the name of the contig to
                                        # ensure duplicates are not entered into the list
                                        if genename not in self.contignames:
                                            self.sequences.append(record)
                                            self.contignames.append(genename)
        # Write the list of records to file
        with open(os.path.join(self.multifasta), 'w') as fasta:
            SeqIO.write(self.sequences, fasta, 'fasta')

    def blast(self):
        """
        Perform remote BLAST requests for each protein sequence created above
        """
        from Bio.Blast import NCBIWWW
        printtime('Running remote BLAST analyses', self.start)
        make_path(self.blastpath)
        self.records = list(SeqIO.parse(self.multifasta, "fasta"))
        # Iterate through all the records
        for record in self.records:
            # Add all the names of the ORFs of interest to a set
            self.ctg.add(record.id)
            # Create the name of the
            outputfile = os.path.join(self.blastpath, record.id + '.xml')
            # Create a list of all the BLAST output files
            self.blastoutputs.append(outputfile)
            # Allow exceptions to pass - don't want the script to stop if there is a timeout error
            # noinspection PyBroadException
            try:
                if not os.path.isfile(outputfile):
                    # Set up the BLAST parameters
                    wwwblast = NCBIWWW.qblast('blastp',
                                              # Use the 'nr' database
                                              'nr',
                                              record.format('fasta'),
                                              # Only return the top hit
                                              alignments=1,
                                              hitlist_size=1,
                                              # Mildly restrictive evalue
                                              expect='1e-25')
                    # Write the remote BLAST result to the local file
                    with open(outputfile, "w") as out_handle:
                        out_handle.write(wwwblast.read())
                    wwwblast.close()
            except:
                pass
        # Parse the results, and create a report
        self.strainassociator()

    def strainassociator(self):
        """
        Parse the BLAST results, and associate the original strain ID with the ORFs of interest
        """
        from Bio.Blast import NCBIXML
        printtime('Loading .gbk cluster files', self.start)
        # Create a list of the .gbk files from anti-SMASH
        self.clusters = glob(os.path.join(self.clusterpath, '*.gbk'))
        # Iterate through the ORF files to associate the original strain name with the ORF name
        for filename in self.clusters:
            # Open the .gbk file
            with open(filename) as genbank:
                # Use SeqIO to parse the file
                for record in SeqIO.parse(genbank, 'genbank'):
                    # The locus_tag feature stores the ORF name
                    for feature in record.features:
                        try:
                            # Set the ORF name - will trigger the try/except loop if it is missing
                            ctg = feature.qualifiers['locus_tag'][0]
                            # Only populate the dictionaries if the ORF is in the list of ORFs of interest
                            if ctg in self.ctg:
                                # Populate two dictionaries - one with the strain name; split the contig information
                                # from the strain on an underscore. The second dictionary has the full contig name
                                self.strains[feature.qualifiers['locus_tag'][0]] = record.description.split('_')[0]
                                self.contighits[feature.qualifiers['locus_tag'][0]] = record.description
                        except KeyError:
                            pass
        printtime('Creating report', self.start)
        # Open the report
        with open(self.report, 'w') as report:
            # Start the string storing the results with the header information
            string = 'strain\torf\tcontig\tdescription\taccession\n'
            # Iterate through the BLAST outputs to extract the description of the BLAST hit
            for gbk in self.blastoutputs:
                # Open the GenBank file
                with open(gbk) as genbank:
                    # Parse the BLAST result with BioPython
                    blast_record = NCBIXML.read(genbank)
                    # Iterate through the alignments
                    for record in blast_record.alignments:
                        try:
                            # Set the name of the contig
                            orf = blast_record.query
                            # Extract the strain name from the dictionary using the ORF name as the key - triggers
                            # the try/except loop
                            strain = self.strains[orf]
                            # Extract the contig name from the dictionary
                            contig = self.contighits[orf]
                            # Split the hit definition on '>' symbols - save the first part of the string
                            hitdef = record.hit_def.split('>')[0]
                            try:
                                # Set the hit accession as the second part of the split string
                                hitaccession = '>' + record.hit_def.split('>')[1]
                            except IndexError:
                                # If no accession is included (missing a '>' symbol), set the accession to 'NA'
                                hitaccession = 'NA'
                            # Populate the string with the appropriate variables
                            string += '{}\t{}\t{}\t{}\t{}\n'.format(strain, orf, contig, hitdef, hitaccession)
                        except KeyError:
                            pass
            # Write the results string to the report
            report.write(string)

    def __init__(self, args):
        self.path = args.path
        self.file = os.path.join(self.path, args.file)
        # Read in all the accessions from the file
        try:
            with open(self.file, 'r') as contigs:
                self.contigs = set(contigs.read().splitlines())
        except IOError:
            print('Cannot find file of accessions as provided: {}'.format(self.file))
            quit()
        self.biosyn = os.path.join(self.path, 'biosynML.xml')
        assert os.path.isfile(self.biosyn), 'Cannot locate biosynML.xml file in the path. Please ensure that the file' \
                                            'is present'
        self.start = args.start
        self.multifasta = os.path.join(self.path, 'sequences.fasta')
        self.blastpath = os.path.join(self.path, 'outputs')
        self.clusterpath = os.path.join(args.clusterpath, '')
        self.report = os.path.join(self.path, 'report.csv')
        # Initialise class variables
        self.sequences = list()
        self.contignames = list()
        self.records = list()
        self.clusters = list()
        self.blastoutputs = list()
        self.contighits = dict()
        self.strains = dict()
        self.ctg = set()
        # Run the script
        if not os.path.isfile(self.multifasta):
            self.parse()
        else:
            self.blast()

if __name__ == '__main__':
    # Argument parser for user-inputted values, and a nifty help menu
    from argparse import ArgumentParser
    import time
    # Parser for arguments
    parser = ArgumentParser(description='Parses anti-SMASH output file biosyn.xml to extract protein sequences'
                                        'of genes of interest')
    parser.add_argument('path',
                        help='Specify input directory')
    parser.add_argument('-f', '--file',
                        help='Name of file of contig numbers. This is created by running the following command in '
                             'your antiSMASH_run folder: '
                             'find -name "*_gene_info.txt" -exec grep "ctg*" {} \; | cut -f1 | uniq > contigs.txt')
    parser.add_argument('-c', '--clusterpath',
                        required=True,
                        help='Path to .gbk cluster files created by antiSMASH. Required to associate renamed contigs '
                             'with original files')
    # Get the arguments into an object
    arguments = parser.parse_args()

    # Define the start time
    arguments.start = time.time()
    # Run it
    ParseResults(arguments)
    # Print a bold, green exit statement
    print('\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time.time() - arguments.start) + '\033[0m')
