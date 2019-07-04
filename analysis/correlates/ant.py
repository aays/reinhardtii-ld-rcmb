'''
Early draft of a parser (modeled heavily after pyVCF) to read in a custom annotation table.

Uses generator objects + lazy loading of records.

Usage examples:
import ant
parser = ant.Reader('file.txt.gz') # file should have equivalent .tbi index
chr15_genic = [r for r in parser.fetch('chromosome_15') if r.is_genic]

mt_locus = parser.fetch('chromosome_6', start = 3e5, end = 8e5)
mt_locus_genic = len([r for r in mt_locus if r.is_genic]
mt_locus_intergenic = len([r for r in mt_locus if r.is_intergenic])
gen_to_intergen = mt_locus_genic / mt_locus_intergenic
print(gen_to_intergen)

Quick and dirty way to get all available methods/attributes of a record:
parser = ant.Reader('file.txt.gz')
dir(next(parser.fetch('chromosome_1'))

AH - 10/2017
'''

import gzip
import ast

try:
    import pysam
except ImportError:
    pysam = None

class _Record(object):
    '''A record object - stores all information at a single row in the annotation table.
    '''
    def __init__(self, chromosome, position, reference_base, genic, exonic, intronic, intergenic, utr5,
        utr3, fold0, fold4, fold2, fold3, CDS, mRNA, rRNA, tRNA, feature_names, feature_types,
        feature_ID, cds_position, strand, frame, codon, aa, degen, FPKM, rho, FAIRE, recombination, mutability, quebec_alleles):

        def type_make(ob, ob_type):
            if ob_type == 'bool':
                if ob in [0, '0', '.']:
                    return False
                elif ob == 1 or ob == '1':
                    return True
            elif ob_type == 'int':
                if ob == '.':
                    ob = 'NA'
                else:
                    try:
                        ob = int(ob) # redundant?
                    except ValueError:
                        ob = 'NA'
                return ob
            elif ob_type == 'float':
                try:
                    ob = float(ob)
                except ValueError:
                    ob = 'NA'
                return ob
            elif ob_type == 'str':
                try:
                    ob = str(ob)
                except:
                    ob = 'NA'
                return ob

        self.chrom = type_make(chromosome, 'str')
        self.pos = type_make(position, 'int')
        self.ref = type_make(reference_base, 'str')
        self.is_genic = type_make(genic, 'bool')
        self.is_exonic = type_make(exonic, 'bool')
        self.is_intronic = type_make(intronic, 'bool')
        self.is_intergenic = type_make(intergenic, 'bool')
        self.is_utr5 = type_make(utr5, 'bool')
        self.is_utr3 = type_make(utr3, 'bool')
        self.is_fold0 = type_make(fold0, 'bool')
        self.is_fold4 = type_make(fold4, 'bool')
        self.is_fold2 = type_make(fold2, 'bool')
        self.is_fold3 = type_make(fold3, 'bool')
        self.is_in_CDS = type_make(CDS, 'bool')
        self.is_in_mRNA = type_make(mRNA, 'bool')
        self.is_rRNA = type_make(rRNA, 'bool')
        self.is_tRNA = type_make(tRNA, 'bool')
        self.feature_names = ast.literal_eval(feature_names)
        self.feature_types = ast.literal_eval(feature_types)
        self.feature_ID = type_make(feature_ID, 'str')
        self.cds_position = type_make(cds_position, 'int')
        self.strand = type_make(strand, 'str')
        self.frame = type_make(frame, 'int')
        self.codon = type_make(codon, 'str')
        self.aa = type_make(aa, 'bool')
        self.degeneracy = type_make(degen, 'int')
        self.FPKM = type_make(FPKM, 'float')
        self.rho = type_make(rho, 'float')
        self.FAIRE = type_make(FAIRE, 'float')
        self.map_rho = type_make(recombination, 'float')
        self.mutability = type_make(mutability, 'float')
        self.quebec_alleles = list(quebec_alleles)


class Reader(object):
    '''The actual parser.
    
    Usage: 
    import ant
    parser = ant.Reader([annotation table filename])
    records = [r for r in parser] 
    
    If file is compressed + tabix-indexed, can also use .fetch():
    chr1_start = [r for r in parser.fetch('chromosome_1', start = 0, end = 50)]
    
    Fetch uses tabix's half-open (?) indexing. (ie start = 0 and end = 3 corresponds to getting 1, 2, and 3).
    
    Can set raw = True when initializing parser object to get raw lines instead of _Record objects.
    '''
    def __init__(self, filename = None, compressed = None, raw = False):
        
        super(Reader, self).__init__
        
        if not filename:
            raise Exception('Error: filename not provided.')
        elif filename:
            if compressed is None:
                compressed = filename.endswith('.gz')
            self._reader = open(filename, 'rb' if compressed else 'rt')
        self.filename = filename
        if compressed:
            self._reader = gzip.GzipFile(fileobj = self._reader)
            self.reader = (line.decode('utf-8') for line in self._reader) # gzipped obj returns lines as bytes objs
            self._tabix = pysam.TabixFile(self.filename)
        else:
            self.reader = (line for line in self._reader) # init generator
            self._tabix = None

        # 'burn' header from generator + set aside if user needs
        line = next(self.reader)
        header = []
        header.append(line)
        while line.startswith('##'):
            line = next(self.reader)
            if line.startswith('#chromosome'):
                break
            else:
                header.append(line)

        assert line.startswith('#chromosome') # make sure header has been completely read in

        def _line_to_rec(line):
            '''Converts lines in annotation table to Record objects.
            Helper function to ensure easy fetching of actual Records and not just
            tab-split lines.
            '''
            row = line.rstrip().split('\t')
            assert len(row) == 32
            
            chromosome, position, reference_base, genic, exonic, intronic, intergenic, \
            utr5, utr3, fold0, fold4, fold2, fold3, CDS, mRNA, rRNA, tRNA, feature_names, \
            feature_types, feature_ID, cds_position, strand, frame, codon, aa, degen, \
            FPKM, rho, FAIRE, recombination, mutability, quebec_alleles = row # unpacking list

            record = _Record(chromosome, position, reference_base, genic, exonic, intronic, intergenic, utr5,
            utr3, fold0, fold4, fold2, fold3, CDS, mRNA, rRNA, tRNA, feature_names, feature_types,
            feature_ID, cds_position, strand, frame, codon, aa, degen, FPKM, rho, FAIRE, recombination, mutability, quebec_alleles)

            return record

        # generator without header
        if not raw:
            self.reader = (_line_to_rec(line) for line in self.reader)
        elif raw:
            self.reader = (line for line in self.reader) # lines as strings, not records

        self.cols = line.split('#')[1].split('\t') # get column names
        self.header = list(header)

    def __iter__(self):
        return self.reader

    def metadata(self): # parser.metadata returns column names
        return self.cols

    def head(self): # parser.head returns entire header
        return self.header

    def next(self): # very similar to _line_to_rec above
        '''Return next record in file.'''
        line = next(self.reader)
        row = line.rstrip().split('\t')
        assert len(row) == 32

        chromosome, position, reference_base, genic, exonic, intronic, intergenic, \
        utr5, utr3, fold0, fold4, fold2, fold3, CDS, mRNA, rRNA, tRNA, feature_names, \
        feature_types, feature_ID, cds_position, strand, frame, codon, aa, degen, \
        FPKM, rho, FAIRE, recombination, mutability, quebec_alleles = row # unpacking list

        record = _Record(chromosome, position, reference_base, genic, exonic, intronic, intergenic, utr5,
        utr3, fold0, fold4, fold2, fold3, CDS, mRNA, rRNA, tRNA, feature_names, feature_types, feature_ID, cds_position, strand,
        frame, codon, aa, degen, FPKM, rho, FAIRE, recombination, mutability, quebec_alleles) # most args I've ever written...

        return record
    
    def fetch(self, chrom, start = None, end = None, raw = False):
        '''Returns an iterable of _Record instances. Tabix file needs to have been made using
        the vcf preset. If raw = True, returns raw lines.'''
        
        if not pysam:
            raise Exception('Error: pysam not installed.')
        if not self._tabix:
            self._tabix = pysam.TabixFile(self.filename)
        self.reader = self._tabix.fetch(chrom, start, end)
        
        def _line_to_rec(line):
            '''Repeated from above. I know this is very inefficient and that there's probably a better way
            to be doing this, but I don't know what that better way is for the life of me.
            '''
            row = line.rstrip().split('\t')
            assert len(row) == 32

            chromosome, position, reference_base, genic, exonic, intronic, intergenic, \
            utr5, utr3, fold0, fold4, fold2, fold3, CDS, mRNA, rRNA, tRNA, feature_names, \
            feature_types, feature_ID, cds_position, strand, frame, codon, aa, degen, \
            FPKM, rho, FAIRE, recombination, mutability, quebec_alleles = row # unpacking list

            record = _Record(chromosome, position, reference_base, genic, exonic, intronic, intergenic, utr5,
            utr3, fold0, fold4, fold2, fold3, CDS, mRNA, rRNA, tRNA, feature_names, feature_types, feature_ID, cds_position, 
            strand, frame, codon, aa, degen, FPKM, rho, FAIRE, recombination, mutability, quebec_alleles)

            return record
        
        if not raw:
            self.reader = (_line_to_rec(line) for line in self.reader)
        elif raw:
            self.reader = (line for line in self.reader)
        
        return self.reader
    
        
