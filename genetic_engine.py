import textwrap

# --- BASE DE DATOS BIOQUÍMICA ---
class BioConstants:
    BASES = ['A', 'T', 'C', 'G']
    PAIRS = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    CODON_TABLE = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }

# --- ESTRUCTURA BIOLÓGICA ---
class DNASequence:
    def __init__(self, sequence_id, sequence_str):
        self.id = sequence_id
        # Limpieza básica: mayúsculas y eliminar saltos de línea
        self.strand_5_3 = sequence_str.upper().replace('\n', '').strip()
        self.validate_sequence()
        self.strand_3_5 = self._generate_complementary()
        
    def validate_sequence(self):
        for base in self.strand_5_3:
            if base not in BioConstants.BASES:
                raise ValueError(f"Base ilegal detectada: {base}")

    def _generate_complementary(self):
        return ''.join([BioConstants.PAIRS[base] for base in self.strand_5_3])

    def get_codons(self):
        return textwrap.wrap(self.strand_5_3, 3)

# --- LABORATORIO DE SIMULACIÓN ---
class GeneticLab:
    @staticmethod
    def read_fasta(file_content):
        """Parsea un archivo FASTA simple (string)."""
        lines = file_content.strip().split('\n')
        seq_id = "Unknown"
        sequence = []
        
        for line in lines:
            line = line.strip()
            if line.startswith(">"):
                seq_id = line[1:]
            else:
                sequence.append(line)
        
        return seq_id, "".join(sequence)

    @staticmethod
    def synthesize_protein(dna_obj):
        codons = dna_obj.get_codons()
        protein_chain = []
        for codon in codons:
            if len(codon) == 3:
                amino_acid = BioConstants.CODON_TABLE.get(codon, '?')
                if amino_acid == '_':
                    protein_chain.append("*") # Stop codon visual
                    break
                protein_chain.append(amino_acid)
        return "-".join(protein_chain)

    @staticmethod
    def analyze_gc_content(dna_obj):
        g_count = dna_obj.strand_5_3.count('G')
        c_count = dna_obj.strand_5_3.count('C')
        total = len(dna_obj.strand_5_3)
        return round(((g_count + c_count) / total) * 100, 2) if total > 0 else 0

    @staticmethod
    def mutate_snp(dna_obj, position, new_base):
        if position < 0 or position >= len(dna_obj.strand_5_3):
            raise IndexError("Posición fuera del rango.")
        
        seq_list = list(dna_obj.strand_5_3)
        seq_list[position] = new_base.upper()
        dna_obj.strand_5_3 = "".join(seq_list)
        dna_obj.strand_3_5 = dna_obj._generate_complementary()
        return dna_obj