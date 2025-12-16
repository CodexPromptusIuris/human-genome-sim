import streamlit as st
import textwrap
from Bio import Entrez

# ==========================================
# 1. MOTOR GENÉTICO (BACKEND LÓGICO)
# ==========================================

class BioConstants:
    """Constantes universales del dogma central."""
    BASES = ['A', 'T', 'C', 'G']
    PAIRS = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    # Tabla de Codones Estándar (ADN -> Aminoácido)
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

class DNASequence:
    """Objeto que representa la doble hélice."""
    def __init__(self, sequence_id, sequence_str):
        self.id = sequence_id
        # Limpieza: Solo mayúsculas y quitamos espacios/saltos de línea
        self.strand_5_3 = sequence_str.upper().replace('\n', '').replace(' ', '').strip()
        self.validate_sequence()
        self.strand_3_5 = self._generate_complementary()
        
    def validate_sequence(self):
        # Filtramos caracteres extraños silenciosamente o lanzamos error
        valid_seq = ""
        for base in self.strand_5_3:
            if base in BioConstants.BASES:
                valid_seq += base
        # Si queremos ser estrictos, descomenta la siguiente línea:
        # if len(valid_seq) != len(self.strand_5_3): raise ValueError("Bases inválidas detectadas")
        self.strand_5_3 = valid_seq

    def _generate_complementary(self):
        return ''.join([BioConstants.PAIRS[base] for base in self.strand_5_3])

    def get_codons(self):
        return textwrap.wrap(self.strand_5_3, 3)

class NCBIConnector:
    """Conector para descargar genes reales."""
    @staticmethod
    def fetch_gene(gene_name):
        Entrez.email = "estudiante@genolab.com" # Necesario para NCBI
        try:
            # Buscamos ARNm humano para evitar intrones y facilitar traducción
            query = f"{gene_name}[Gene Name] AND Homo sapiens[Organism] AND biomol_mrna[PROP]"
            
            # 1. Buscar ID
            handle = Entrez.esearch(db="nucleotide", term=query, retmax=1, sort="relevance")
            record = Entrez.read(handle)
            handle.close()
            
            if not record["IdList"]:
                return None, "No se encontró el gen o ARNm adecuado."
            
            gi_id = record["IdList"][0]
            
            # 2. Descargar Secuencia FASTA
            fetch_handle = Entrez.efetch(db="nucleotide", id=gi_id, rettype="fasta", retmode="text")
            fasta_data = fetch_handle.read()
            fetch_handle.close()
            
            # Parseo simple del encabezado
            header = fasta_data.split('\n')[0]
            return fasta_data, f"Descarga Exitosa: {header[:30]}..."
            
        except Exception as e:
            return None, f"Error de conexión: {str(e)}"

class GeneticLab:
    """Herramientas de análisis y manipulación."""
    @staticmethod
    def read_fasta(file_content):
        lines = file_content.strip().split('\n')
        seq_id = "Desconocido"
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
                aa = BioConstants.CODON_TABLE.get(codon, '?')
                if aa == '_':
                    protein_chain.append("*") # Stop
                    break
                protein_chain.append(aa)
        return "-".join(protein_chain)

    @staticmethod
    def analyze_gc_content(dna_obj):
        g = dna_obj.strand_5_3.count('G')
        c = dna_obj.strand_5_3.count('C')
        total = len(dna_obj.strand_5_3)
        return round(((g + c) / total) * 100, 2) if total > 0 else 0

    @staticmethod
    def mutate_snp(dna_obj, position, new_base):
        if position < 0 or position >= len(dna_obj.strand_5_3):
            return dna_obj # No hacer nada si está fuera de rango
        seq_list = list(dna_obj.strand_5_3)
        seq_list[position] = new_base.upper()
        # Reconstruimos el objeto con la nueva secuencia
        return DNASequence(dna_obj.id, "".join(seq_list))

# ==========================================
# 2. INTERFAZ GRÁFICA (STREAMLIT FRONTEND)
# ==========================================

def main():
    st.set_page_config(page_title="Genomic Lab", page_icon="🧬", layout="wide")

    # Estilos CSS para el ADN
    st.markdown("""
    <style>
        .base-A { color: #e74c3c; font-weight: bold; }
        .base-T { color: #3498db; font-weight: bold; }
        .base-C { color: #f1c40f; font-weight: bold; }
        .base-G { color: #2ecc71; font-weight: bold; }
        .protein { font-family: 'Courier New', monospace; color: #8e44ad; font-weight: bold; word-break: break-all;}
        .dna-seq { font-family: 'Courier New', monospace; word-break: break-all; font-size: 14px; }
    </style>
    """, unsafe_allow_html=True)

    def colorize_dna(sequence):
        html = ""
        for base in sequence:
            html += f"<span class='base-{base}'>{base}</span>"
        return html

    st.title("🧬 Laboratorio de Simulación Genética")
    st.markdown("Investigación In-Silico: Edición genética, traducción y análisis NCBI.")

    # --- BARRA LATERAL (CONTROLES) ---
    with st.sidebar:
        st.header("1. Cargar Muestra")
        mode = st.radio("Fuente:", ["Manual", "Archivo FASTA", "Descarga NCBI 🌐"])
        
        seq_name = "Muestra_Manual"
        input_seq = ""

        if mode == "Manual":
            input_seq = st.text_area("Secuencia (5'->3'):", "ATGCGTAGCTAG")
        
        elif mode == "Archivo FASTA":
            f = st.file_uploader("Sube .fasta", type=["fasta", "txt"])
            if f:
                content = f.getvalue().decode("utf-8")
                seq_name, input_seq = GeneticLab.read_fasta(content)
        
        elif mode == "Descarga NCBI 🌐":
            st.info("Busca genes oficiales (ej: MYH7, INS, TP53)")
            gene = st.text_input("Nombre del Gen:", "MYH7")
            if st.button("📡 Buscar en NCBI"):
                with st.spinner("Conectando con NIH..."):
                    data, status = NCBIConnector.fetch_gene(gene)
                    if data:
                        n, s = GeneticLab.read_fasta(data)
                        st.session_state['temp_ncbi'] = (n, s)
                        st.success(status)
                    else:
                        st.error(status)
            
            if 'temp_ncbi' in st.session_state:
                seq_name, input_seq = st.session_state['temp_ncbi']
                st.caption(f"Listo para cargar: {seq_name[:20]}...")

        st.divider()
        if st.button("🧪 Inicializar Laboratorio"):
            if input_seq:
                try:
                    dna = DNASequence(seq_name, input_seq)
                    st.session_state['dna'] = dna
                    st.session_state['log'] = [f"Carga inicial: {seq_name} ({len(input_seq)} pb)"]
                    st.success("ADN inyectado correctamente.")
                except Exception as e:
                    st.error(f"Error: {e}")
            else:
                st.warning("No hay secuencia para cargar.")

    # --- PANEL PRINCIPAL ---
    if 'dna' in st.session_state:
        dna = st.session_state['dna']

        # Métricas
        c1, c2, c3 = st.columns(3)
        c1.metric("Longitud (pb)", len(dna.strand_5_3))
        c2.metric("Contenido GC", f"{GeneticLab.analyze_gc_content(dna)}%")
        c3.metric("ID", dna.id[:20])

        st.divider()

        # Visualizador ADN
        st.subheader("🔬 Estructura Molecular")
        limit = 300 # Límite visual para no saturar navegador
        
        st.markdown(f"<div class='dna-seq'>5' {colorize_dna(dna.strand_5_3[:limit])} ... 3'</div>", unsafe_allow_html=True)
        st.markdown(f"<div class='dna-seq'>&nbsp;&nbsp;&nbsp;{'|' * min(len(dna.strand_5_3), limit)}</div>", unsafe_allow_html=True)
        st.markdown(f"<div class='dna-seq'>3' {colorize_dna(dna.strand_3_5[:limit])} ... 5'</div>", unsafe_allow_html=True)
        
        if len(dna.strand_5_3) > limit:
            st.caption(f"*Visualización truncada a las primeras {limit} bases por rendimiento.*")

        st.divider()

        # Traducción
        st.subheader("🧬 Proteína Resultante (Traducción)")
        protein = GeneticLab.synthesize_protein(dna)
        st.markdown(f"<div class='protein'>{protein[:500]} {'...' if len(protein)>500 else ''}</div>", unsafe_allow_html=True)

        st.divider()

        # Experimentos
        st.subheader("⚠️ Zona de Edición (CRISPR / Mutagénesis)")
        
        tab1, tab2 = st.tabs(["Mutación Puntual", "Historial de Cambios"])
        
        with tab1:
            col_a, col_b, col_c = st.columns([1,1,1])
            pos = col_a.number_input("Posición (Índice 0):", 0, len(dna.strand_5_3)-1, 0)
            base_new = col_b.selectbox("Nueva Base:", BioConstants.BASES)
            
            if col_c.button("Aplicar Mutación"):
                old_base = dna.strand_5_3[pos]
                # Ejecutar mutación
                new_dna = GeneticLab.mutate_snp(dna, pos, base_new)
                
                # Actualizar Estado
                st.session_state['dna'] = new_dna
                st.session_state['log'].append(f"Mutación en pos {pos}: {old_base} -> {base_new}")
                st.rerun()

        with tab2:
            for log_entry in st.session_state['log']:
                st.text(f"- {log_entry}")

    else:
        st.info("👈 Selecciona una fuente en la barra lateral y pulsa 'Inicializar Laboratorio' para comenzar.")

if __name__ == "__main__":
    main()
