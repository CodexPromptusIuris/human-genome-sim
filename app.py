import streamlit as st
from genetic_engine import DNASequence, GeneticLab, BioConstants

# Configuración de la página
st.set_page_config(page_title="Genomic Lab Simulator", page_icon="🧬", layout="wide")

# Estilos CSS personalizados para visualizar ADN
st.markdown("""
<style>
    .base-A { color: #FF5733; font-weight: bold; }
    .base-T { color: #3375FF; font-weight: bold; }
    .base-C { color: #FFC300; font-weight: bold; }
    .base-G { color: #28B463; font-weight: bold; }
    .protein { font-family: 'Courier New'; color: #8E44AD; font-weight: bold;}
</style>
""", unsafe_allow_html=True)

def colorize_dna(sequence):
    """Genera HTML para colorear las bases."""
    html = ""
    for base in sequence:
        html += f"<span class='base-{base}'>{base}</span>"
    return html

# --- HEADER ---
st.title("🧬 Simulador Genómico Humano")
st.markdown("Plataforma de investigación in-silico para análisis y edición genética (CRISPR/SNP).")

# --- SIDEBAR: CARGA DE DATOS ---
with st.sidebar:
    st.header("1. Cargar Muestra")
    data_source = st.radio("Fuente de ADN:", ["Secuencia Manual", "Archivo FASTA"])
    
    input_seq = ""
    seq_name = "Muestra_Manual"

    if data_source == "Secuencia Manual":
        input_seq = st.text_area("Introduce secuencia (5' -> 3'):", "ATGCGTAGCTAG")
    else:
        uploaded_file = st.file_uploader("Sube archivo .fasta", type=["fasta", "txt"])
        if uploaded_file is not None:
            stringio = uploaded_file.getvalue().decode("utf-8")
            seq_name, input_seq = GeneticLab.read_fasta(stringio)
            st.success(f"Cargado: {seq_name}")

    # Botón para inicializar/resetear
    if st.button("Inicializar Laboratorio"):
        try:
            dna = DNASequence(seq_name, input_seq)
            st.session_state['dna_object'] = dna
            st.session_state['history'] = ["Secuencia Original Cargada"]
            st.success("ADN cargado correctamente.")
        except Exception as e:
            st.error(f"Error: {e}")

# --- AREA PRINCIPAL ---
if 'dna_object' in st.session_state:
    dna = st.session_state['dna_object']

    # Panel de Información
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Longitud (pb)", len(dna.strand_5_3))
    with col2:
        gc_content = GeneticLab.analyze_gc_content(dna)
        st.metric("Contenido GC", f"{gc_content}%")
    with col3:
        st.metric("ID Muestra", dna.id)

    st.divider()

    # Visualizador de Secuencia
    st.subheader("🔬 Visualización de la Doble Hélice")
    
    # Mostramos solo los primeros 100 caracteres para no saturar
    display_limit = 60
    st.markdown(f"**5'** {colorize_dna(dna.strand_5_3[:display_limit])} ... **3'**", unsafe_allow_html=True)
    st.text(f"   {'|' * len(dna.strand_5_3[:display_limit])}")
    st.markdown(f"**3'** {colorize_dna(dna.strand_3_5[:display_limit])} ... **5'**", unsafe_allow_html=True)

    st.divider()

    # Traducción
    st.subheader("🧪 Expresión de Proteínas (Traducción)")
    protein = GeneticLab.synthesize_protein(dna)
    st.markdown(f"**Cadena de Aminoácidos:**")
    st.markdown(f"<p class='protein'>{protein}</p>", unsafe_allow_html=True)

    st.divider()

    # --- ZONA DE EXPERIMENTACIÓN ---
    st.subheader("⚠️ Zona de Edición Genética")
    
    tab1, tab2 = st.tabs(["Mutación Puntual (SNP)", "Reporte"])

    with tab1:
        c1, c2, c3 = st.columns([1, 1, 1])
        with c1:
            pos = st.number_input("Posición (Índice)", min_value=0, max_value=len(dna.strand_5_3)-1, value=0)
        with c2:
            base_new = st.selectbox("Nueva Base", BioConstants.BASES)
        with c3:
            st.write("") # Espaciador
            st.write("") 
            if st.button("Aplicar Mutación"):
                old_base = dna.strand_5_3[pos]
                dna = GeneticLab.mutate_snp(dna, pos, base_new)
                st.session_state['dna_object'] = dna # Actualizar estado
                st.session_state['history'].append(f"Mutación en pos {pos}: {old_base} -> {base_new}")
                st.rerun() # Recargar app para ver cambios
    
    with tab2:
        st.write("Historial de cambios en esta sesión:")
        for item in st.session_state['history']:
            st.write(f"- {item}")

else:
    st.info("👈 Por favor, carga una secuencia de ADN en la barra lateral para comenzar.")
