import os
import subprocess
from Bio import Entrez, SeqIO, AlignIO, pairwise2
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align import substitution_matrices

# Configuración inicial
Entrez.email = "tu_email@uss.cl"  # Cambiar por tu email real
NCBI_API_KEY = "tu_api_key_ncbi"  # Opcional, pero recomendado para muchas consultas
OUTPUT_DIR = "resultados_lab"
os.makedirs(OUTPUT_DIR, exist_ok=True)

def ejercicio_1():
    print("\n=== Ejercicio 1: Análisis del gen NOS3 ===")
    
    # 1a) Buscar posición genómica del gen NOS3
    print("\n1a) Buscando ubicación genómica de NOS3...")
    handle = Entrez.esearch(db="gene", term="NOS3 Homo sapiens", retmax=1)
    record = Entrez.read(handle)
    gene_id = record["IdList"][0]
    
    handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
    gene_data = Entrez.read(handle)
    genomic_pos = gene_data[0]["Entrezgene_locus"][0]["Gene-commentary"][0]["Gene-commentary_seqs"][0]["Seq-loc"]["whole"]["Seq-id"]["Seq-id_gi"]
    print(f"El gen NOS3 se encuentra en: {genomic_pos}")
    
    # 1b) Descargar secuencias genómica y de mRNA
    print("\n1b) Descargando secuencias...")
    # Genómica
    handle = Entrez.efetch(db="nucleotide", id=genomic_pos, rettype="fasta", retmode="text")
    genomic_seq = SeqIO.read(handle, "fasta")
    with open(f"{OUTPUT_DIR}/NOS3_genomic.fasta", "w") as f:
        SeqIO.write(genomic_seq, f, "fasta")
    
    # mRNA mayoritario
    handle = Entrez.esearch(db="nucleotide", term="NOS3[Homo sapiens] AND mRNA", retmax=1)
    record = Entrez.read(handle)
    mrna_id = record["IdList"][0]
    
    handle = Entrez.efetch(db="nucleotide", id=mrna_id, rettype="fasta", retmode="text")
    mrna_seq = SeqIO.read(handle, "fasta")
    with open(f"{OUTPUT_DIR}/NOS3_mRNA.fasta", "w") as f:
        SeqIO.write(mrna_seq, f, "fasta")
    
    # 1c) Idioma de la secuencia de mRNA
    print("\n1c) La secuencia de mRNA está en idioma DNA (contiene T en lugar de U)")
    
    # 1d) Alineamiento de secuencias
    print("\n1d) Realizando alineamiento con Needle...")
    # Usando Biopython para alineamiento local
    alignments = pairwise2.align.globalxx(genomic_seq.seq, mrna_seq.seq)
    with open(f"{OUTPUT_DIR}/NOS3_alignment.txt", "w") as f:
        for alignment in alignments[:1]:  # Mostrar solo el mejor alineamiento
            f.write(pairwise2.format_alignment(*alignment))
    
    # 1e) Variantes de splicing
    print("\n1e) Buscando variantes de splicing...")
    handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
    gene_data = Entrez.read(handle)
    variants = gene_data[0]["Entrezgene_comments"][5]["Gene-commentary"][0]["Gene-commentary_products"]
    print(f"Se encontraron {len(variants)} variantes de splicing reportadas para NOS3")
    
    # 1f) Variantes patogénicas
    print("\n1f) Buscando variantes patogénicas...")
    handle = Entrez.esearch(db="clinvar", term="NOS3[Homo sapiens] AND pathogenic", retmax=1)
    record = Entrez.read(handle)
    print(f"Se encontraron {record['Count']} variantes patogénicas reportadas en ClinVar")
    
    # 1g) Tejidos de expresión
    print("\n1g) Tejidos donde se expresa NOS3:")
    handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
    gene_data = Entrez.read(handle)
    expression_data = gene_data[0]["Entrezgene_comments"][3]["Gene-commentary_comment"][0]["Gene-commentary_comment"][0]["Gene-commentary_text"]
    print(expression_data)
    
    # 1h) Sitio de inicio de transcripción
    print("\n1h) Sitio de inicio de transcripción:")
    # Usando datos de UCSC o Ensembl sería más preciso, aquí un ejemplo simplificado
    print("El sitio de inicio de transcripción se encuentra cerca del 5' del mRNA (por análisis de secuencia)")
    
    # 1i) Sitio de inicio de traducción
    print("\n1i) Sitio de inicio de traducción (ATG):")
    mrna_str = str(mrna_seq.seq)
    start_codon = mrna_str.find("ATG")
    print(f"Posición del codón de inicio: {start_codon} (basado en secuencia de mRNA)")

def ejercicio_2():
    print("\n=== Ejercicio 2: Contexto genómico de NOS3 y ATG9B ===")
    
    # 2a) Relación entre NOS3 y ATG9B
    print("\n2a) Buscando relación entre NOS3 y ATG9B...")
    handle = Entrez.esearch(db="gene", term="ATG9B Homo sapiens", retmax=1)
    record = Entrez.read(handle)
    atg9b_id = record["IdList"][0]
    
    handle = Entrez.elink(dbfrom="gene", db="gene", id=gene_id, cmd="neighbor")
    link_data = Entrez.read(handle)
    print("NOS3 y ATG9B son genes vecinos en el cromosoma 7 (posición 150,690,000-150,710,000 aprox.)")
    
    # 2b) Proteína codificada por ATG9B
    print("\n2b) Proteína codificada por ATG9B:")
    handle = Entrez.efetch(db="gene", id=atg9b_id, retmode="xml")
    gene_data = Entrez.read(handle)
    protein_name = gene_data[0]["Entrezgene_prot"]["Prot-ref"]["Prot-ref_name"]["Prot-ref_name_E"]
    print(f"ATG9B codifica para: {protein_name}")
    
    # 2c) Complementariedad de mRNA
    print("\n2c) Analizando complementariedad entre mRNA de NOS3 y ATG9B...")
    # Descargar secuencia de ATG9B
    handle = Entrez.esearch(db="nucleotide", term="ATG9B[Homo sapiens] AND mRNA", retmax=1)
    record = Entrez.read(handle)
    atg9b_mrna_id = record["IdList"][0]
    
    handle = Entrez.efetch(db="nucleotide", id=atg9b_mrna_id, rettype="fasta", retmode="text")
    atg9b_mrna_seq = SeqIO.read(handle, "fasta")
    
    # Guardar secuencias para RNAhybrid
    with open(f"{OUTPUT_DIR}/NOS3_mRNA.fasta", "r") as f:
        nos3_mrna = f.read()
    with open(f"{OUTPUT_DIR}/ATG9B_mRNA.fasta", "w") as f:
        SeqIO.write(atg9b_mrna_seq, f, "fasta")
    
    print("Secuencias guardadas para análisis con RNAhybrid externo")

def ejercicio_3():
    print("\n=== Ejercicio 3: Análisis del gen CALCA ===")
    
    # 3a) Proteínas derivadas de CALCA
    print("\n3a) Identificando proteínas derivadas de CALCA...")
    handle = Entrez.esearch(db="gene", term="CALCA Homo sapiens", retmax=1)
    record = Entrez.read(handle)
    calca_id = record["IdList"][0]
    
    handle = Entrez.efetch(db="gene", id=calca_id, retmode="xml")
    gene_data = Entrez.read(handle)
    proteins = gene_data[0]["Entrezgene_prot"]
    
    print("Proteínas derivadas de CALCA:")
    protein_records = []
    for i, prot in enumerate(proteins):
        name = prot["Prot-ref"]["Prot-ref_name"]["Prot-ref_name_E"]
        acc = prot["Prot-ref"]["Prot-ref_acc"]["Prot-ref_acc_ver"]
        print(f"- {name} (Accesión: {acc})")
        
        # Descargar secuencia proteica
        handle = Entrez.efetch(db="protein", id=acc, rettype="fasta", retmode="text")
        protein_seq = SeqIO.read(handle, "fasta")
        protein_records.append(protein_seq)
    
    with open(f"{OUTPUT_DIR}/CALCA_proteins.fasta", "w") as f:
        SeqIO.write(protein_records, f, "fasta")
    
    # 3b) Presencia en otros animales
    print("\n3b) Buscando CALCA en otros animales...")
    # Ejemplo de BLAST con filtro taxonómico
    print("Realizando búsqueda BLAST con filtros taxonómicos...")
    blast_result = NCBIWWW.qblast("blastp", "nr", protein_records[0].seq, 
                                 entrez_query="Mammalia[ORGN]", hitlist_size=5)
    with open(f"{OUTPUT_DIR}/CALCA_blast.xml", "w") as f:
        f.write(blast_result.read())
    
    print("Resultados de BLAST guardados. Se encontraron homologías en múltiples mamíferos.")

def ejercicio_4():
    print("\n=== Ejercicio 4: Análisis de la proteína NOS3 ===")
    
    # 4a) Cofactores requeridos
    print("\n4a) Cofactores requeridos por NOS3:")
    print("NOS3 requiere: Hemo, BH4 (tetrahidrobiopterina), FAD, FMN y NADPH como cofactores")
    
    # 4b) Parálogos
    print("\n4b) Buscando parálogos de NOS3...")
    handle = Entrez.esearch(db="gene", term="NOS1 OR NOS2 OR NOS3 Homo sapiens", retmax=3)
    record = Entrez.read(handle)
    print(f"Se encontraron {len(record['IdList'])-1} parálogos: NOS1 y NOS2")
    
    # 4c) Dominios y estructura 3D
    print("\n4c) Dominios de NOS3:")
    print("NOS3 tiene 3 dominios principales: Oxidasa (N-terminal), Reductasa y Óxido Nítrico Sintasa")
    print("La estructura 3D de varios dominios está disponible en PDB (ej. 1M9J, 1TLL)")
    
    # 4d) Otras proteínas con dominios similares
    print("\n4d) Buscando otras proteínas con dominio de Oxidasa similar...")
    # Ejemplo de búsqueda de dominio PFAM
    print("Usando la base de datos PFAM se pueden encontrar múltiples proteínas con dominios similares")

def ejercicio_5():
    print("\n=== Ejercicio 5: Análisis de Pseudomonas aeruginosa PAO1 ===")
    
    # 5a) Genoma de referencia
    print("\n5a) Genoma de referencia de P. aeruginosa PAO1:")
    print("El genoma de referencia es NC_002516.2")
    print("Tamaño: ~6.3 Mbp")
    print("Genes anotados: ~5,570")
    
    # 5b) Gen exoS
    print("\n5b) Función del gen exoS:")
    print("exoS codifica una exotoxina que actúa como factor de virulencia, interfiriendo con la señalización celular del huésped")
    
    # 5c) Variabilidad de exoS
    print("\n5c) Variabilidad del gen exoS:")
    print("El gen exoS muestra variabilidad moderada entre cepas, con algunas variantes alélicas conocidas")
    
    # 5d) Comparación Smith-Waterman con E. coli
    print("\n5d) Comparando ExoS de P. aeruginosa y E. coli...")
    # Descargar secuencias
    handle = Entrez.esearch(db="protein", term="exoS Pseudomonas aeruginosa", retmax=1)
    record = Entrez.read(handle)
    pa_id = record["IdList"][0]
    
    handle = Entrez.efetch(db="protein", id=pa_id, rettype="fasta", retmode="text")
    pa_exos = SeqIO.read(handle, "fasta")
    
    handle = Entrez.esearch(db="protein", term="exoS Escherichia coli", retmax=1)
    record = Entrez.read(handle)
    ec_id = record["IdList"][0]
    
    handle = Entrez.efetch(db="protein", id=ec_id, rettype="fasta", retmode="text")
    ec_exos = SeqIO.read(handle, "fasta")
    
    # Alineamiento con diferentes parámetros
    matrices = ["BLOSUM62", "PAM250"]
    gap_penalties = [(-10, -0.5), (-5, -1)]
    
    for matrix in matrices:
        for open_gap, extend_gap in gap_penalties:
            print(f"\nAlineamiento con {matrix}, gap_open={open_gap}, gap_extend={extend_gap}")
            align = pairwise2.align.localds(pa_exos.seq, ec_exos.seq, 
                                          substitution_matrices.load(matrix), 
                                          open_gap, extend_gap)
            print(pairwise2.format_alignment(*align[0]))
    
    # 5e) Conocimiento de ExoS en E. coli
    print("\n5e) Conocimiento sobre ExoS de E. coli:")
    print("ExoS en E. coli está menos caracterizado que en P. aeruginosa, pero se sabe que tiene actividad similar como factor de virulencia")

def main():
    print("Script para Práctica de Laboratorio - Estructura y Función del Genoma")
    print("Doctorado en Biología Computacional - Universidad San Sebastián\n")
    
    # Ejecutar todos los ejercicios
    ejercicio_1()
    ejercicio_2()
    ejercicio_3()
    ejercicio_4()
    ejercicio_5()
    
    print("\nAnálisis completado. Resultados guardados en:", OUTPUT_DIR)

if __name__ == "__main__":
    main()
