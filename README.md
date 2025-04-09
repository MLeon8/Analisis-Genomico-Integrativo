# Analisis-Genomico-Integrativo

Este informe presenta un análisis exhaustivo de los genes nos3, atg9B y calca en humanos, así como del genoma de la bacteria Pseudomonas aeruginosa (var. PAO1). Se explora la localización genómica, variantes de empalme, expresión tisular y estructuras funcionales asociadas a nos3, incluyendo su relación con el gen atg9B mediante análisis de complementariedad de ARN. También se estudia el gen calca, del cual se generan múltiples proteínas funcionales por splicing alternativo, con énfasis en su conservación evolutiva. Finalmente, se aborda la virulencia bacteriana mediante el análisis del gen exoS de P. aeruginosa, comparándolo con su posible homólogo en E. coli, empleando algoritmos de alineamiento local como Smith-Waterman. Este trabajo integra herramientas bioinformáticas como BLAST, Biopython y RNAhybrid, proporcionando una aproximación práctica a la caracterización funcional de genes en distintos contextos biológicos.

# README.md

```markdown
# Práctica de Laboratorio - Estructura y Función del Genoma

Script para el análisis bioinformático del gen NOS3 y otros elementos genómicos como parte del Doctorado en Biología Computacional.

## Requisitos

- Python 3.7 o superior
- Biopython (`pip install biopython`)
- Conexión a internet (para consultas a NCBI)
- Opcional: API key de NCBI (para mayor cantidad de consultas)

## Instalación

1. Clonar el repositorio o descargar los archivos:
   ```bash
   git clone https://github.com/tuusuario/practica-genoma.git
   cd practica-genoma
   ```

2. Instalar dependencias:
   ```bash
   pip install -r requirements.txt
   ```

## Configuración

Editar el archivo `config.py` con tus datos:
```python
EMAIL = "tu_email@uss.cl"  # Email requerido por NCBI
NCBI_API_KEY = "tu_api_key"  # Opcional pero recomendado
```

## Uso

Ejecutar el script principal:
```bash
python practica_genoma.py
```

Los resultados se guardarán en la carpeta `resultados_lab/` que incluye:
- Secuencias en formato FASTA
- Alineamientos
- Resultados de BLAST
- Archivos temporales para análisis externos

## Análisis Externos Requeridos

1. **RNAhybrid**: Para análisis de complementariedad de mRNA
   - Instalar desde: [https://bibiserv.cebitec.uni-bielefeld.de/rnahybrid](https://bibiserv.cebitec.uni-bielefeld.de/rnahybrid)
   - Usar los archivos generados en `resultados_lab/NOS3_mRNA.fasta` y `resultados_lab/ATG9B_mRNA.fasta`

2. **Visualización de Estructuras 3D**:
   - PyMOL o Chimera para visualizar dominios proteicos

## Generación del Informe LaTeX

Ejecutar:
```bash
pdflatex informe.tex
```

Se requiere una instalación funcional de LaTeX con los paquetes:
- biber
- biblatex
- graphicx
- hyperref


## Estructura de archivos recomendada

```
/practica-genoma
│
├── README.md               # Instrucciones
├── practica_genoma.py      # Script principal
├── config.py               # Configuración
├── informe.tex            # Plantilla de informe
├── requirements.txt        # Dependencias
│
├── /resultados_lab         # Resultados generados
│   ├── secuencias/         # Archivos FASTA
│   ├── alineamientos/      # Resultados de alineamiento
│   └── imagenes/           # Gráficos para el informe
│
└── /informe                # Archivos de LaTeX compilados
    ├── informe.pdf         # Informe final
    └── referencias.bib     # Bibliografía
```

## Para compilar el informe LaTeX

1. Instalar una distribución LaTeX completa (TeX Live, MiKTeX)
2. Ejecutar en orden:
   ```bash
   pdflatex informe.tex
   bibtex informe
   pdflatex informe.tex
   pdflatex informe.tex
   ```

Este conjunto de archivos proporciona una solución completa para:
- Ejecutar los análisis
- Documentar el proceso
- Generar un informe académico profesional
- Organizar los resultados
