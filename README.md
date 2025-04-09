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
```

# informe.tex

```latex
\documentclass[12pt, a4paper]{article}
\usepackage[utf8]{inputenc}
\usepackage[spanish]{babel}
\usepackage{graphicx}
\usepackage{hyperref}
\usepackage{float}
\usepackage{booktabs}
\usepackage{longtable}

\title{Informe de Práctica: Estructura y Función del Genoma}
\author{Doctorando en Biología Computacional \\ Universidad San Sebastián}
\date{\today}

\begin{document}

\maketitle

\section{Introducción}
Este informe documenta el análisis bioinformático del gen \textit{NOS3} y elementos relacionados como parte de la práctica de laboratorio del curso troncal 2025 del Doctorado en Biología Computacional.

\section{Metodología}
Se utilizaron herramientas bioinformáticas y bases de datos públicas mediante:

\begin{itemize}
\item Consultas a NCBI via Entrez API
\item Alineamientos de secuencias con Biopython
\item Búsqueda de dominios proteicos
\item Análisis de contexto genómico
\end{itemize}

\section{Resultados}

\subsection{Análisis del gen NOS3}

\subsubsection{Ubicación genómica}
El gen \textit{NOS3} se encuentra ubicado en...

\begin{figure}[H]
\centering
\includegraphics[width=0.8\textwidth]{images/nos3_genomic_location.png}
\caption{Ubicación genómica del gen NOS3}
\label{fig:nos3_loc}
\end{figure}

\subsubsection{Secuencias obtenidas}
Se obtuvieron las siguientes secuencias:

\begin{longtable}{lll}
\toprule
Tipo & Accesión & Longitud \\
\midrule
Genómica & NC\_000007.14 & 21,026 bp \\
mRNA & NM\_000603.4 & 4,052 bp \\
\bottomrule
\end{longtable}

\subsubsection{Alineamiento}
El alineamiento entre las secuencias genómicas y de mRNA mostró:

\begin{verbatim}
[Resultados del alineamiento aquí]
\end{verbatim}

\subsection{Variantes de splicing}
Se identificaron X variantes de splicing:

\begin{itemize}
\item Variante 1:...
\item Variante 2:...
\end{itemize}

\subsection{Expresión tisular}
El gen \textit{NOS3} se expresa principalmente en:

\begin{itemize}
\item Tejido endotelial
\item Células musculares lisas
\end{itemize}

\section{Discusión}

Los resultados obtenidos demuestran que... [Análisis crítico de los hallazgos]

\section{Conclusiones}

\begin{itemize}
\item El gen \textit{NOS3} presenta...
\item Se confirmó la relación con \textit{ATG9B}...
\item Las proteínas derivadas de \textit{CALCA}...
\end{itemize}

\section*{Anexos}

\subsection*{Código Fuente}
El código utilizado está disponible en: \url{https://github.com/tuusuario/practica-genoma}

\subsection*{Archivos Generados}
Todos los archivos resultantes se encuentran en la carpeta \texttt{resultados\_lab/}

\end{document}
```

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
