# BFOA_Seminario2
Optimización paralelizada del algoritmo de forrajeo bacteriano (BFOA) para resolver problemas de alineamiento múltiple de secuencias biológicas.
# Descripción para el Repositorio GitHub

**Parallel Bacterial Foraging Optimization for Multiple Sequence Alignment**

# Algoritmo BFO Paralelo para Alineamiento de Secuencias

Optimización paralelizada del algoritmo de forrajeo bacteriano (BFOA) para resolver problemas de alineamiento múltiple de secuencias biológicas.

## Características Principales

- 🚀 **Procesamiento paralelo** utilizando multiprocessing
- 🧬 **Optimización para MSA** (Multiple Sequence Alignment)
- ⚡ **Versión optimizada** con ORBC (Operador de Refinamiento Basado en Consenso)
- 📊 **Balance perfecto** entre calidad de solución y costo computacional
- 🧪 **Integración con matrices BLOSUM** para evaluación biológica

## Requisitos

- Python 3.8+
- Bibliotecas:
  - NumPy
  - Biopython (opcional para análisis avanzado)

## Uso Básico

```python
from parallel_BFOA_optimizado import ejecutar_optimizacion

resultados = ejecutar_optimizacion(
    archivo_fasta="secuencias.fasta",
    num_bacterias=10,
    iteraciones=20
)
```

## Estructura del Proyecto

```
bfoa-msa/
├── parallel_BFOA.py             # Algoritmo principal
├── BFOAMejorado.py              # Lógica mejorada
├── bacteria.py                  # Clase base bacteriana
├── evaluadorBlosum.py           # Evaluación con BLOSUM
├── fastaReader.py               # Lector de archivos FASTA
└── examples/                    # Ejemplos de prueba
```

## Contribuciones

¡Bienvenidas! Por favor abra un issue o pull request para:
- 🐛 Reportar errores
- 💡 Sugerir mejoras
- 📚 Añadir documentación

## Licencia

Puedes clonar este repositorio y ejecutar el algoritmo siguiendo las instrucciones en [especificar archivo o carpeta, si es necesario].
```

```markdown
# Parallel BFOA 

Optimización paralela del algoritmo de forrajeo bacteriano para alineamiento múltiple de secuencias. Incluye:

- 🚀 Versión paralelizada con multiprocessing
- 🧬 Operador ORBC para refinamiento inteligente
- ⚡ Balance óptimo calidad/rendimiento
- 📊 Evaluación con matrices BLOSUM

```bash
pip install numpy biopython
python parallel_BFOA.py -f secuencias.fasta
```
