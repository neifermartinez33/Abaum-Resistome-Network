# Abaum Resistome Network

Análisis de la red de co-ocurrencia del resistoma de *Acinetobacter baumannii*
bajo un marco One Health.

## Entornos conda

| Archivo | Entorno | Etapa |
|---|---|---|
| envs/qc.yml | abaum_qc | Descarga y QC |
| envs/annotation.yml | abaum_annotation | Anotación ARGs |
| envs/network.yml | abaum_network | Redes y análisis |
| envs/r_analysis.yml | abaum_r | Análisis en R |

## Instalación
```bash
for env in envs/*.yml; do mamba env create -f $env; done
```
