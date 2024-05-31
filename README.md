# Amplicon sequence data processing

## Description
These scripts are for processing amplicon sequencing data (16S, 18S, ITS and rbcL) for the EA biofilm project following the [DADA2 workflow](https://benjjneb.github.io/dada2/index.html) and the [diatom DADA2 pipeline](https://github.com/fkeck/DADA2_diatoms_pipeline/).
A script including preliminary analyses for the 16S dataset is included as an example.

## R packages
The `dada2` package can be installed from Bioconductor

```r
if (!requireNamespace("BiocManager", quietly = TRUE))  
  install.packages("BiocManager")  
BiocManager::install("dada2", version = "3.16")
```

The `diatbarcode` package can be installed from GitHub

```r
if (!requireNamespace("devtools", quietly = TRUE))  
  install.packages("devtools")  
devtools::install_github("fkeck/diatbarcode")
```

## Databases
Taxonomy assignment is performed using the following databases
- [SILVA v138.1](https://www.arb-silva.de/) reference database for 16S sequences
- [PR2 v5.0.0](https://pr2-database.org/) reference database for 18S sequences
- [UNITE all eukaryotes v9.0](https://unite.ut.ee/) reference database for ITS sequences
- [Diat.barcode v10](https://github.com/fkeck/diatbarcode?tab=readme-ov-file) reference database for rbcL sequences

## References

QUAST, C., PRUESSE, E., YILMAZ, P., GERKEN, J., SCHWEER, T., YARZA, P., PEPLIES, J. AND GLÖCKNER, F.O. (2012) The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Nucleic Acids Research, 41 (D1), D590-D596. https://doi.org/10.1093/nar/gks1219

GUILLOU, L., BACHAR, D., AUDIC, S., BASS, D., BERNEY, C., BITTNER, L., BOUTTE, C., BURGAUD, G., DE VARGAS, C., DECELLE, J., DEL CAMPO, J., DOLAN, J.R., DUNTHORN, M., EDVARDSEN, B., HOLZMANN, M., KOOISTRA, W.H.C.F., LARA, E., LE BESCOT, N., LOGARES, R., MAHÉ, F., MASSANA, R., MONTRESOR, M., MORARD, R., NOT, F., PAWLOWSKI, J., PROBERT, I., SAUVADET, A., SIANO, R., STOECK, T., VAULOT, D., ZIMMERMANN, P. AND CHRISTEN, R. (2012) The Protist Ribosomal Reference database (PR2): a catalogue of unicellular eukaryote small sub-unit rRNA sequences with curated taxonomy. Nucleic Acids Research, 41 (D1), D597-D604. https://doi.org/10.1093/nar/gks1160

ABARENKOV, K., ZIRK, A., PIIRMANN, T., PÖHÖNEN, R., IVANOV, F., NILSSON, R.H., KÕLJALG, U. (2023): UNITE general FASTA release for eukaryotes. UNITE Community. https://doi.org/10.15156/BIO/2938069

RIMET, F., GUSEV, E., KAHLERT, M., KELLY, M.G., KULIKOVSKIY, M., MALTSEV, Y., MANN, D.G., PFANNKUCHEN, M., TROBAJO, R., VASSELON, V. AND ZIMMERMANN, J. (2019). Diat.barcode, an open-access curated barcode library for diatoms. Scientific Reports, 9 (1), 15116. https://doi.org/10.1038/s41598-019-51500-6

KECK, F., RIMET, F., VASSELON, V. AND BOUCHEZ, A. (2019) A ready-to-use database for DADA2: Diat.barcode_rbcL_312bp_DADA2. Portail Data Inra, V1. https://doi.org/10.15454/HNI1EK
