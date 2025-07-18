# LoReC - Comparator

### LoReC – Comparator is a part of LoReC (Long-Read-Checker) toolkit.

LoReC – Comparator tool enables to compare structural variants (SVs) from long-read-sequencing (LRS) datasets obtained by all available technologies as well as short-read sequencing (SRS) datasets and optical genome mapping (OGM). 

The tool enables to analyse targeted or whole-genome sequencing datasets using reference genome of interest. The SVs are compared based on defined parameters such as distance variance, intersection variance and size differences between individual SVs from the different technologies or variant callers. 

It offers several filtering options for particular SVs in base technology/dataset, including the closest SV in compared dataset(s), overlapping genes, variant types, and distance variance, intersection variance and size differences between datasets etc. Additionally, the tool provides possibility to annotate SVs using various databases such as the dbVar_Common, ClinVar or custom annotation file(s).

The tool allows the analysis of multiple input files from the same technology or variant caller at once. Simply enter multiple file paths, separated by semicolons.

This tool is available as a platform-independent CLI application and is part of the Long-Read-Checker (LoReC) toolkit.

## Current stable version
<b>Current stable version is 1.0 - branch 1.0.</b> Main/master branch is dedicated for further development and should not be used in production environment.

## Requirements
Java Runtime Environment 8 or higher.

## Command line arguments
| Parameter | Long                               | Type    | Default | Description                                                                                                                                                                | Required |
|-----------|------------------------------------|---------|---------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------|
| -b        | --bionano_input                    | String  |         | Bionano Genomics analysis pipeline result SMAP file path.                                                                                                                  |          |
| -a        | --annotsv_input                    | String  |         | AnnotSV analysis result TSV file paths delimited by semicolon.                                                                                                             |          |
| -s        | --samplot_input                    | String  |         | Samplot csv variants file paths delimited by semicolon.                                                                                                                    |          |
| -vl       | --vcf_longranger_input             | String  |         | Longranger vcf variants file paths delimited by semicolon.                                                                                                                 |          |
| -vs       | --vcf_sniffles_input               | String  |         | Sniffles vcf variants file paths delimited by semicolon.                                                                                                                   |          |
| -vm       | --vcf_manta_input                  | String  |         | Manta vcf variants file paths delimited by semicolon.                                                                                                                      |          |
| -vi       | --vcf_iclr_input                   | String  |         | Illumina Dragen ICLR wgs vcf variants file paths delimited by semicolon.                                                                                                   |          |
| -vc       | --vcf_clinvar_input                | String  |         | Clinvar vcf variants file paths delimited by semicolon.                                                                                                                    |          |
| -vd       | --vcf_dbvar_input                  | String  |         | dbVar vcf variants file paths delimited by semicolon.                                                                                                                      |          |
| -mi       | --main_input                       | String  |         | Main variant file path used to determine main technology and input between other inputs.                                                                                   |          |
| -vfp      | --vcf_filter_pass                  |         | Off     | Process only structural variants with filter value PASS.                                                                                                                   |          |
| -g        | --gene_intersection                |         | Off     | Overlapping genes filter (i.e. variants with non-overlapping genes are filtered out).                                                                                      |          |
| -svt      | --prefer_base_svtype               |         | Off     | Whether to prefer base variant type (SVTYPE) in case of BND variant and 10x/TELL-Seq (default off i.e. preferring SVTYPE2).                                                |          |
| -d        | --distance_variance                | Integer |         | Distance variance filter - number of bases difference between compared variants.                                                                                           |          |
| -i        | --intersection_variance            | Double  |         | Intersection variance filter - threshold difference between compared variants.                                                                                             |          |
| -mp       | --minimal_proportion               | Double  |         | Minimal proportion filter - minimal proportion of target variant within query variant (0.0 - 1.0).                                                                         |          |
| -t        | --variant_type                     | String  |         | Variant type filter, any combination of [BND,CNV,DEL,INS,DUP,INV,UNK], delimited by semicolon, only variant types listed will processed.                                   |          |
| -rff      | --region_filter_file               | String  |         | List of regions to be excluded from analysis (bed format, tab separated).                                                                                                  |          |
| -gf       | --gene_file                        | String  |         | File containing gene information list (i.e. gene symbol, chromosome, start, end) - gene annotation.                                                                        |          |
| -dvs      | --distance_variance_statistics     | String  |         | Distance variance statistics - bases counts delimited by semicolon (e.g. 10000;50000;100000). Number of variants having distance variance at least 1000, 50000 bases, etc. |          |
| -ivs      | --intersection_variance_statistics | String  |         | Intersection variance statistics - thresholds delimited by semicolon (e.g. 0.1;0.3;0.5). Number of variants having intersection variance score at least 0.1, 0.3, etc.     |          |
| -so       | --statistics_output                | String  |         | Output structural variants statistics csv file path.                                                                                                                       |          |
| -o        | --output                           | String  |         | Output result file path.                                                                                                                                                   | \*       |


## Example of usage
Some basic example usage of structural variant comparator follows. More detailed usage with sample data and results are presented in sample package in <b>./example</b> directory in this repository. Each directory in sample package contains README.txt file where can be found detailed description of each file. Runnable binary version of the application is presented in <b>./bin</b> directory of the repository.

### Basic usage
In basic setup, application compares all SVs contained in input files. No filters are applied here.

```console
java -jar lorec-comparator.jar -vl "longranger_01.vfc;longranger_02.vcf" -vm "manta_01.vcf;manta_02.vcf;manta_03.vcf" -a annotsv.tsv -b bionano.smap -mi manta_01.vcf -o result.csv
```

### Distance Variance filter
Following command filters out variants which have distance variance greater than 50000 bases.

```consolev
java -jar lorec-comparator.jar -vl longranger.vcf -b bionano.smap -d 50000 -o result.csv 
```

### Overlapping genes filter
Following command filters out variants which have distance variance greater than 50000 bases and have no genes in overlap.

```console
java -jar lorec-comparator.jar -vl longranger.vcf -b bionano.smap -d 50000 -g -o result.csv 
```

### Variant type filter
Following command will analyze only translocations (BND), deletions (DEL) and insertions (INS). Other variant types are ignored.

```console
java -jar lorec-comparator.jar -vl longranger.vcf -b bionano.smap -t "BND,DEL,INS" -o result.csv 
```

## Contact
If you have any problem or questions about the software tool, please contact us.

Tomáš Novosád (tomas.novosad@vsb.cz)

Jakub Savara (jakub.savara@vsb.cz)

