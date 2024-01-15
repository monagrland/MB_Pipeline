# Configuration

To configure the `MB_Pipeline` metabarcoding workflow, modify the configuration
file `config.yaml` and input file list `reads.tsv` in this folder with your
desired options and the paths to your input files and reference databases.


## Input file table

Sample names and paths to the input files for the pipeline should be listed in
`reads.tsv` in tab-separated format, with the following fields:

 * `sample` - Short identifier the sample (alphanumeric and underscore only)
 * `lib` - If there is more than one read file/library for a given sample, disambiguate them here.
 * `fwd` - Path to forward read file, relative to Snakemake working directory. Absolute paths are not recommended.
 * `rev` - Path to reverse read file if reads are paired-end; leave blank if reads are single-end


## Config file structure

Pipeline parameters are specified in the `config.yaml` file in this folder.

Paths in the config file are relative to the working directory that Snakemake
is run from; alternatively, absolute paths can be specified but are not
recommended.


### Path to input file table

```yaml
reads_table: config/reads.tsv
```


### Paired- or single-end reads
If paired end reads are used, specify this at the `paired` key. This key
accepts only `true` or `false`.

```yaml
paired: true
```


### Adapter Trimming

Adapter trimming options passed to
[Cutadapt](https://cutadapt.readthedocs.io/en/stable/) are under
`adapter_trimming_options`.

The 5' primer sequence (`5p`), 3' primer sequence (`3p`) and minimum number of
overlapping bases to be trimmed (`min_overlap`) are specified with keywords.

Other Cutadapt options can be specified with their command line flags under
`others_common`, the same way one would usually specify them in a Cutadapt
command.

```yaml
adapter_trimming_options:
  5p: "CCNGAYATRGCNTTTYCCNCG"
  3p: "TANACYTCNGGRTGNCCRAARAAYCA"
  min_overlap: 23
  others_common: "--minimum-length 1"
```


### Merging paired reads

To merge the forward and reverse reads in a read pair, the `--fastq_mergepairs`
argument of the [VSEARCH](https://github.com/torognes/vsearch) tool is used;
other parameters can also be specified and will be passed verbatim to VSEARCH.

```yaml
merge_options:
  - "--fastq_allowmergestagger"
  - "--fastq_minovlen 100"
  - "--fastq_maxdiffs 15"
  - "--fastq_eeout"
```

### Quality Filtering

Quality filtering parameters will also be passed verbatim to VSEARCH.

```yaml
filter_options:
  - "--fastq_maxee 1.0"
  - "--fastq_minlen 200"
  - "--fastq_maxlen 500"
  - "--fastq_maxns 0"
  - "--fasta_width 0"
```


### Dereplication

Exactly identical sequences are removed with VSEARCH `--derep_fulllength`.


### Protein coding sequences

Protein-coding sequences (e.g. the mitochondrial cytochrome oxidase I marker
sequence) can be processed differently from non-coding sequences (e.g. tRNA or
rRNA markers). Instead of filtering for chimeras with UCHIME (see below),
putative pseudogenes are removed that have excessive in-frame stop codons or
where the translation does not match a HMM of the target protein.

Activate the coding sequence-specific subworkflow with:

```yaml
protein_coding: true
```

Genetic code and a HMM file of the target protein should be supplied to screen
translated sequences for pseudogenes; PCR chimeras should also be filtered out
in this step. The approach is adapted from [Porter & Hajibabaei,
2021](https://doi.org/10.1186/s12859-021-04180-x).

```yaml
coding:
  frame: 3 # default for the Leray fragment of mtCOI
  code: 5 # Genetic code, must not be a stopless code
  hmm: null # path to the HMM file, relative to workdir
```

Entropy ratio-based distance denoising with DnoisE (see below) is only
available for coding sequences; the reading frame must be specified. However it
is possible to denoise with DnoisE but still use UCHIME to remove chimeras.


### Denoising

Denoising is performed with Unoise ([Edgar,
2016](https://doi.org/10.1101/081257 )) implemented in
[VSEARCH](https://github.com/torognes/vsearch) by default, but DnoisE (see
below) is an option for coding sequences.

Specify either `unoise` or `dnoise` to the key `method` under `denoising`.

```yaml
denoising:
  method: 'unoise'
  alpha: 5
  minsize: 8
```

Both methods use the parameters `alpha` and `minsize`.

Parameter `alpha` controls the tradeoff between "sensitivity to small
differences against an increase in the number of bad sequences which are
wrongly predicted to be good." Higher values of `alpha` retain more sequences
(more sensitive, more bad sequences), whereas lower values retain fewer (less
sensitive, fewer bad sequences).

`minsize` is the minimum number of sequences represented by a cluster after
denoising.


#### Denoising with DnoisE

This will perform entropy-based distance denoising with
[DnoisE](https://github.com/adriantich/DnoisE/) ([Antich et al.,
2022](https://doi.org/10.7717/peerj.12758)) instead of VSEARCH UNOISE.

The expected reading frame of the amplified metabarcoding fragment should be
known, based on the PCR primers used, and denote the codon position (1, 2, or
3) of the first base in the fragment.

DnoisE can calculate the entropy ratio of codon positions 2 and 3 to help set
values of the denoising parameter alpha and the minimum cluster size. Specify
the range of alpha and minsize values to test:

```yaml
alpha_range: [1,2,3,4,5,6,7,8,9,10]
minsize_range: [2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100]
```

The pipeline first runs with the default alpha and minsize values (see below),
and also produces plots of entropy ratio vs. alpha and minsize for the
specified ranges. After reviewing the plots, the user can update the values for
alpha and minsize and rerun the pipeline if necessary.


### Screening for non-target sequences and artefacts

Non-target sequences (e.g. pseudogenes, non-target amplicons) and technical
artefacts (e.g. PCR chimeras) should be removed from the sequences.

Coding sequences are screened with an HMM of the target protein, and for
excessive in-frame stop codons (see above). This procedure should also
indirectly remove PCR chimeras.

For non-coding sequences are screened for PCR chimeras with UCHIME de-novo
implemented in VSEARCH.


### Community Table Creation

Reads are aligned to the final set of ASVs to calculate abundance per ASV per
sample.

```yaml
community_table_options:
  - "--id 0.97"
  - "--strand plus"
```


### Taxonomic classification


#### Program for phylogenetic inference

The phylogenetic tree can be built with either
[FastTree](http://www.microbesonline.org/fasttree/) or
[IQ-TREE](http://www.iqtree.org/).

```yaml
treeprog: "fasttree" # Options: fasttree, iqtree
```


#### Reference databases for taxonomic classification

The database should contain unaligned reference sequences for the target gene
in Fasta format, along with the taxonomic classification in SINTAX format in
the sequence header, which is described in the VSEARCH manual:

> The reference database must contain taxonomic information in the header of
> each sequence in the form of a string starting with ";tax=" and followed by a
> comma-separated list of up to eight taxo- nomic identifiers. Each taxonomic
> identifier must start with an indication of the rank by one of the letters d
> (for domain) k (kingdom), p (phylum), c (class), o (order), f (family), g
> (genus), or s (species). The letter is followed by a colon (:) and the name
> of that rank. Commas and semicolons are not allowed in the name of the rank.
>
> Example:
> ">X80725_S000004313;tax=d:Bacteria,p:Proteobacteria,c:Gammaproteobacteria,o:Enterobacteriales,f:Enterobacteriaceae,g:Escherichia/Shigella,s:Escherichia_coli".

The database can also be already pre-formatted into UDB format with the
`--makeudb_usearch` option in VSEARCH.  For large databases, preparing a UDB
file is faster, as it avoids re-indexing every time the pipeline is run.

```yaml

direct_dbs: # Path is relative to workdir
  - "testdata/db/bold_coi-5p_test.fasta"

hierarchical_db:
  "testdata/db/bold_coi-5p_test.fasta"
```


#### Classification Thresholds

```yaml
classification_threshold: 0.90

hierarchical_threshold: 0.8
```


