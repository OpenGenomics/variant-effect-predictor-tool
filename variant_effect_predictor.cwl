cwlVersion: v1.0
class: CommandLineTool
doc: |
  For more info see -> http://uswest.ensembl.org/info/docs/tools/vep/script/vep_options.html
  VEP determines the effect of your variants (SNPs, insertions, deletions, 
  CNVs or structural variants) on genes, transcripts, and protein sequence, 
  as well as regulatory regions. Simply input the coordinates of your variants 
  and the nucleotide changes to find out -
    * genes and transcripts affected by the variants
    * location of the variants (e.g. upstream of a transcript, in coding sequence, in non-coding RNA, in regulatory regions)
    * consequence of your variants on the protein sequence (e.g. stop gained, missense, stop lost, frameshift)
    * known variants that match yours, and associated minor allele frequencies from the 1000 Genomes Project
    * SIFT and PolyPhen scores for changes to protein sequence
    * And more!

requirements:
  DockerRequirement:
    dockerPull: "opengenomics/variant-effect-predictor"
        
baseCommand: 
  - "variant_effect_predictor.pl"

arguments:
  - "--no_progress"
  - "--no_stats"
  - "--offline"

inputs:
  # basic options
  everything:
    type: boolean
    doc: |
      Shortcut flag to switch on all of the following -
      --sift b, --polyphen b, --ccds, --uniprot, --hgvs, --symbol, --numbers, 
      --domains, --regulatory, --canonical, --protein, --biotype, --uniprot, 
      --tsl, --appris, --gene_phenotype --af, --af_1kg, --af_esp, --af_exac, 
      --max_af, --pubmed, --variant_class
    inputBinding:
      prefix: "--everything"

  fork:
    type: int?
    doc: |
      Enable forking, using the specified number of forks. Forking can 
      dramatically improve the runtime of the script. Not used by default
    inputBinding:
      prefix: "--fork"

  # input options
  inputFile:
    type: File
    doc: "Input file name"
    inputBinding:
      prefix: "--input_file"

  format:
    type: string?
    doc: |
      Input file format - one of ensembl, vcf, hgvs, id. 
      By default, the script auto-detects the input file format. Using this 
      option you can force the script to read the input file as Ensembl, VCF, 
      IDs or HGVS. Auto-detects format by default
    inputBinding:
      prefix: "--format"

  outputFile:
    type: File
    doc: "Output file name"
    default: "variant_effect_output.txt"
    inputBinding:
      prefix: "--output_file"

  species:
    type: string?
    doc: |
      Species for your data. This can be the latin name e.g. homo_sapiens or 
      any Ensembl alias e.g. mouse. Specifying the latin name can speed up 
      initial database connection as the registry does not have to load all 
      available database aliases on the server. Default = "homo_sapiens"
    inputBinding:
      prefix: "--species"

  assembly:
    type: string?
    doc: |
      Select the assembly version to use if more than one available. If using 
      the cache, you must have the appropriate assemblys cache file installed. 
      If not specified and you have only 1 assembly version installed, this will 
      be chosen by default. Default = use found assembly version
    inputBinding:
      prefix: "--assembly"

  # cache options
  dirCache:
    type: Directory
    doc: "Specify the cache directory to use"
    inputBinding:
      prefix: "--dir_cache"

  dirPlugins:
    type: Directory?
    doc: "Specify the plugin directory to use"
    inputBinding:
      prefix: "--dir_plugins"

  fasta:
    type: 
      - "null"
      - File
      - Directory
    doc: |
      Specify a FASTA file or a directory containing FASTA files to use to look 
      up reference sequence. The first time you run the script with this 
      parameter an index will be built which can take a few minutes. 
      This is required if fetching HGVS annotations (--hgvs) or checking 
      reference sequences (--check_ref) in offline mode (--offline), 
      and optional with some performance increase in cache mode (--cache). 
      See documentation for more details. Not used by default
    inputBinding:
      prefix: "--fasta"

  refseq:
    type: boolean
    doc: |
      Specify this option if you have installed the RefSeq cache in order for 
      VEP to pick up the alternate cache directory. This cache contains transcript 
      objects corresponding to RefSeq transcripts (to include CCDS and Ensembl ESTs 
      also, use --all_refseq). Consequence output will be given relative to these 
      transcripts in place of the default Ensembl transcripts (see documentation)
    inputBinding:
      prefix: "--refseq"

  merged:
    type: boolean
    doc: |
      Use the merged Ensembl and RefSeq cache. Consequences are flagged with the 
      SOURCE of each transcript used.
    inputBinding:
      prefix: "--merged"
  
  cacheVersion:
    type: int?
    doc: |
      Use a different cache version than the assumed default (the VEP version). 
      This should be used with Ensembl Genomes caches since their version numbers 
      do not match Ensembl versions. For example, the VEP/Ensembl version may be 
      88 and the Ensembl Genomes version 35. Not used by default
    inputBinding:
      prefix: "--cache_version"

  bufferSize:
    type: int?
    doc: |
      Sets the internal buffer size, corresponding to the number of variants that 
      are read in to memory simultaneously. Set this lower to use less memory at 
      the expense of longer run time, and higher to use more memory with a faster run time. 
      Default = 5000
    inputBinding:
      prefix: "--buffer_size"
  
  # Other annotation sources
  plugin:
    doc: |
      Use named plugin. Plugin modules should be installed in the Plugins 
      subdirectory of the VEP cache directory (defaults to $HOME/.vep/). 
      Multiple plugins can be used by supplying the --plugin flag multiple times. 
      See plugin documentation. Not used by default
    type:
      - "null"
      - type: array
        items: string
        inputBinding:
          prefix: "--plugin"
  
  gff:
    type: File?
    doc: |
      Use GFF transcript annotations in [filename] as an annotation source. 
      Requires a FASTA file of genomic sequence. Not used by default
    inputBinding:
      prefix: "--gff"

  gtf:
    type: File?
    doc: |
      Use GTF transcript annotations in [filename] as an annotation source. 
      Requires a FASTA file of genomic sequence. Not used by default
    inputBinding:
      prefix: "--gtf"

  bam:
    type: File?
    doc: |
      ADVANCED Use BAM file of sequence alignments to correct transcript models 
      not derived from reference genome sequence. Used to correct RefSeq transcript 
      models. Not used by default
    inputBinding:
      prefix: "--bam"

  # Output options
  variantClass:
    type: boolean
    doc: "Output the Sequence Ontology variant class. Not used by default"
    inputBinding:
      prefix: "--variant_class"

  sift:
    type:
      - "null"
      - type: enum
        symbols: 
          - "p"
          - "s"
          - "b"
    doc: |
      Species limited SIFT predicts whether an amino acid substitution affects 
      protein function based on sequence homology and the physical properties of 
      amino acids. VEP can output the prediction term (p), score (s) or both (b). 
      Not used by default
    inputBinding:
      prefix: "--sift"

  polyphen:
    type:
      - "null"
      - type: enum
        symbols: 
          - "p"
          - "s"
          - "b"
    doc: |
      Human only PolyPhen is a tool which predicts possible impact of an amino acid 
      substitution on the structure and function of a human protein using straightforward 
      physical and comparative considerations. VEP can output the 
      prediction term (p), score (s) or both (b). VEP uses the humVar score by 
      default - use --humdiv to retrieve the humDiv score. Not used by default
    inputBinding:
      prefix: "--polyphen"
    
  nearest:
    type:
      - "null"
      - type: enum
        symbols:
          - "transcript"
          - "gene"
          - "symbol"
    doc: |
      Retrieve the transcript or gene with the nearest protein-coding transcription 
      start site (TSS) to each input variant. Use transcript to retrieve the 
      transcript stable ID, gene to retrieve the gene stable ID, or symbol to 
      retrieve the gene symbol. Note that the nearest TSS may not belong to a 
      transcript that overlaps the input variant, and more than one may be reported 
      in the case where two are equidistant from the input coordinates.
    inputBinding:
      prefix: "--nearest"

  humandiv:
    type: boolean
    doc: |
      Human only Retrieve the humDiv PolyPhen prediction instead of the default 
      humVar. Not used by default
    inputBinding:
      prefix: "--humandiv"

  genePhenotype:
    type: boolean
    doc: |
      Indicates if the overlapped gene is associated with a phenotype, disease 
      or trait. See list of phenotype sources. Not used by default
    inputBinding:
      prefix: "--gene_phenotype"

  regulatory:
    type: boolean
    doc: |
      Look for overlaps with regulatory regions. The script can also call if a 
      variant falls in a high information position within a transcription factor 
      binding site. Output lines have a Feature type of RegulatoryFeature or 
      MotifFeature. Not used by default
    inputBinding:
      prefix: "--regulatory"

  cellType:
    type: string[]?
    doc: |
      Report only regulatory regions that are found in the given cell type(s). 
      Can be a single cell type or a comma-separated list. The functional type in 
      each cell type is reported under CELL_TYPE in the output. Not used by default
    inputBinding:
      prefix: "--cell_type"
      itemSeparator: ","
    
  individual:
    type: string[]?
    doc: |
      Consider only alternate alleles present in the genotypes of the specified 
      individual(s). May be a single individual, a comma-separated list or all
      to assess all individuals separately. Individual variant combinations 
      homozygous for the given reference allele will not be reported. Each 
      individual and variant combination is given on a separate line of output. 
      Only works with VCF files containing individual genotype data; individual 
      IDs are taken from column headers. Not used by default
    inputBinding:
      prefix: "--individual"

  phased:
    type: boolean
    doc: |
      Force VCF genotypes to be interpreted as phased. For use with plugins that 
      depend on phased data. Not used by default
    inputBinding:
      prefix: "--phased"

  alleleNumber:
    type:
      - "null"
      - type: enum
        symbols:
          - "1"
          - "2"
    doc: |
      Identify allele number from VCF input, where 1 = first ALT allele, 
      2 = second ALT allele etc. Useful when using --minimal Not used by default
    inputBinding:
      prefix: "--allele_number"

  totalLength:
    type: boolean
    doc: |
      Give cDNA, CDS and protein positions as Position/Length. Not used by default
    inputBinding:
      prefix: "--total_length"

  numbers:
    type: boolean
    doc: |
      Adds affected exon and intron numbering to to output. 
      Format is Number/Total. Not used by default
    inputBinding:
      prefix: "--numbers"

  domains:
    type: boolean
    doc: |
      Adds names of overlapping protein domains to output. Not used by default
    inputBinding:
      prefix: "--domains"

  noEscape:
    type: boolean
    doc: "Don't URI escape HGVS strings. Default = escape"
    inputBinding:
      prefix: "--no_escape"

  keepCSQ:
    type: boolean
    doc: "Don't overwrite existing CSQ entry in VCF INFO field. Overwrites by default"
    inputBinding:
      prefix: "--keep_csq"

  vcfInfoField:
    type: string?
    doc: |
      Change the name of the INFO key that VEP write the consequences to in its 
      VCF output. Use ANN for compatibility with other tools such as snpEff. 
      Default = CSQ
    inputBinding:
      prefix: "--vcf_info_field"

  terms:
    type:
      - "null"
      - type: enum
        symbols:
          - "ensembl"
          - "so"
    doc: |
      The type of consequence terms to output. The Ensembl terms are described here. 
      The Sequence Ontology is a joint effort by genome annotation centres to 
      standardise descriptions of biological sequences. Default = SO
    inputBinding:
      prefix: "--terms"

  # Identifiers
  hgvs:
    type: boolean
    doc: |
      Add HGVS nomenclature based on Ensembl stable identifiers to the output. Both 
      coding and protein sequence names are added where appropriate. To generate HGVS 
      identifiers when using --cache or --offline you must use a FASTA file and --fasta.
      HGVS notations given on Ensembl identifiers are versioned. Not used by default
    inputBinding:
      prefix: "--hgvs"

  hgvsg:
    type: boolean
    doc: |
      Add genomic HGVS nomenclature based on the input chromosome name. To generate 
      HGVS identifiers when using --cache or --offline you must use a FASTA file and 
      --fasta. Not used by default
    inputBinding:
      prefix: "--hgvsg"

  shiftHgvs:
    type:
      - "null"
      - type: enum
        symbols:
          - "0"
          - "1"
    doc: |
      Enable or disable 3-prime shifting of HGVS notations. When enabled, this causes 
      ambiguous insertions or deletions (typically in repetetive sequence tracts) to be 
      shifted to their most 3-prime possible coordinates (relative to the transcript 
      sequence and strand) before the HGVS notations are calculated; the flag 
      HGVS_OFFSET is set to the number of bases by which the variant has shifted, 
      relative to the input genomic coordinates. Disabling retains the original input 
      coordinates of the variant. Default = 1 (shift)
    inputBinding:
      prefix: "--shift_hgvs"

  protein:
    type: boolean
    doc: "Add the Ensembl protein identifier to the output where appropriate. Not used by default"
    inputBinding:
      prefix: "--protein"

  symbol:
    type: boolean
    doc: "Adds the gene symbol (e.g. HGNC) (where available) to the output. Not used by default"
    inputBinding:
      prefix: "--symbol"

  ccds:
    type: boolean
    doc: "Adds the CCDS transcript identifer (where available) to the output. Not used by default"
    inputBinding:
      prefix: "--ccds"

  uniprot:
    type: boolean
    doc: |
      Adds best match accessions for translated protein products from three 
      UniProt-related databases (SWISSPROT, TREMBL and UniParc) to the output. 
      Not used by default
    inputBinding:
      prefix: "--uniprot"

  tsl:
    type: boolean
    doc: |
      Adds the transcript support level for this transcript to the output. 
      NB - not available for GRCh37. Not used by default
    inputBinding:
      prefix: "--tsl"

  appris:
    type: boolean
    doc: |
      Adds the APPRIS isoform annotation for this transcript to the output. 
      NB - not available for GRCh37. Not used by default
    inputBinding:
      prefix: "--appris"

  canonical:
    type: boolean
    doc: |
      Adds a flag indicating if the transcript is the canonical transcript for 
      the gene. Not used by default
    inputBinding:
      prefix: "--canonical"

  biotype:
    type: boolean
    doc: |
      Adds the biotype of the transcript or regulatory feature. Not used by default
    inputBinding:
      prefix: "--biotype"

  xrefRefseq:
    type: boolean
    doc: |
      Output aligned RefSeq mRNA identifier for transcript. NB - theRefSeq and 
      Ensembl transcripts aligned in this way MAY NOT, AND FREQUENTLY WILL NOT, 
      match exactly in sequence, exon structure and protein product. Not used by default
    inputBinding:
      prefix: "--xref_refseq"

  synonyms:
    type: File?
    doc: |
      Load a file of chromosome synonyms. File should be tab-delimited with the primary 
      identifier in column 1 and the synonym in column 2. Synonyms are used bi-directionally 
      so columns may be switched. Synoyms allow different chromosome identifiers to be 
      used in the input file and any annotation source (cache, database, GFF, custom 
      file, FASTA file). Not used by default
    inputBinding:
      prefix: "--synonyms"

  # Co-located variants
  checkExisting:
    type: boolean
    doc: |
      Checks for the existence of known variants that are co-located with your input. 
      By default the alleles are compared - to compare only coordinates, use --no_check_alleles. 
    inputBinding:
      prefix: "--check_existing"

  excludeNullAlleles:
    type: boolean
    doc: |
      Do not include variants with unknown alleles when checking for co-located variants. 
      The human variation database contains variants from HGMD and COSMIC for which the 
      alleles are not publically available; by default these are included when using 
      --check_existing, use this flag to exclude them. Not used by default
    inputBinding:
      prefix: "--exclude_null_alleles"

  noCheckAlleles:
    type: boolean
    doc: |
      When checking for existing variants, by default VEP only reports a co-located 
      variant if none of the input alleles are novel. For example, if the user input 
      has alleles A/G, and an existing co-located variant has alleles A/C, the 
      co-located variant will not be reported. Strand is also taken into account - 
      in the same example, if the user input has alleles T/G but on the negative 
      strand, then the co-located variant will be reported since its alleles match 
      the reverse complement of user input. Use this flag to disable this behaviour 
      and compare using coordinates alone. Not used by default
    inputBinding:
      prefix: "--no_check_alleles"

  af:
    type: boolean
    doc: |
      Add the global allele frequency (AF) from 1000 Genomes Phase 3 data for any 
      known co-located variant to the output. For this and all --af_* flags, the 
      frequency reported is for the input allele only, not necessarily the 
      non-reference or derived allele. Not used by default
    inputBinding:
      prefix: "--af"

  maxAf:
    type: boolean
    doc: |
      Report the highest allele frequency observed in any population from 1000 
      genomes, ESP or ExAC. Not used by default
    inputBinding:
      prefix: "--max_af"

  af1kg:
    type: boolean
    doc: |
      Add allele frequency from continental populations (AFR,AMR,EAS,EUR,SAS) of 
      1000 Genomes Phase 3 to the output. Must be used with --cache. Not used by default
    inputBinding:
      prefix: "--af_1kg"

  afExac:
    type: boolean
    doc: |
      Include allele frequency from ExAC project populations. Must be used with --cache Not used by default
    inputBinding:
      prefix: "--af_exac"

  afEsp:
    type: boolean
    doc: |
      Include allele frequency from NHLBI-ESP populations. Must be used with --cache Not used by default
    inputBinding:
      prefix: "--af_esp"

  pubmed:
    type: boolean
    doc: |
      Report Pubmed IDs for publications that cite existing variant. Must be used with --cache. Not used by default
    inputBinding:
      prefix: "--pubmed"

  failed:
    type:
      - "null"
      - type: enum
        symbols:
          - "0"
          - "1"
    doc: |
      When checking for co-located variants, by default the script will exclude 
      variants that have been flagged as failed. Set this flag to include such 
      variants. Default = 0 (exclude)
    inputBinding:
      prefix: "--failed"

  # Data format options
  vcf:
    type: boolean
    doc: "Write output in VCF format"
    inputBinding:
      prefix: "--vcf"

  tab:
    type: boolean
    doc: "Write output in tab-delimited format"
    inputBinding:
      prefix: "--tab"

  json:
    type: boolean
    doc: "Write output in JSON format"
    inputBinding:
      prefix: "--json"

  compressOutput:
    type:
      - "null"
      - type: enum
        symbols:
          - "gzip"
          - "bgzip"
    doc: "Writes output compressed using either gzip or bgzip. Not used by default"
    inputBinding:
      prefix: "--compress_output"

  fields:
    type: string[]?
    doc: |
      Configure the output format using a comma separated list of fields. Fields 
      may be those present in the default output columns, or any of those that appear 
      in the Extra column (including those added by plugins or custom annotations). 
      Output remains tab-delimited. Can only be used with tab or VCF format output. 
      Not used by default
    inputBinding:
      prefix: "--fields"
      itemSeparator: ","

  minimal:
    type: boolean
    doc: |
      Convert alleles to their most minimal representation before consequence 
      calculation i.e. sequence that is identical between each pair of reference 
      and alternate alleles is trimmed off from both ends, with coordinates adjusted 
      accordingly. Note this may lead to discrepancies between input coordinates 
      and coordinates reported by VEP relative to transcript sequences; to avoid 
      issues, use --allele_number and/or ensure that your input variants have unique 
      identifiers. The MINIMISED flag is set in the VEP output where relevant. 
      Not used by default
    inputBinding:
      prefix: "--minimal"

  # Filtering and QC options
  gencodeBasic:
    type: boolean
    doc: |
      Limit your analysis to transcripts belonging to the GENCODE basic set. This 
      set has fragmented or problematic transcripts removed. Not used by default
    inputBinding:
      prefix: "--gencode_basic"

  allRefseq:
    type: boolean
    doc: |
      When using the RefSeq or merged cache, include e.g. CCDS and Ensembl EST 
      transcripts in addition to those from RefSeq (see documentation). 
      Only works when using --refseq or --merged
    inputBinding:
      prefix: "--all_refseq"

  excludePredicted:
    type: boolean
    doc: |
      When using the RefSeq or merged cache, exclude predicted transcripts (i.e. 
      those with identifiers beginning with XM_ or XR_).
    inputBinding:
      prefix: "--exclude_predicted"

  checkRef:
    type: boolean
    doc: |
      Force the script to check the supplied reference allele against the sequence 
      stored in the Ensembl Core database or supplied FASTA file. Lines that do 
      not match are skipped. Not used by default
    inputBinding:
      prefix: "--check_ref"

  dontSkip:
    type: boolean
    doc: "Don't skip input variants that fail validation, e.g. those that fall on unrecognised sequences"
    inputBinding:
      prefix: "--dont_skip"

  allowNonVariant:
    type: boolean
    doc: |
      When using VCF format as input and output, by default VEP will skip 
      non-variant lines of input (where the ALT allele is null). Enabling this 
      option the lines will be printed in the VCF output with no consequence data added.
    inputBinding:
      prefix: "--allow_non_variant"

  chr:
    type: string[]?
    doc: |
      Select a subset of chromosomes to analyse from your file. Any data not on 
      this chromosome in the input will be skipped. The list can be comma separated, 
      with - characters representing an interval. For example, to include chromosomes 
      1, 2, 3, 10 and X you could use --chr 1-3,10,X. Not used by default
    inputBinding:
      prefix: "--chr"
      itemSeparator: ","

  codingOnly:
    type: boolean
    doc: "Only return consequences that fall in the coding regions of transcripts. Not used by default"
    inputBinding:
      prefix: "--coding_only"

  noIntergenic:
    type: boolean
    doc: "Do not include intergenic consequences in the output. Not used by default"
    inputBinding:
      prefix: "--no_intergenic"

  pick:
    type: boolean
    doc: |
      Pick once line or block of consequence data per variant, including transcript-specific 
      columns. Consequences are chosen according to the criteria described here, and the 
      order the criteria are applied may be customised with --pick_order. This is the 
      best method to use if you are interested only in one consequence per variant. 
      Not used by default
    inputBinding:
      prefix: "--pick"

  pickAllele:
    type: boolean
    doc: |
      Like --pick, but chooses one line or block of consequence data per variant allele. 
      Will only differ in behaviour from --pick when the input variant has multiple 
      alternate alleles. Not used by default
    inputBinding:
      prefix: "--pick_allele"

  perGene:
    type: boolean
    doc: |
      Output only the most severe consequence per gene. The transcript selected is 
      arbitrary if more than one has the same predicted consequence. Uses the same 
      ranking system as --pick. Not used by default
    inputBinding:
      prefix: "--per_gene"

  pickAlleleGene:
    type: boolean
    doc: |
      Like --pick_allele, but chooses one line or block of consequence data per 
      variant allele and gene combination. Not used by default
    inputBinding:
      prefix: "--pick_allele_gene"

  flagPick:
    type: boolean
    doc: |
      As per --pick, but adds the PICK flag to the chosen block of consequence data 
      and retains others. Not used by default
    inputBinding:
      prefix: "--flag_pick"

  flagPickAllele:
    type: boolean
    doc: |
      As per --pick_allele, but adds the PICK flag to the chosen block of consequence 
      data and retains others. Not used by default
    inputBinding:
      prefix: "--flag_pick_allele"

  flagPickAlleleGene:
    type: boolean
    doc: |
      As per --pick_allele_gene, but adds the PICK flag to the chosen block of 
      consequence data and retains others. Not used by default
    inputBinding:
      prefix: "--flag_pick_allele_gene"

  pickOrder:
    type: string[]?
    doc: |
      Customise the order of criteria applied when choosing a block of annotation 
      data with e.g. --pick. Valid criteria are - canonical,appris,tsl,biotype,ccds,rank,length
    inputBinding:
      prefix: "--pick_order"
      itemSeparator: ","

  mostSevere:
    type: boolean
    doc: |
      Output only the most severe consequence per variant. Transcript-specific
      columns will be left blank.  Not used by default
    inputBinding:
      prefix: "--most_severe"

  summary:
    type: boolean
    doc: |
      Output only a comma-separated list of all observed consequences per variant. 
      Transcript-specific columns will be left blank. Not used by default
    inputBinding:
      prefix: "--summary"

  filterCommon:
    type: boolean
    doc: |
      Shortcut flag for the filters below - this will exclude variants that have a 
      co-located existing variant with global AF > 0.01 (1%). May be modified using 
      any of the following freq_* filters. Not used by default
    inputBinding:
      prefix: "--filter_common"

  checkFrequency:
    type: boolean
    doc: |
      Turns on frequency filtering. Use this to include or exclude variants based 
      on the frequency of co-located existing variants in the Ensembl Variation database. 
      You must also specify all of the --freq_* flags below. Frequencies used in 
      filtering are added to the output under the FREQS key in the Extra field. 
      Not used by default
    inputBinding:
      prefix: "--check_frequency"

  freqPop:
    type:
      - "null"
      - type: enum
        symbols:
          - "1KG_ALL"
          - "1KG_AFR"
          - "1KG_AMR"
          - "1KG_EAS"
          - "1KG_EUR"
          - "1KG_SAS"
          - "ESP_AA"
          - "ESP_EA"
          - "ExAC"
          - "ExAC_Adj"
          - "ExAC_AFR"
          - "ExAC_AMR"
          - "ExAC_EAS"
          - "ExAC_FIN"
          - "ExAC_NFE"
          - "ExAC_SAS"
          - "ExAC_OTH"
    doc: "Name of the population to use in frequency filter."
    inputBinding:
      prefix: "--freq_pop"
  
  freqFreq:
    type: float?
    doc: "Allele frequency to use for filtering. Must be a float value between 0 and 1"
    inputBinding:
      prefix: "--freq_freq"

  freqGtLt:
    type: 
      - "null"
      - type: enum
        symbols:
          - "gt"
          - "lt"
    doc: |
      Specify whether the frequency of the co-located variant must be greater 
      than (gt) or less than (lt) the value specified with --freq_freq
    inputBinding:
      prefix: "--freq_gt_lt"

  freqFilter:
    type: 
      - "null"
      - type: enum
        symbols:
          - "exclude"
          - "include"
    doc: "Specify whether to exclude or include only variants that pass the frequency filter"
    inputBinding:
      prefix: "--freq_filter"

outputs:
  annotatedVCF:
    type: File
    outputBinding:
      glob: $(inputs.outputFile)
