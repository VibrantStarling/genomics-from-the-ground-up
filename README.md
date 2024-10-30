# Eukaryotic Gene Prediction and Annotation
Automated gene prediction and annotation is an imperfect process. It can be messy, and it is not always better than manual annotation. Its result still need to be checked manually to make sure they make sense. Perhaps one day there will be a trust worthy quality metric, but today is not that day.

I always recommend loading your genome into a genome viewer like Artemis, Apollo, or IGV **and** checking summary stats make sense for your organism (e.g. average gene length or number of gene).


## Software you'll need
The singularity image for `dfam-TEtool` (dfam-tetools-latest.sif) in your home directory. This allows access to `Repeatmasker` and `RepeatModeller`

`hisat2`, the aligner we will be using

`trimmomatic`, to clean up messy RNAseq data

`samtools`, to handle SAM and BAM files (the raw aligned data)

`seqkit`, to quickly and cleanly handle and edit fastq files




```python
singularity pull dfam-tetools-latest.sif docker://dfam/tetools:latest
conda install -c bioconda hisat2 trimmomatic samtools seqkit
```

# The steps to annotate a genome *in silico*

## 1. Prepare your genome and RNAseq data

### **Retrieve your genome fasta** 

Download your genome from NCBI or ENA etc. with `wget` or assemble it yourself from your own data.

### **Download your SRA files**

You can mass download SRA data using Michael Gerth's [SRA download script](https://github.com/gerthmicha/perlscripts/blob/master/sra_download.pl)

Create a text file list of SRA identifiers with one identifier per line like this:

```
SRR000001
SRR000002
SRR000003
SRR000004
```

Then pass it to the script and check out all of your download stats with `seqkit stats` to make sure everything has downloaded correctly. If any pairs of RNA do not match in the number of reads (num_seqs), re-download them.

### **Repeat mask your genome** 

Annotating and masking your repeats is a vital step in gene annotation. It cuts down on spurious annotations and highlights important drivers of evolution in your genome such as transposable elements.

Here we use TETools `RepeatModeler2` to identify repeats *de novo*, followed by `RepeatMasker` to softmask repeats in the genome.

**Soft masking** makes repeats appear as lowercase letters in the DNA sequence and means you do not lose the genetic information in those zones.
It looks like this:


```
ATGCCGCAAAAAAATTTTTAGGC --> ATGCCGCaaaaaaatttttAGGC
```

Hard masking replaces repeat regions with Ns and means you lose information. Avoid doing this unless you absolutely need to. It looks like this:
```
ATGCCGCAAAAAAATTTTTAGGC --> ATGCCGCNNNNNNNNNNNNAGGC
```

*NOTE*

Maskers can struggle if memory is not fast enough and data may need to be stored on compute HPC nodes


```python
DB=output_name_for_repeat_db
GENOME=genome.fa

tput setaf 5; echo "softmasking genome"; tput sgr0
singularity exec -B ${PWD}:${PWD},${HOME} ${HOME}/dfam-tetools-latest.sif BuildDatabase -name ${DB} ${GENOME}
singularity exec -B ${PWD}:${PWD},${HOME} ${HOME}/dfam-tetools-latest.sif RepeatModeler -database ${DB} -threads 32 -LTRStruct
singularity exec -B ${PWD}:${PWD},${HOME} ${HOME}/dfam-tetools-latest.sif RepeatMasker -pa 32 -lib ${DB}-families.fa -xsmall ${GENOME}
tput setaf 2; echo "softmasking complete"; tput sgr0
```

## 2a. Preparing your files for, and running BRAKER3

[Here's a video](https://www.youtube.com/watch?v=UXTkJ4mUkyg) on how to run all versions of BRAKER and GALBA (feat. OMArk and BUSCO). I highly reccomend taking the time to watch it.

#### **Clean your RNAseq data**

You don't know where your RNAseq data has come from. Quality check (QC) it and give it a clean.

`fastQC` can be used to QC and will rpoduce some nice visuals for you. A tutorial video can be found [here](https://www.youtube.com/watch?v=bz93ReOv87Y).

Next we need to clean it. We will be using `Trimmomatic`. Trimmomatic can be a little tempermental but it does its job. If you want to know more, read [this](https://www.biocomputix.com/post/trimming-ngs-data-trimmomatic#viewer-6fi0o)

You'll need to pass trimmomatic the location of an adaptor file. These can be found on their Github. In the code below replace `TruSeq3-PE.fa` with the path/name.fa of your adpaters (probably also TruSeq3-PE.fa).



```python
# define your SRA list file
SRA_LIST=sra-list.txt

# for PAIRED data
for RNA_PREFIX in $(cat ${SRA_LIST})
do
    RNASEQ_FWD=${RNA_PREFIX}_1_clean.fastq.gz
    RNASEQ_REV=${RNA_PREFIX}_2_clean.fastq.gz
    
    # trim the fastq file with trimmomatic
    tput setaf 6; echo "------START of trimming for ${RNA_PREFIX}"; tput sgr0
    trimmomatic PE -phred33 -threads 32 ${RNASEQ_FWD} ${RNASEQ_REV} ${RNA_PREFIX}_fpaired.fq.gz ${RNA_PREFIX}_funpaired.fq.gz ${RNA_PREFIX}_rpaired.fq.gz ${RNA_PREFIX}_runpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
    echo "Number of trimmed forward paired reads: " 
    echo $(zcat ${RNA_PREFIX}_rpaired.fq.gz |wc -l)/4|bc
    echo "Number of trimmed reverse paired reads: " 
    echo $(zcat ${RNA_PREFIX}_fpaired.fq.gz |wc -l)/4|bc
    tput setaf 2; echo "------END of  trimming for ${RNA_PREFIX}------"; tput sgr0
    echo
done
```


```python
# for UNPAIRED data
for RNA_PREFIX in $(cat ${SRA_LIST})
do
    RNASEQ=${RNA_PREFIX}_clean.fastq.gz
    trimmomatic SE -phred33 -threads 32 ${RNASEQ} ${RNA_PREFIX}_trimmed.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
done
```

#### **Align your RNAseq data**

Next, use `Hisat2` to align your RNAseq data to your genome. 

Another aligner you could use for this is STAR, but for most purposes Hisat2 is good enough. STAR is more computationally expensive and runs more slowly, but can be better with draft genome and poor quality genomes. If you're interested in knowing more, start [here first](https://www.biostars.org/p/288726/), then [here](https://pmc.ncbi.nlm.nih.gov/articles/PMC5792058/), and then [here](https://pmc.ncbi.nlm.nih.gov/articles/PMC7084517/).


```python
# name your file prefixes
NAME=name-of-your-alignment-run
IDX=prefix-for-your-database

# build the database for hisat2. This only needs to be done once per genome.
hisat2-build ${GENOME} ${IDX}
```


```python
# gather your forward and reverse reads from trimmomatic into comma delimited lists with no white spaces
FWD_FILES=$(ls -m SRR*fpaired.fq.gz)
FWD_FILES=$(sed -s 's/ //g' $FWD_FILES)
REV_FILES=$(ls -m SRR*rpaired.fq.gz | sed -s 's/ //g')
REV_FILES=$(sed -s 's/ //g' $REV_FILES)

# run hisat2 for PAIRED data
hisat2 -p 32 -q -x ${IDX} -1 ${FWD_FILES} -2 ${REV_FILES} > ${NAME}-hisat2-rnaseq.sam  2> ${NAME}-hisat2-align.err
```


```python
# run hisat2 for PAIRED data
hisat2 -p 32 -q -x ${IDX} -U ${RNASEQ} > ${NAME}-hisat2-rnaseq.sam  2> ${NAME}-hisat2-align.err
```


```python
# turn your SAM files into sorted BAM files
samtools view -bS -@ 12 ${NAME}-hisat2-rnaseq.sam -o ${NAME}-hisat2-rnaseq.bam
samtools sort -@ 12 ${NAME}-hisat2-rnaseq.bam -o ${NAME}-hisat2-rnaseq_sorted.bam
rm ${NAME}-hisat2-rnaseq.sam; rm ${NAME}-hisat2-rnaseq.bam
```

### **Run BRAKER3**

[BRAKER3](https://github.com/Gaius-Augustus/BRAKER) requires three main inputs to run 1) your genome 2) your aligned reads for that genome 3) a reference protein database for your organism. You can make your own protein database, but the BRAKER3 team have complied some [prepartitioned collections](https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/) from OrthoDB for you. 

**We want to run BRAKER3 with several additional flags:**

`--threads`, define the number of threads to use

`--gff`, outputs results as a gff

`--workingdir`, names the output directory


```python
T=32
SORTED_BAM="rnaseq_sorted.bam"

# This combines the input bam name and potein database used into one name
# This will be your output directory
WD=$(basename -s .bam ${SORTED_BAM})_$(basename -s .fa ${PROT_DB})

# run BRAKER3
singularity exec -B ${PWD}:${PWD},${HOME} ${HOME}/braker3.sif braker.pl --genome=${GENOME} --prot_seq=${HOME}/BRAKER-DB/Alveolata.fa --bam=${SORTED_BAM} --threads=${T} --gff
```

#### **Getting BUSCO scores**

You may also want to run BUSCO to check your outputs. You can do this in BRAKER via `compleasm` by adding `--busco_lineage=lineage` and naming a [BUSCO lineage](https://busco.ezlab.org/list_of_lineages.html) like 'alveolata_odb10' or 'fungi_odb10'.

Alternatively, run BUSCO through [compleasm](https://github.com/huangnengCSU/compleasm) yourself.

## 2b. Run Funnannotate

Perhaps your BRAKER3 run failed. Perhaps you just want to run funannotate. Here is a how to do that.

Funannotate is very talkative, and provides many more supporting files so you can really go to town figuring out what has happened during the gene prediction process. It can also predict UTR regions.
However. For most genomes, Funannotate simply does not perform as well as BRAKER3. One expection to this seems to be Fungi (shocking that it performs best on the organisms it was designed for, huh?)

Anyway, running funannotate is arguably simpler than BRAKER and you can add a whole host of supporting data to help it annotate gene models, such as proten databases, other gene predictor results, and bam files. The software will also point you in the right direction for functional annotation with antiSMASH and interproscan (don't run their commands for it if you are on HPC clusters though, at least one uses docker which is not and will not be install for security reasons).

All you need to give funannotate is:

- the `softmasked genome` 
- all the cleaned (but *NOT* trimmed) SRA fastq files. Funannotate will run trimmomatic itself unless you tell it not to.

You may want to edit `--max_intronlen` to be appropriate for your species.

You may also want to turn off `--repeats2evm` if you are dealing with a fungi genome with high gene density. This option is best for large genomes and genomes with high repeat content

If you want to learn more about funannotate's gene prediction options, read [this](https://funannotate.readthedocs.io/en/stable/predict.html#gene-prediction).


```python
GENOME=softmasked-genome.fasta
OUT_DIR=funannotate-output-name
SPECIES="SpeciesNameWithNoSpaces"
STRAIN="StrainName"

#TRAIN
singularity exec -B ${PWD}:${PWD},${HOME} ${HOME}/funannotate.sif funannotate train -i ${GENOME} -o ${OUT_DIR} --left [ES]RR*_1.fastq.gz --right [ES]RR*_2.fastq.gz  --species ${SPECIES} --strain ${STRAIN} --cpus 48 --max_intronlen 100000
#PREDICT
singularity exec -B ${PWD}:${PWD},${HOME} ${HOME}/funannotate.sif funannotate predict -i ${GENOME} -o ${OUT_DIR} -s ${SPECIES} --strain ${STRAIN} --cpus 48 --repeats2evm --organism other
#UPDATE
singularity exec -B ${PWD}:${PWD},${HOME} ${HOME}/funannotate.sif funannotate update -i ${OUT_DIR} -g ${GENOME} --cpus 48 --max_intronlen 100000
```

## 3. What next? 


### Predict UTRs
BRAKER3 does not have a stable way to predict UTR regions, but Funannotate does. Luckily, you can use `funannotate update` on its own to predict UTR regions through PASA. See the documentation [here](https://funannotate.readthedocs.io/en/latest/update.html).

There are many other programmes designed to predict 3' and 5' UTRs for different organisms, and like gene models, UTR prediction may not always be biologically correct. However, UTRs are important for evolution and gene function, so it is worth giving a go and assessing for yoursel if the UTRs predicted make sense.

remember to set your `--max_intronlen` to something sensible for your species


```python
# predict UTRs with funannotate update
singularity exec -B ${PWD}:${PWD},${HOME} ${HOME}/funannotate.sif funannotate update --fasta genome.fa --gff braker.gff -o output-dir-name --species "Species name" --cpus 48 --max_intronlen 100000 --left [ES]RR*_1.fastq.gz --right [ES]RR*_2.fastq.gz 
```

### Functional annotation
So far all you have is *structural annotation*. If you want to predict what a gene *might* do, you need to do *functional annotation*. There are several things you can do:
- GO term annotation with [EggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) or [interproscan](https://interproscan-docs.readthedocs.io/en/latest/HowToRun.html)
- interpro domain, pfam, PANTHER, and Gene3D annotation through [interproscan](https://interproscan-docs.readthedocs.io/en/latest/HowToRun.html)
- secondary metabolite prediction with [antiSMASH](https://docs.antismash.secondarymetabolites.org/intro/) for fungi, plants, and microbes.

Following functional annotation, you may want to do a [pathway enrichment analysis](https://geneontology.org/docs/go-enrichment-analysis/) to assess what genes are most prevelant in you sample.

I highly recommend running interproscan with the `--goterms` and `--iprlookup` flags, plus `-appl` with any additional analyses you're interested in.

### Compare gff files

Do you have multiple annotations for one genome, including a reference which you would like to compare? Try using [GffCompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml).

Make a textfile list of all the gffs you want to compare. They must all be on the same genome for it to work.

The output `.tracking` file will have a serioes of columns, one for each gff. They will not be named, but they are in order of your input.


```python
gffcompare -r reference-annotation.gff -i input-gffs.txt -o gffcompare-output-prefix
```

You may visualise your `.tracking` file as stacked bar charts with this python code:


```python
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def plot_stats(file):      
    palette = ['#9e9e9e',
                 '#ffa200',
                 '#502F4C',
                 '#70587c',
                 '#C8B8DB',
                 '#56EEF4',
                 '#74B3CE',
                 '#034E21',
                 '#DDCC77',
                 '#AA4499',
                 '#FF5376',
                 '#882255',
                 '#09BC8A',
                 '#005C69']
    sns.set(font_scale=1, rc={"figure.figsize":(10,7)})
    sns.set_palette(palette=palette)
    df3 = df.iloc[:,4:-1]!='-'
    df3['Class code'] = df['Class code']
    df3 = df3.sort_values(by='Class code')
    fig, ax = plt.subplots(figsize=(8, 6))
    df3.groupby('Class code').sum().T.plot(kind='bar', stacked=True, ax=ax)
    plt.title(outname)
    plt.legend(loc='upper left',bbox_to_anchor=(1, 1))
    return df


file = 'gffcompare.tracking'
df = pd.read_csv(file, sep='\t', header=None)
outname = Path(file).stem
df.rename(columns={0:'Query transfrag id',
                        1:'Query locus id',
                        2:'Reference gene id',
                        3:'Class code',
                        4:'gff 1 name',
                        5:'gff 2 name',
                        6:'gff 3 name'
                        }, inplace=True)
code_dict = {'=':'Exact match of intron chain', 
                'c':'Contained in reference (intron compatible)', 
                'm':'All introns matched or retained', 
                'o':'Other same strand overlap with reference exons', 
                'j':'Multi-exon with atleast one junction match', 
                'i':'Fully contained within a reference intron', 
                'n':'Not all introns matched or retained', 
                'e':'Single exon transfrag partially covering an intron', 
                'x':'Exonic overlap on opposite strand', 
                'p':'Possible polymerase run-on (no actual overlap)', 
                'u':'Unknown, intergenic', 
                'k':'Containment of reference (reverse containment)'}
df['Class name'] = df['Class code'].map(code_dict)

plot_stats(df)
```
