# Genome assembly 
Protocol I follow to assemble the reference genomes of _Calliphora hilli, Calliphora quadrimaculata, Calliphora stygia,_ and __Calliphora vicina_. The data is ONT and was sequenced at Bragata on a PromethION platform. Three libraries were made, resulting in three datasets (R0406, R0408, and R0410). Basecalling was done individually on each library. I received files in POD5 format. The first step is to basecall using a high accuracy model. 

# 1: basecalling using Dorado
I used DORADO (https://github.com/nanoporetech/dorado). You first need to create a sample file (in .csv format) that contains metadata for your samples like so:

**run R0406**

| Kit            | Experiment ID | Flow Cell ID | Barcode   | Alias      |
|----------------|---------------|--------------|-----------|------------|
| SQK-NBD114-96  | R0406         | PAY72012     | barcode63 | 3_stygia   |
| SQK-NBD114-96  | R0406         | PAY72012     | barcode64 | 4_hilli    |

**run R0408**

| Kit            | Experiment ID | Flow Cell ID | Barcode   | Alias               |
|----------------|---------------|--------------|-----------|----------------------|
| SQK-NBD114-96  | R0408         | PAY72012     | barcode61 | 1_quadrimaculata     |
| SQK-NBD114-96  | R0408         | PAY72012     | barcode62 | 2_vicina             |
| SQK-NBD114-96  | R0408         | PAY72012     | barcode63 | 3_stygia             |
| SQK-NBD114-96  | R0408         | PAY72012     | barcode64 | 4_hilli              |

**run R0410**

| Kit            | Experiment ID | Flow Cell ID | Barcode   | Alias               |
|----------------|---------------|--------------|-----------|----------------------|
| SQK-NBD114-96  | R0410         | PAY78695     | barcode61 | 1_quadrimaculata     |
| SQK-NBD114-96  | R0410         | PAY78695     | barcode62 | 2_vicina             |
| SQK-NBD114-96  | R0410         | PAY78695     | barcode63 | 3_stygia             |
| SQK-NBD114-96  | R0410         | PAY78695     | barcode64 | 4_hilli              |


NOTE: you will need the kit name and experiment ID to match what is in your data. If you didn't recieve this information from your sequencer, it is possible to retrieve this information from the POD5 files using the POD5 tool (https://github.com/nanoporetech/pod5-file-format).

The CSV file needs to be in the same folder as the POD5 files. If you have two libraries to basecall, then you need to have two folders and run each library seperately. Here is an example code for basecalling (do for every library):

```
module load Dorado/0.9.1

dorado basecaller sup /nesi/nobackup/uow03920/BLOWFLY_ASSEMBLY_DATA/02_basecalling/run_0410/pod5 --recursive --device 'cuda:all' --kit-name SQK-NBD114-96 --sample-sheet 0410.csv > /nesi/nobackup/uow03920/BLOWFLY_ASSEMBLY_DATA/03_BAM/0410_sup_calls_rsumd.bam
```

# 2: demultiplexing using DEMUX
I used the DORADO demux tool (https://github.com/nanoporetech/dorado), running all three libraries seperately. The output is four BAM files (i.e. one BAM file per sample).

```
module load Dorado/0.9.1
dorado demux --output-dir nesi/nobackup/uow03920/05_blowfly_assembly_march/04_demux/0406 --no-classify 0406_sup_calls_rsumd.bam
```

# 3: convert bam files to FASTQ files
We used SAMtools because it is the best for insect genomes. Other tools are highly optimised for humans etc. SAMtools is also the fastest. I have provided an example of one library (0408) here, but you will need to do it for every library (e.g., 0406, 0410).

```
module load SAMtools

samtools fastq 48305005d92ac294badde036140df6e8ce951184_1_quadrimaculata.bam > /nesi/nobackup/uow03920/05_blowfly_assembly_march/05_fastq/0408_quadrimaculata.fastq
samtools fastq 48305005d92ac294badde036140df6e8ce951184_2_vicina.bam > /nesi/nobackup/uow03920/05_blowfly_assembly_march/05_fastq/0408_vicina.fastq
samtools fastq 48305005d92ac294badde036140df6e8ce951184_3_stygia.bam > /nesi/nobackup/uow03920/05_blowfly_assembly_march/05_fastq/0408_stygia.fastq
samtools fastq 48305005d92ac294badde036140df6e8ce951184_4_hilli.bam > /nesi/nobackup/uow03920/05_blowfly_assembly_march/05_fastq/0408_hilli.fastq
```

# 4: merge the libraries together so that we are only working off one file
I have provided one example here, but do it for every species. 

```
cat 0406/0406_hilli.fastq  0408/0408_hilli.fastq 0410/0410_hilli.fastq > 01_hilli.fastq
```

# 5: QC prior to filtering using Nanoplot

```
module load NanoPlot

for i in 01_hilli 02_quadrimaculata 03_stygia 04_vicina; do
  NanoPlot --verbose -t 8 --fastq /nesi/nobackup/uow03920/05_blowfly_assembly_march/06_concatenate/${i}.fastq -o /nesi/nobackup/uow03920/05_blowfly_assembly_march/07_pre_assembly_QC/01_NanoPlot/${i}
done
```

# 6: filtering using CHOPPER
Based off the NanoPlot results, we decided to filter reads that were less than 1000 bases long and with quality scores lower than 10.

```
module load chopper

for i in 01_hilli 02_quadrimaculata 03_stygia 04_vicina; do
chopper --threads 8 -q 10 -l 1000 < /nesi/nobackup/uow03920/05_blowfly_assembly_march/06_concatenate/${i}.fastq > /nesi/nobackup/uow03920/05_blowfly_assembly_march/08_filtered/${i}.fastq ;
done
```

# 7: Rerun Nanoplot to check that the quality improved post-filtering

# 8 : assemble genome using Flye

used a nextflow script that Meeran created:

```
#!/usr/bin/env nextflow

// Define the list of IDs
params.ids = [ '01_hilli', '02_quadrimaculata', '03_stygia', '04_vicina']

// Create a channel from the list of IDs
Channel.from(params.ids)
    .set { id_ch }

// Define the process for running flye
process runFlye {
    cpus 24
    memory '60 GB'
    publishDir("/nesi/nobackup/uow03920/05_blowfly_assembly_march/09_assembly", mode: 'move')
    // Define the input
    input:
    val id from id_ch

    // Define the output
    output:
    file("${id}_flye") into flye_out_ch

    // Define the command to be executed
    script:
    """
    flye --nano-hq /nesi/nobackup/uow03920/05_blowfly_assembly_march/08_filtered/${id}.fastq \
         --out-dir ${id}_flye --threads ${task.cpus} -g 700m
    """
}

// Define what to do with the output (optional)
flye_out_ch.view()
```

and ran it using a slurm script: 

```
ml Flye/2.9.3-gimkl-2022a-Python-3.11.3 
ml Nextflow/22.10.7

nextflow run flye.nf
```


# 9: Initial post-assembly QC using QUAST

```
ml QUAST

quast.py -t 16 -o /nesi/nobackup/uow03920/05_blowfly_assembly_march/10_post_assembly_QC -l '01_hilli 02_quadrimaculata 03_stygia 04_vicina'  /nesi/nobackup/uow03920/05_blowfly_assembly_march/09_assembly/01_hilli/01_hilli.fasta /nesi/nobackup/uow03920/05_blowfly_assembly_march/09_assembly/02_quadrimaculata/02_quadrimaculata.fasta /nesi/nobackup/uow03920/05_blowfly_assembly_march/09_assembly/03_stygia/03_stygia.fasta /nesi/nobackup/uow03920/05_blowfly_assembly_march/09_assembly/04_vicina/04_vicina.fasta
```

# 10: Initial post-assembly completeness QC using compleasm (BUSCO)

```
module load compleasm/0.2.5-gimkl-2022a

for i in 01_hilli 02_quadrimaculata 03_stygia 04_vicina;
do
compleasm.py run -a /nesi/nobackup/uow03920/05_blowfly_assembly_march/09_assembly/${i}_flye/${i}.fasta -o /nesi/nobackup/uow03920/05_blowfly_assembly_march/10_post_assembly_QC/02_busco/${i}/${i}_compleasm -l diptera_odb10 -t 8 ;
done
```

# 11: purge haplotigs 
Compelasm revealed that there were quite a few duplicated genes, so we will purge haplotigs to try and remove them.

First align the filtered .fastq files to the .fasta sequences to generate .paf files. Make a folder for each sample and run the code in those folders. 

```
module load minimap2

for i in 01_hilli 02_quadrimaculata 03_stygia 04_vicina; do
minimap2 -c /nesi/nobackup/uow03920/05_blowfly_assembly_march/09_assembly_${i}_flye/${i}.fasta /nesi/nobackup/uow03920/05_blowfly_assembly_march/08_filtered/${i}.fastq > /nesi/nobackup/uow03920/05_blowfly_assembly_march/11_purge_haplotigs/${i}/${i}_alignment.paf
done
```

Second, generate statistics for each .paf file: 

```
ml purge_dups
```
```
pbcstat 01_hilli_alignment.paf
```

Third, calculate the cutoffs:

```
calcuts PB.stat > cutoffs 2> calcuts_01_hilli.log
```

Fourth, split the consensus fasta file:

```
split_fa /nesi/nobackup/uow03920/05_blowfly_assembly_march/09_assembly/01_hilli_flye/01_hilli.fasta > con_split_01_hilli.fa
```

Fifth, generate a self mapping .paf file:

```
minimap2 -xasm5 -DP con_split_01_hilli.fa con_split_01_hilli.fa | gzip -c - > con_split_01_hilli.self.paf.gz
```

purge duplicates:

```
purge_dups -2 -T cutoffs -c PB.base.cov con_split_01_hilli.self.paf.gz > dups_01_hilli.bed 2> purge_dups_01_hilli.log
```

extract sequences:

```
get_seqs -e dups_01_hilli.bed 01_hilli.fasta
```

# 12: Repeat QC (Compleasm & QUAST) to check that contiguity and completeness improved.

# 13: polishing using MEDAKA

```
for i in 01_hilli 02_quadrimaculata 03_stygia 04_vicina; do
    medaka_consensus -t $SLURM_CPUS_PER_TASK -i /nesi/nobackup/acc_name/PX024_Parasitoid_wasp/02_Basecaller_sup/03_fastqfiles/02_filtered/MO_${i}_cat_fil.fastq -d MO_${i}_flye/assembly.fasta -o MO_${i}_flye_medaka ;
done
```

# 14: scaffold using LongStitch
I used LongStitch for de novo scaffolding in all species.

```
longstitch run draft=draft-assembly reads=reads G=gsize
``` 

I also did additional scaffolding for _C. vicina_ (because there is a high quality reference genome available; Sivell, 2024) using RagTag

```
    ragtag.py scaffold 04_vicina_reference.fa \
    04_vicina_longstitch.fa \
    -o 04_vicina_ragtag.fa
```

# 15: QC of all assemblies again using QUAST and BUSCO to check improvements

# 16: BlobToolKit for QC and detecting contamination

BlobToolKit is used here for quality control and to detect potential contamination in the genome assembly. Several files need to be generated for BlobToolKit to function effectively:

-Blast output (to identify taxonomic affiliations)
-Aligned Nanopore reads (BAM file with .csi index)
-Assembly (FASTA format)
-BUSCO report (CSV format)

## BLAST the assembly against the NT database

```
module load BLASTDB/2024-07 BLAST/2.13.0-GCC-11.3.0

# Run BLAST for each scaffolded genome
    blastn \
    -query 04_vicina.ragtag.scaffold.fasta \
    -task megablast \
    -db nt \
    -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
    -culling_limit 10 \
    -num_threads 42 \
    -evalue 1e-3 \
    -out 04_vicina_megablast.out;
```

## blob to check for contamination

Blob uses the older version of BUSCO, so need to first run that first:

```
ml miniprot
ml busco

busco -i /nesi/nobackup/uow03920/05_blowfly_assembly_march/19_scaffold/01_hilli/01_hilli_scaffold.fasta -l diptera_odb10 -o 01_hilli -m genome --out_path ./01_hilli
```

Then run blobtools. You need to download the taxdump file from NCBI and put it in the directory with a link to it. 

```
blobtools create --busco /nesi/nobackup/uow03920/05_blowfly_assembly_march/23_busco_for_blob/01_hilli/01_hilli/run_diptera_odb10/full_table.tsv --fasta /nesi/nobackup/uow03920/05_blowfly_assembly_march/19_scaffold/01_hilli/01_hilli_scaffold.fasta --cov /nesi/nobackup/uow03920/05_blowfly_assembly_march/20_contamination/01_minimap/01_hilli.bam --hits /nesi/nobackup/uow03920/05_blowfly_assembly_march/20_contamination/01_hilli_megablast.out --replace --taxdump /nesi/nobackup/uow03920/05_blowfly_assembly_march/22_blob/new_taxdump 01_hilli

```

To view this: 

download the blobtoolkit-api and blobtoolkit-viewer from https://github.com/genomehubs/blobtoolkit/releases. Put them in your working directory. Make them executable (```chmod +x blobtoolkit-api blobtoolkit-viewer```). Use nesi's ondemand virtual desktop. In the virtual desktop, open a terminal and cd to the working directory (with the api and viewer and data needed; `cd /nesi/nobackup/{id}/05_blowfly_assembly_march/22_blob`). Then add the api to your path: `export=$PATH:$(pwd)`. Then create the viewer `blobtools host --port 8081 --api-port 9001 --hostname localhost /01_hilli`. Then open the link it sends you in the browser within virtual environment. 


download the taxonomy table to use to filter out contaminants!!!

This script removes contaminants and filters out contigs that are less than 1000bp in length (script provided at https://github.com/meeranhussain/Genome_assembly_AND_annotation). I altered the script because i had long read data and Meeran's script used qualimap (i used seqkit instead to filter out reads less than 1000bp). 

This is to remove low quality, small and contaminated contigs.

```
#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: $0 -i <fasta_file> -b <busco_file> -f <blob_file> -o <output_name>"
    echo
    echo "Options:"
    echo "  -i <fasta_file>      Path to the FASTA file"
    echo "  -b <busco_file>      Path to the BUSCO full_table.tsv file"
    echo "  -f <blob_file>       Path to the BlobToolKit CSV file"
    echo "  -o <output_name>     Desired name for the filtered output FASTA"
    echo "  -h                   Display this help message"
    exit 1
}

# Check if seqkit is installed
if ! command -v seqkit &> /dev/null; then
    echo "Error: seqkit is not installed. Please install seqkit to proceed."
    exit 1
fi

# Check if no arguments were passed and display usage
if [ $# -eq 0 ]; then
    usage
fi

# Read command-line arguments
while getopts "i:b:f:o:h" flag; do
    case "${flag}" in
        i) fasta=${OPTARG};;
        b) busco_file=${OPTARG};;
        f) blob_file=${OPTARG};;
        o) output_name=${OPTARG};;
        h) usage;;
        *) usage;;
    esac
done

# Check if all required arguments are provided
if [ -z "${fasta}" ] || [ -z "${busco_file}" ] || [ -z "${blob_file}" ] || [ -z "${output_name}" ]; then
    echo "Error: Missing required arguments."
    usage
fi

############## Logging braces open #####################
{
###### Filtering assembly based on contigs size

contigs_siz_fil() {
    local fasta_file=$1

    # Extract contig names with lengths < 1000bp
    seqkit fx2tab -n -l "$fasta_file" | awk -F'\t' '$2 < 1000 {print $1}'
}

######### Filtering assembly based on contaminated contigs

contigs_cont_fil() {
    local fileA=$1
    # Extract contaminated contigs
    awk -F "," 'NR > 1 {
        gsub(/"/, "", $6);
        gsub(/"/, "", $7);
        if ($6 != "Arthropoda" && $6 != "no-hit") {
            print $7;
        }
    }' "$fileA"
}

echo -e "File $fasta is loaded\n"
echo -e "File $busco_file is loaded\n"
echo -e "File $blob_file is loaded\n"

# Call the size filter function
siz_fil_contigs=$(contigs_siz_fil "$fasta")

# Call the contamination filter function
cont_fil_contigs=$(contigs_cont_fil "$blob_file")

# Combine both sets of contigs to remove
echo -e "$siz_fil_contigs\n$cont_fil_contigs" > contigs_rmvd.txt

# SeqKit to filter scaffold
seqkit grep -v -f contigs_rmvd.txt "$fasta" > "$output_name"

echo -e "Filtered fasta is stored in $output_name\n"

# Report contig counts
ori_contig_num=$(grep -c "^>" "$fasta")
fil_contig_num=$(grep -c "^>" "$output_name")

echo -e "Number of contigs\nBefore filtering: $ori_contig_num\nAfter filtering : $fil_contig_num\n"

} 2>&1 | tee -a file.log
```

Run by using a slurm script: 

```
ml SeqKit

for i in 01_hilli 02_quadrimaculata 03_stygia 04_vicina; do 
./filter.sh -i /nesi/nobackup/uow03920/05_blowfly_assembly_march/19_scaffold/${i}/${i}_scaffold.fasta -b /nesi/nobackup/uow03920/05_blowfly_assembly_march/23_busco_for_blob/${i}/${i}/run_diptera_odb10/full_table.tsv -f ${i}.csv -o ${i}_filtered_contaminants;
done
```

# 16: QC again (QUAST and BUSCO)

# ANNOTATION

First we need to softmask using EDTA

## rename the contigs with a sequential naming scheme that adheres to EDTA's requirements
First - you need to rename the contig IDs because EDTA requires the input FASTA file to have contig IDs that are no longer than 13 characters. 

use the following script (make a .sh file)

```
#!/bin/bash

# contig_name.sh: Rename contigs for strict EDTA compatibility (≤13 chars, no underscores)
# Usage: ./contig_name.sh input.fasta output.fasta

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 input.fasta output.fasta"
  exit 1
fi

input="$1"
output="$2"
mapfile="${output}.contig_map.tsv"

awk -v map="$mapfile" 'BEGIN {
  i = 0;
  OFS = "\t";
  print "Old_ID", "New_ID" > map
}
/^>/ {
  old = substr($0, 2);
  new = sprintf("contig%05d", ++i);
  print old, new >> map;
  print ">" new;
  next;
}
{ print }' "$input" > "$output"
```

make executable  `chmod +x contig_name.sh`

and run on each assembly by doing 

` ./contig_name.sh ../24_filtered_contamination/01_hilli_filtered_contamination.fasta 01_hilli_final.fasta ```

This script will also produce a mapping file in case we need to see what the original names of the contigs were. 

## Identify the repeat regions in the genomes
with anno --1, we are soft masking, so changing repeat regions to lowercase letters. 

```
ml EDTA/2.1.0

# Main loop
for i in 01_hilli 02_quadrimaculata 03_stygia 04_vicina; do
  echo "Starting EDTA for $i"

  mkdir -p "$i/EDTA"
  cp "${i}_final.fasta" "$i/EDTA/"
  cd "$i/EDTA"

  EDTA.pl --genome "${i}_final.fasta" --threads 32 --sensitive 1 --anno 1

  cd ../../  # Return to starting directory
  echo "Finished EDTA for $i"
done

```

# Repeat masker
This code masks the repeats that we identified from the above code in a .fasta genome assembly. This will avoid the prediction of false positive gene structures in repetitive and low complexity regions.

```
ml RepeatMasker/4.1.0-gimkl-2020a

# Loop through selected strains and run RepeatMasker
for i in 01_hilli 03_stygia; do
    cd ${i}
    RepeatMasker -pa 32 -xsmall -lib ${i}_final.fasta.mod.EDTA.TElib.fa ${i}_filtered_contaminants.fasta
    cd ../
done
```

# Gene prediction using BRAKER
We will annotate genomes on the repeat-masked genome using homology-based evidence from protein databases, leveraging BRAKER3 pipelines.

Protein evidence used for annotation was a concatenated database comprising:
Uniprot-SwissProt proteins
Busco proteins

## combine the two proteins into one database
```
cat Arthropoda.fasta Uniprot_sprot.fasta > proteins.fasta
```

## Run BRAKER: 
```
ml Apptainer
ml AUGUSTUS

cd /nesi/nobackup/uow03920/05_blowfly_assembly_march/28_annotation/01_hilli

#link in the genome assembly with repetitive regions masked and also the proteins fasta we generated from ortho and busco
ln -s /nesi/nobackup/uow03920/05_blowfly_assembly_march/27_repeat_masker/01_hilli/01_hilli_filtered_contaminants.fasta.masked ./01_hilli_filtered_contaminants.fasta.masked 
ln -s /nesi/nobackup/uow03920/05_blowfly_assembly_march/28_annotation/proteins.fasta ./proteins.fasta

#Run BRAKER3
apptainer exec \
  --bind /nesi/nobackup/uow03920/05_blowfly_assembly_march/28_annotation/01_hilli/augustus_config:/opt/Augustus/config \
  /nesi/nobackup/uow03920/genemark/braker3.sif braker.pl \
    --threads=12 \
    --genome=01_hilli_filtered_contaminants.fasta.masked \
    --prot_seq=proteins.fasta \
    --species=01_hilli \
    --gff3 \
    --AUGUSTUS_ab_initio \
    --crf

cd ../
```

# Check the quality of the annotation using busco with the protein function 

```
ml BUSCO
module load compleasm/0.2.5-gimkl-2022a

LINEAGE=/nesi/nobackup/uow03920/05_blowfly_assembly_march/29_annotation_qc/diptera_odb10

# Run BUSCO on BRAKER3-predicted protein sequences for each strain
compleasm.py protein -p /nesi/nobackup/uow03920/05_blowfly_assembly_march/28_annotation/01_hilli/braker/braker.aa -o /nesi/nobackup/uow03920/05_blowfly_assembly_march/29_annotation_qc/01_hilli -l $LINEAGE -t 12
```

my proteome assemblies had quite a bit of duplications so I ran a bit of extra code to keep only the longest isoforms, and then repeat busco to check that it made duplication go down. 

```
#!/bin/bash -e
#SBATCH --account=uow03920
#SBATCH --job-name=busco_clean_proteins
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paige.matheson14@gmail.com
#SBATCH --output busco_clean_%j.out
#SBATCH --error busco_clean_%j.err

module purge

# Load required tools (adjust if different module names exist on NeSI)
module load AGAT
module load gffread
module load compleasm/0.2.5-gimkl-2022a

# ---- Input files ----
GENOME=/nesi/nobackup/uow03920/05_blowfly_assembly_march/27_repeat_masker/01_hilli/01_hilli_filtered_contaminants.fasta.masked
GFF=/nesi/nobackup/uow03920/05_blowfly_assembly_march/28_annotation/01_hilli/braker/braker.gff3
LINEAGE=/nesi/nobackup/uow03920/05_blowfly_assembly_march/29_annotation_qc/diptera_odb10
OUTDIR=/nesi/nobackup/uow03920/05_blowfly_assembly_march/30_clean_duplicate_proteins/01_hilli

mkdir -p $OUTDIR

# ---- Step 1: Keep only the longest isoform per gene ----
agat_sp_keep_longest_isoform.pl \
  --gff $GFF \
  -o $OUTDIR/annotation_longest.gff3

# ---- Step 2: Extract protein sequences from cleaned GFF ----
gffread $OUTDIR/annotation_longest.gff3 \
  -g $GENOME \
  -y $OUTDIR/annotation_longest.aa

# ---- Step 3: Run BUSCO (protein mode) on the cleaned proteome ----
compleasm.py protein \
  -p $OUTDIR/annotation_longest.aa \
  -o $OUTDIR/busco_protein_longest \
  -l $LINEAGE \
  -t 8
```

# functional annotation using egg nog

```
module load eggnog-mapper/2.1.12-gimkl-2022a

# Loop through ecotypes for functional annotation
for i in 01_hilli 02_quadrimaculata 03_stygia 04_vicina; do
  mkdir ${i}
  cd ${i}

  # Optional: Create softlink to protein FASTA file
 ln -s /nesi/nobackup/uow03920/05_blowfly_assembly_march/30_clean_duplicate_proteins/${i}/annotation_longest.aa ${i}_protein.fa

  # Copy GFF3 file for functional annotation decoration
  cp /nesi/nobackup/uow03920/05_blowfly_assembly_march/30_clean_duplicate_proteins/${i}/annotation_longest.gff3 ./${i}_braker.gff3

  # Run EggNOG-mapper
  emapper.py \
    -i ${i}_protein.fa \
    -o ${i}_eggnog \
    --data_dir /nesi/nobackup/uow03920/05_blowfly_assembly_march/34_functional/eggnog_db \
    --decorate_gff ${i}_braker.gff3 \
    --cpu 22 \
```

# interproscan
```
module load InterProScan/5.66-98.0-gimkl-2022a 
ml Perl
ml Python 
ml PCRE2/10.42-GCCcore-12.3.0 
ml GCC

for i in 05_dros_mel 06_cuprina 07_chryso 08_musca; do

mkdir ${i}
cd ${i}

# Run InterProScan
interproscan.sh \
  -i ../../30_clean_duplicate_proteins/${i}/annotation_longest.aa \
  -t p \
  -dp \
  -appl Pfam,SMART,PrositeProfiles,Gene3D,SUPERFAMILY,PANTHER \
  --goterms \
  --pathways \
  --iprlookup \
  -cpu 22 \
  -T temp/

cd ../
done
```

# COMPARATIVE TRANSPOSABLE ELEMENT ANALYSIS
This section details how to run the transposable element/repeat element analyses. Also note this website: https://cmdcolin.github.io/awesome-genome-visualization/?latest=true&tag=Comparative that has a bunch of genome visualisation packages.

# 1: Identify TEs
I identified transposable elements from the un-masked genome assemblies using EarlGrey - https://github.com/TobyBaril/EarlGrey. Earl grey generates longer, non-redundant TE consensus sequences than older software using de novo and homology-based approaches. It feeds the input genome into RepeatMasker and compares it against the Dfam database to identify well-characterised repeat families. It also runs RepeatModeler2 on any remaining unmasked genomic regions to find novel elements unique to the organism and builds initial bases consensus sequences for newly discovered families. It uses CD-HIT to cluster and collapse duplicates, and merges overlapping definitions into a single, high-fidelity, non-redundant TE library. It finally performans a final RepeatMasker pass over to resolve any lingering overlapping annotations and has a lot of 'paper-ready' tools such as TE landscape plots and GC%. 

I ran this code for each of the six blowfly species: 

```
source /home/pm101/miniconda3/etc/profile.d/conda.sh
conda activate earlgrey
for i in 01_hilli; do
cd ${i}
earlGrey -g ${i}.fa -s 01_hilli -o ${i}_earlgrey -r Diptera -t $SLURM_CPUS_PER_TASK;
done
```

To further improve the quality of the TE library and annotation, I ran TEtrimmer on the EarlGrey outputs. TEtrimmer relies on a core module called TEstrainer to build its library. TEtrimmer fixes divergent family over-simplification, cleans multiple sequence alignments using a sliding window strategy (EarlGrey can sometimes leave behind messy, poorly conseved regions in its final consensus sequences), and pinpoints exact terminal coordinates of the TE (which can sometimes be fuzzy following Earl Grey's elongation loop). 

Run this for each species: 

```
TEtrimmer \
  --input_file 01_hilli_combined_library.fasta \
  --genome_file ../../37_TE_PROPER/02_earlgrey/01_hilli/01_hilli.fa \
  --output_dir 01_hilli_TEtrimmer_out \
  --num_threads 16 \
  --classify_all
```

And then, run RepeatMasker to get an updated TE summary (make sure to run with the -a flag if you are wanting to calculate kimura divergence, as you need the aligned file to do this). 

```
module purge 
module load RepeatMasker 

RepeatMasker \
  -pa 16 \
  -a \
  -gff \
  -lib 01_hilli_TEtrimmer_out/TEtrimmer_consensus_merged.fasta \
  ../../37_TE_PROPER/02_earlgrey/01_hilli/01_hilli.fa
  ```

The above code will produce .tbl files that contain summary information on the % of genome comprised of various transposable element families (which I used to create Fig. 2A). 

# 2: Calculate TE length
To calculate the TE length, I used the following R code (this is an example for one species, need to repeat for each of them):

```
library(dplyr)
library(ggplot2)
library(readr)

# skip the 3 header lines
rm <- read_table(
  "02_quadrimaculata.fa.out",
  skip = 3,
  col_names = FALSE
)

colnames(rm) <- c(
  "SW_score",
  "perc_div",
  "perc_del",
  "perc_ins",
  "chr",
  "start",
  "end",
  "left_query",
  "strand",
  "repeat_name",
  "class_family",
  "repeat_begin",
  "repeat_end",
  "repeat_left",
  "ID",
  "flag"
)

# inspect
glimpse(rm)

rm <- rm %>%
  mutate(length = end - start + 1)

sort(table(rm$class_family), decreasing = TRUE)[1:50]

rm <- rm %>%
  mutate(
    TE_order = case_when(
      grepl("^LTR", class_family, ignore.case = TRUE) ~ "LTR",
      grepl("^LINE", class_family, ignore.case = TRUE) ~ "LINE",
      grepl("^SINE", class_family, ignore.case = TRUE) ~ "SINE",
      grepl("^DNA", class_family, ignore.case = TRUE) ~ "DNA",
      grepl("Helitron|RC", class_family, ignore.case = TRUE) ~ "Rolling_circle",
      grepl("Unknown", class_family, ignore.case = TRUE) ~ "Unclassified",
      TRUE ~ "Other"
    )
  )

table(rm$TE_order)

rm %>%
  group_by(TE_order) %>%
  summarise(
    n = n(),
    mean_length = mean(length),
    median_length = median(length)
  ) %>%
  arrange(desc(mean_length))


#plot all species together on a figure
library(tidyverse)

# Create dataframe
df <- data.frame(
  Species = c("HIL", "QUA", "STY",
              "VIC", "SER", "MEG"),
  LTR = c(280, 653, 227, 235, 187, 437),
  LINE = c(200, 170, 169, 202, 179, 199),
  DNA = c(175, 248, 184, 189, 122, 191),
  Other = c(40, 41, 43, 44, 52, 47),
  Rolling_circle = c(234, 236, 214, 145, 202, 84),
  Unclassified = c(163, 132, 142, 152, 91, 215)
)

# Convert to long format
df_long <- df %>%
  pivot_longer(
    cols = -Species,
    names_to = "TE_Order",
    values_to = "Median_Length"
  )

# Set species order
df_long$Species <- factor(
  df_long$Species,
  levels = c("QUA", "STY", "HIL", "VIC", "SER", "MEG")
)

# Set order of TE categories
df_long$TE_Order <- factor(
  df_long$TE_Order,
  levels = c("LTR", "LINE", "DNA", "Other", "Rolling_circle", "Unclassified")
)


# Plot
length <- ggplot(df_long,
       aes(x = Species,
           y = Median_Length,
           fill = TE_Order)) +
  
  geom_col(position = position_dodge(width = 0.8),
           width = 0.7) + 
  
  scale_fill_manual(values = c(
    "LTR" = "#9385ED",
    "LINE" = "#85EDC7",
    "DNA" = "#85DFED",
    "Other" = "#fab27f",
    "Rolling_circle" = "#E0ED85",
    "Unclassified" = "#ED85DF"
  )) +
  
  theme_classic() +
  
  labs(
    x = element_blank(),
    y = "Median RE Length (bp)",
    fill = "TE Order"
  ) +
  
  theme(
    axis.text.x = element_text(
      size = 18, colour = "black"
    ),
    axis.text.y = element_text(size = 18, colour = "black"),
    axis.title = element_text(size = 18),
    legend.position = 'right',
    axis.line = element_line(linewidth = 0.5)
  )
length

ggsave("length.png", plot = length, dpi = 400)
``` 

This script calculates both mean and median TE length. I chose to present median, since mean was likely influenced by a few elements with extremely long/short length. 

## 3: TE-associated BUSCOs
I got the idea for this analysis from https://pmc.ncbi.nlm.nih.gov/articles/PMC12590894/. The idea is that you can use BUSCO genes as a proxy for detecting genome-wide TE-gene associations. More scripts/information can be found at https://pubmed.ncbi.nlm.nih.gov/35217860/ and https://pubmed.ncbi.nlm.nih.gov/37739812/.  

First, I ran BUSCO on the un-masked genomes: 
```
module load BUSCO/5.8.2-gimkl-2022a
ml miniprot

for i in 01_hilli 02_quadrimaculata 03_stygia 04_vicina 05_sericata 06_megacephala; do 

cd ${i}

busco \
  -i ${i}.fa \
  -o ${i}_busco \
  -l diptera_odb10 \
  -m genome \
  --long \
  -c 24
cd ../
done
```

I then ran the following bash script. This script BLASTs the single copy, complete, BUSCO sequences against it's own genome to see how often they occur.  

```
ml BLAST
ml BEDTools

for i in 02_quadrimaculata; do 

cd ${i}

#grab coordinates of all complete BUSCOs
awk -F'\t' 'NR>1 && $2=="Complete" {
    id=$1
    if (!(id in count)) {
        chr[id]=$3
        s=$4
        e=$5
        if (s=="" || e=="") next
        start[id] = (s < e) ? s : e
        end[id] = (s < e) ? e : s
    }
    count[id]++
}
END {
    for (id in count)
        if (count[id]==1) {
            out = start[id] - 1
            if (out < 0) out = 0
            print chr[id], out, end[id], id
        }
}' OFS='\t' ${i}_full_table.tsv > ${i}_all_single_copy.bed


#extract these sequences and turn into a fasta
bedtools getfasta \
-fi ${i}_masked.fa \
-bed ${i}_all_single_copy.bed \
-name \
-fo ${i}_all_single_copy.fa


#create the BLAST database 
makeblastdb \
  -in ${i}_masked.fa \
  -dbtype nucl \
  -parse_seqids \
  -blastdb_version 5 \
  -out ${i}_softmasked_db


#blast the BUSCO sequences against the genome to see how often they occur
blastn \
-query ${i}_all_single_copy.fa \
-db ${i}_softmasked_db \
-outfmt 6 \
-max_target_seqs 50000 \
-num_threads 32 \
-out ${i}_vs_softmasked.out


cd ../

done

```

Then, need to parse the BLAST output into a table that shows the copy number of each busco gene within an assembly. Sort this from lowest to highest. Calculate the mean number of BLAST hits for the first three quartiles of the ordered data so that any buscos with unexpectedly high hits (i.e., those in the fourth quartile) are removed from the average calculation -- this is the baseline for expected BLAST hits. as per Sproul, we will use any busco with 10X the expected BLAST hits as putative te-associated buscos. I then used the coordinates of the TEs from gff file and the coordinates of BUSCOs to see if they intersect with TE annotations.

```
#!/usr/bin/env Rscript

# List of species directories/file paths
file_paths <- c(
  "01_hilli/01_hilli_vs_softmasked.out",
  "02_quadrimaculata/02_quadrimaculata_vs_softmasked.out",
  "03_stygia/03_stygia_vs_softmasked.out",
  "04_vicina/04_vicina_vs_softmasked.out",
  "05_sericata/05_sericata_vs_softmasked.out",
  "06_megacephala/06_megacephala_vs_softmasked.out"
)

analyze_file <- function(file_path) 
  
  #Read raw BLAST lines
  blast_lines <- readLines(file_path)
  if (length(blast_lines) == 0)

  #Extract full query IDs (Column 1)
  qseqids <- sub("\t.*", "", blast_lines)
  
  #Extract uniquely represented coordinates per query header
  unique_queries <- unique(qseqids)
  
  busco_ids_all <- sub("::.*", "", qseqids)
  
  #Parse coordinates mapping table from unique query headers
  coord_df <- data.frame(
    query = unique_queries,
    busco_id = sub("::.*", "", unique_queries),
    chrom    = sub(".*::([^:]+):.*", "\\1", unique_queries),
    chromStart = as.numeric(sub(".*:([0-9]+)-[0-9]+$", "\\1", unique_queries)),
    chromEnd   = as.numeric(sub(".*-[0-9]+$", "", sub(".*:([0-9]+-[0-9]+)$", "\\1", unique_queries))),
    stringsAsFactors = FALSE
  )
  #Fix chromEnd extraction directly
  coord_df$chromEnd <- as.numeric(sub(".*-", "", sub(".*:", "", unique_queries)))
  
  #Deduplicate to 1 coordinate entry per BUSCO ID
  coord_df <- coord_df[!duplicated(coord_df$busco_id), c("busco_id", "chrom", "chromStart", "chromEnd")]
  
  #Count total copy numbers (BLAST hits)
  counts_table <- table(busco_ids_all)
  busco_df <- data.frame(
    busco_id = names(counts_table),
    copy_number = as.numeric(counts_table),
    stringsAsFactors = FALSE
  )
  
  #Sort ascending, calculate Q3 baseline, and filter >= 10x candidates
  busco_df <- busco_df[order(busco_df$copy_number), ]
  rownames(busco_df) <- NULL
  
  q3_cutoff <- quantile(busco_df$copy_number, 0.75)
  baseline_hits <- mean(busco_df[busco_df$copy_number <= q3_cutoff, ]$copy_number)
  cutoff_10x <- 10 * baseline_hits
  
  te_candidates <- busco_df[busco_df$copy_number >= cutoff_10x, ]
  
  if (nrow(te_candidates) > 0) {
    te_candidates$fold_increase <- round(te_candidates$copy_number / baseline_hits, 2)
    te_candidates <- te_candidates[order(-te_candidates$copy_number), ]
    rownames(te_candidates) <- NULL
  }
  
  cat("Total BUSCOs        :", nrow(busco_df), "\n")
  cat("Baseline Mean (Q1-3):", round(baseline_hits, 2), "hits\n")
  cat("10x Threshold       :", round(cutoff_10x, 2), "hits\n")
  cat("Candidate TE BUSCOs :", nrow(te_candidates), "\n")
  
  #File output
  dir_name <- dirname(file_path)
  species_prefix <- sub("_vs_softmasked\\.out$", "", basename(file_path))
  
  out_all_counts <- file.path(dir_name, paste0(species_prefix, "_copy_numbers_sorted.tsv"))
  out_candidates <- file.path(dir_name, paste0(species_prefix, "_candidate_TE_associated_buscos.tsv"))
  out_bed        <- file.path(dir_name, paste0(species_prefix, "_candidate_TE_buscos.bed"))
  
  # Write TSVs
  write.table(busco_df, file = out_all_counts, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(te_candidates, file = out_candidates, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # 7. Merge candidates with parsed genomic coordinates and write BED
  if (nrow(te_candidates) > 0) {
    bed_df <- merge(te_candidates, coord_df, by = "busco_id")
    
    # 6-column BED format: chrom, chromStart, chromEnd, name, score, strand
    bed_df <- data.frame(
      chrom      = bed_df$chrom,
      chromStart = bed_df$chromStart,
      chromEnd   = bed_df$chromEnd,
      name       = bed_df$busco_id,
      score      = bed_df$copy_number,
      strand     = "."
    )
    
    write.table(bed_df, file = out_bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    cat("BED file generated  :", out_bed, "\n")
  }
}

invisible(lapply(file_paths, analyze_file))
```

convert the Repeatmasker output to gff (this script is in the Repeatmasker repository github) 

```
rmOutToGFF3.pl 01_hilli.fa.out > 01_hilli_repeats.gff3
```

find intersect
```
module load BEDTools

bedtools intersect \
  -a 01_hilli/01_hilli_candidate_TE_buscos.bed \
  -b 01_hilli/01_hilli_repeats.gff3 \
  -wa -wb \
  > 01_hilli/01_hilli_busco_te_overlaps.tsv

cut -f 4 06_megacephala_busco_te_overlaps.tsv | sort -u | wc -l
```

# 4: Kimura divergence/TE landscape plots
This analysis uses the .fa.align files produced by including the -a flag in RepeatMasker earlier. 

First, I used the script `CalcDivergenceFromAlign.pl` script from RepeatMasker (https://github.com/Dfam-consortium/RepeatMasker/blob/master/util/calcDivergenceFromAlign.pl). I put the script in the same working directory as my .fa.align files and ran:

```perl calcDivergenceFromAlign.pl -s 01_hilli.divsum 01_hilli.fa.align```

This generates a .divsum file which can be used in createRepeatLandscape.pl, also implemented in RepeatMaskers github page (https://github.com/Dfam-consortium/RepeatMasker/blob/master/util/createRepeatLandscape.pl). You need the genome size for this. I edited this perl script so that it would generate csv files rather than html files, as I wanted to customise the plots (I've added this customised script in this repository -- it's called 'createRepeatLandscape_generate_csv.pl').

I used this R code to generate the kimura plots for each species:

```
library(tidyverse)

df <- read_csv("01_hilli_landscape.csv")

#make long
long_df <- df %>%
  pivot_longer(
    -Divergence,
    names_to = "Class",
    values_to = "GenomePercent"
  )

#i don't want to plot every single TE family, group by higher level classifications
grouped_df <- long_df %>%
  mutate(
    Group = case_when(
      str_detect(Class, "^LINE") ~ "LINEs",
      str_detect(Class, "^LTR") ~ "LTRs",
      str_detect(Class, "^DNA") ~ "DNA transposons",
      Class == "RC/Helitron" ~ "Rolling circle",
      Class == "Unknown" | Class == "Other" ~ "Unclassified",
      TRUE ~ "Unclassified"
    )
  ) %>%
  group_by(Divergence, Group) %>%
  summarise(
    GenomePercent = sum(GenomePercent),
    .groups = "drop"
  )

#group by levels
grouped_df$Group <- factor(
  grouped_df$Group,
  levels = c(
    "DNA transposons",
    "LINEs",
    "LTRs",
    "Rolling circle",
    "Unclassified"
  )
)

#colours 
fill_colours <- c(
  "DNA transposons" = "#e69f00",
  "LINEs" = "#009e73",
  "LTRs" = "#0072b2",
  "Rolling circle" = "#f0e442",
  "Unclassified" = "#cc79a7"
)


ggplot(grouped_df, aes(x = Divergence, y = GenomePercent, fill = Group)) +
  geom_col(width = 1, color = "black", linewidth = 0.05) +
  scale_x_continuous(
    limits = c(-0.5, 50.5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    expand = c(0, 0)
  ) +
  labs(
    x = "Kimura substitution level (CpG adjusted)",
    y = "Percent of genome",
    fill = NULL
  ) +
  scale_fill_manual(values = fill_colours) +
  theme_minimal() +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 18, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.3),
    axis.ticks = element_line(colour = "black", linewidth = 0.3),
    axis.ticks.length = unit(0.2, "cm"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )
```

# 5: Expanding/contracting gene families
Using the protein fastas we generated earlier using BRAKER. 

## Run OrthoFinder
I provided OrthoFinder with my own Newick Tree because orthofinder wasn't producing an accurate species tree. I believe this is because the blowflies are too closely related, so it gets a bit messy. 

Orthofinder is used to infer orthogroups (genes descended from a single gene in the last common ancestor), orthologs, rooted gene trees, and rooted species trees. It is essential for identifying gene duplication events, analyse gene family evolution (what I used it for), and compare protein sets. 

My Newick tree was called tree.nwk and was in the same folder as OrthoFinder:
```(06_megacephala_protein:30,(05_sericata_protein:19.75186453,(((01_hilli_protein:5.066546484,03_stygia_protein:5.066546484)1:4.01645853,02_quadrimaculata_protein:9.083005013)1:4.536153609,04_vicina_protein:13.61915862)1:6.132705906)1:10.24813547);
```

Run OrthoFinder. 'proteins' is a folder containing all of the protein files generated in previous steps (e.g., 01_hilli_longest_isoforms_modified.faa)

```
ml OrthoFinder/2.5.2 MAFFT/7.505 IQ-TREE/2.2.2.2

# Run:
orthofinder -f proteins -M msa -T iqtree -t 8 -a 4
```

I mainly used OrthoFinder to generate files that could be used in CAFE3

## Gene family expansion and contraction using CAFE3
Convert our Newick Tree to Ultramentric tree (chronos method) in R. You need to estimate total root to tip time. 

```
library(ape)
#Load the rooted species tree
tree <- read.tree("SpeciesTree_rooted.txt")

#Check that it's rooted and binary
stopifnot(is.rooted(tree))
stopifnot(is.binary(tree))

#Convert to ultrametric with a relaxed clock
ultra_tree <- chronos(tree, lambda = 1, model = "correlated")

#Scale tree so root-to-tip distance = 30 Mya
tree_age <- max(node.depth.edgelength(ultra_tree))
ultra_tree$edge.length <- ultra_tree$edge.length * (22 / tree_age)

#Plot and save
plot(ultra_tree, main = "Ultrametric Tree (Root Scaled to 22 Mya)")
write.tree(ultra_tree, file = "SpeciesTree_ultrametric_22Mya.txt")  
``` 

### Filter the orthogroup gene count file (from OrthoFinder) and modify it for CAFE3
cafetutorial_clade_and_size_filter.py is available here: https://github.com/hahnlab/CAFE5/blob/master/docs/tutorial/clade_and_size_filter.py. 

```
awk 'OFS="\t" {$NF=""; print}' Orthogroups.GeneCount.tsv > tmp \
&& awk '{print "(null)""\t"$0}' tmp > cafe.input.tsv \
&& sed -i '1s/(null)/Desc/g' cafe.input.tsv \
&& rm tmp

python2 cafetutorial_clade_and_size_filter.py \
    -i cafe.input.tsv -o filtered.cafe.input.tsv -s 2
```
### Run CAFE3 

```
ml GCC/12.3.0 
ml OpenBLAS

# Set variables
CAFE5_EXEC="/nesi/nobackup/uow03920/06_blowfly_assembly_jan/32_cafe/CAFE5/build/cafe5"
CAFE_DIR="/nesi/nobackup/uow03920/06_blowfly_assembly_jan/38_cafeagain/"
ORTHO_DIR="/nesi/nobackup/uow03920/06_blowfly_assembly_jan/36_orthofinder/proteins/OrthoFinder/Results_May07"
TREE="/nesi/nobackup/uow03920/06_blowfly_assembly_jan/38_cafeagain/SpeciesTree_ultrametric_17Mya.txt"
COUNTS="/nesi/nobackup/uow03920/06_blowfly_assembly_jan/38_cafeagain/cafe_input_filtered.tsv"
K_VALS=(1 2 3 4 5 6)
REPLICATES=(1 2)
PVALUE_THRESHOLD=0.05

# Step 1: Run CAFE5 with multiple k values and replicates
for k in "${K_VALS[@]}"; do
    for rep in "${REPLICATES[@]}"; do
        OUTDIR="$CAFE_DIR/k${k}_rep${rep}"
        mkdir -p "$OUTDIR"
        echo "Running CAFE5 for k=$k, replicate=$rep..."
        $CAFE5_EXEC -i "$COUNTS" -t "$TREE" -o "$OUTDIR" -k $k > "$OUTDIR/stdout_log.txt"
    done
done

# Step 2: Extract log-likelihoods, check convergence, find best k
best_k=""
best_lnl=999999

echo -e "\n--- Log-likelihoods and convergence check ---"
for k in "${K_VALS[@]}"; do
    lnls=()
    for rep in "${REPLICATES[@]}"; do
        OUTDIR="$CAFE_DIR/k${k}_rep${rep}"
        LNL=$(grep "Final -lnL:" "$OUTDIR/stdout_log.txt" | sed -E 's/.*Final -lnL: //')
        echo "k=$k, rep=$rep: lnL=$LNL"
        lnls+=("$LNL")
    done

    diff=$(echo "${lnls[0]} - ${lnls[1]}" | bc -l | tr -d '-')
    echo "k=$k: ΔlnL between replicates = $diff"

    # Use lnL from first replicate for best model selection
    cmp=$(echo "${lnls[0]} < $best_lnl" | bc)
    if [[ "$cmp" -eq 1 ]]; then
        best_k=$k
        best_lnl=${lnls[0]}
    fi
done
echo -e "\n Best k value is: $best_k (replicate 1) with lnL=$best_lnl"
BEST_RESULT_PATH="$CAFE_DIR/k${best_k}_rep1"
```

Can then use cafeplotter to generate summaries (https://github.com/moshi4/CafePlotter).

## Annotate the expanding and contraction gene families
Using EggNog and Interproscan annotations from above to annotate the expanding gene families in each of the species. I used scripts that were generated by Meeran Hussain: https://github.com/meeranhussain/Comparative_genomics. This is done in R.  

## Merge the eggnog + interproscan annotations for each species
```
library(tidyverse)
library(dplyr)


samples <- c("01_hilli", "02_quadrimaculata", "03_stygia", "04_vicina", "05_sericata", "06_megacephala")

for (sample in samples) {
  # File paths
  eggnog_gff <- file.path(sample, paste0(sample, "_eggnog.emapper.decorated_modified.gff"))
  interpro_tsv <- file.path(sample, paste0(sample, "_longest_isoforms_modified.faa.tsv"))
  
  cat("Processing:", sample, "\n")
  
  #---- 1. Parse EggNOG-mapper GFF3 ----#
  eggnog <- read_tsv(
    eggnog_gff,
    comment = "##",
    col_names = FALSE
  )
  
  # Filter mRNA rows where em_* annotations are present
  eggnog_mrna <- eggnog %>% filter(X3 == "mRNA")
  
  #---- 2. Extract key fields from attribute column using regex ----#
  eggnog_clean <- eggnog_mrna %>%
    transmute(Protein_ID = str_extract(X9, "(?<=ID=)[^;]+"),
              GO_terms   = str_extract(X9, "(?<=em_GOs=)[^;]+"),
              KEGG_paths = str_extract(X9, "(?<=em_KEGG_Pathway=)[^;]+"),
              KEGG_KO    = str_extract(X9, "(?<=em_KEGG_ko=)[^;]+"),
              BRITE_class= str_extract(X9, "(?<=em_BRITE=)[^;]+"),
              EggNOG_PFAMs = str_extract(X9, "(?<=em_PFAMs=)[^;]+"),
              EggNOG_PFAMs_desc = str_extract(X9, "(?<=em_desc=)[^;]+"))
  
  # Remove rows with all NA or Protein_ID missing
  eggnog_clean <- eggnog_clean %>%
    filter(!is.na(Protein_ID)) %>%
    group_by(Protein_ID) %>%
    summarise(across(everything(), ~ toString(na.omit(.))), .groups = "drop") %>%
    mutate(
      ID_number = as.numeric(str_extract(Protein_ID, "(?<=g)\\d+"))
    ) %>%
    arrange(ID_number) %>%
    select(-ID_number)
  
  #---- 2. Parse InterProScan TSV ----#
  interpro <- read_tsv(interpro_tsv, col_names = FALSE, comment = "#")
  
  # 1. Pfam domains and descriptions (only from source = "Pfam")
  pfam_data <- interpro %>%
    filter(X4 == "Pfam") %>%
    select(Protein_ID = X1, Pfam_ID = X5, Pfam_desc = X6, Pfam_Eval = X9) %>%
    mutate(Pfam_Eval = as.numeric(Pfam_Eval)) %>%
    group_by(Protein_ID, Pfam_ID, Pfam_desc) %>%
    summarise(Min_Eval = min(Pfam_Eval, na.rm = TRUE), .groups = "drop") %>%
    group_by(Protein_ID) %>%
    summarise(
      Pfam_domains = toString(Pfam_ID),
      Pfam_descriptions = toString(Pfam_desc),
      Pfam_Evals = toString(Min_Eval),
      .groups = "drop"
    )
  
  #### Panther fetch
  panther_data <- interpro %>%
    filter(X4 == "PANTHER") %>%
    select(Protein_ID = X1, Panther_ID = X5, Panther_desc = X6, Panther_Eval = X9) %>%
    mutate(Panther_Eval = as.numeric(Panther_Eval)) %>%
    group_by(Protein_ID, Panther_ID, Panther_desc) %>%
    summarise(Min_Eval = min(Panther_Eval, na.rm = TRUE), .groups = "drop") %>%
    group_by(Protein_ID) %>%
    summarise(
      Panther_ID = toString(Panther_ID),
      Panther_descriptions = toString(Panther_desc),
      Panther_Evals = toString(Min_Eval),
      .groups = "drop"
    )
  
  # 2. InterPro IDs and descriptions (from all sources)
  ipr_data <- interpro %>%
    select(Protein_ID = X1, InterPro_ID = X12, InterPro_desc = X13) %>%
    filter(!is.na(InterPro_ID)) %>%
    group_by(Protein_ID) %>%
    summarise(
      InterPro_IDs = toString(unique(na.omit(InterPro_ID))),
      InterPro_descriptions = toString(unique(na.omit(InterPro_desc))),
      .groups = "drop"
    )
  
  # 3. GO terms (column X14), clean up GO IDs
  go_data <- interpro %>%
    select(Protein_ID = X1, GO = X14) %>%
    filter(!is.na(GO)) %>%
    separate_rows(GO, sep = ",\\s*") %>%                              # Split multiple GO term blocks
    mutate(GO = str_extract_all(GO, "GO:\\d{7}")) %>%                 # Extract only GO:#######
  unnest(GO) %>%
    group_by(Protein_ID) %>%
    summarise(
      GO_terms_interpro = toString(unique(GO)),
      .groups = "drop"
    )
  
  # 4. Merge InterProScan pieces
  interpro_clean <- reduce(
    list(pfam_data, panther_data, ipr_data, go_data),
    full_join,
    by = "Protein_ID"
  ) %>%
    mutate(
      ID_number = as.numeric(str_extract(Protein_ID, "(?<=g)\\d+"))
    ) %>%
    arrange(ID_number) %>%
    select(-ID_number)
  
  # 5. Merge with EggNOG annotations (GO terms remain separate)
  merged_annotations <- full_join(eggnog_clean, interpro_clean, by = "Protein_ID") %>%
    select(
      Protein_ID,
      EggNOG_PFAMs,
      EggNOG_PFAMs_desc,
      KEGG_paths,
      KEGG_KO,
      BRITE_class,
      Pfam_domains,
      Pfam_descriptions,
      Pfam_Evals,
      Panther_ID,
      Panther_descriptions,
      Panther_Evals,
      InterPro_IDs,
      InterPro_descriptions,
      GO_terms,               # from EggNOG
      GO_terms_interpro      # cleaned InterProScan GO terms
    )
  
  # 6. Clean NA and export
  merged_annotations <- merged_annotations %>%
    mutate(across(everything(), ~str_replace_all(., "NA|^, | ,|, NA", "")))
  
  # Reorder columns in the desired order
  merged_annotations <- merged_annotations %>%
    select(
      Protein_ID,
      EggNOG_PFAMs,
      EggNOG_PFAMs_desc,
      GO_terms,             
      KEGG_KO,
      KEGG_paths,
      BRITE_class,
      Pfam_domains,
      Pfam_descriptions,
      Pfam_Evals,
      InterPro_IDs,
      InterPro_descriptions,
      GO_terms_interpro,
      Panther_ID,
      Panther_descriptions,
      Panther_Evals,
    )
  colnames(merged_annotations) <- c("Protein_ID", "EggNOG_PFAMs_domains","EggNOG_PFAMs_descriptions", "EggNOG_GO_terms", "EggNOG_KEGG_KO", "EggNOG_KEGG_paths", "EggNOG_BRITE_class", "Interpro_Pfam_Ids",
                                    "Interpro_Pfam_descriptions", "Interpro_Pfam_Evals", "InterPro_IDs", "InterPro_descriptions", "Interpro_GO_terms", "Panther_ID", "Panther_descriptions", "Panther_Evals")
  
  # ---- 4. Write to file ----
  output_file <- file.path(sample, "merged_protein_annotations.tsv")
  write_tsv(merged_annotations, output_file)
}
```

Annotate the merged annotation file with orthogroups from orthofinder:

```
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# === File paths ===
orthogroups_file <- "Orthogroups.tsv"
og_list_dir <- "Orthogroups.tsv"
annotation_dir <- "01_annotated_genes"

# === Load OrthoFinder orthogroups file ===
orthogroups <- read_tsv(orthogroups_file)

# Create a long-format OG → gene mapping
og_gene_map <- orthogroups %>%
  pivot_longer(-Orthogroup, names_to = "Species", values_to = "Genes") %>%
  mutate(Genes = strsplit(Genes, ", ")) %>%
  unnest(Genes) %>%
  distinct(Orthogroup, Species, Genes)


# === Process each species ===
ann_files <- list.files(annotation_dir, pattern = "_annotated.tsv", full.names = TRUE)

for (ann_file in ann_files) {
  species <- str_replace(basename(ann_file), "_annotated.tsv", "")
  message("Processing species: ", species)
  
  # Load OG list
  og_ids <- unique(og_gene_map$Orthogroup)
  og_subset <- og_gene_map %>%
    filter(`Orthogroup` %in% og_ids & str_detect(Genes, species))  # Filter OG + species
  
  
  # Load merged annotation
  annotation_file <- file.path(annotation_dir, paste0(species, "_annotated.tsv"))
  ann <- read_tsv(annotation_file, show_col_types = FALSE)
  
  # Clean GeneID column name
  colnames(ann)[1] <- "GeneID"
  
  
  # Merge with OG data
  annotated <- og_subset %>%
    rename(GeneID = Genes) %>%
    left_join(ann, by = "GeneID") %>%
    select(
      Orthogroup, GeneID,
      EggNOG_PFAMs_domains,
      EggNOG_PFAMs_descriptions,
      EggNOG_GO_terms,
      EggNOG_KEGG_KO,
      EggNOG_KEGG_paths,
      EggNOG_BRITE_class,
      Interpro_Pfam_Ids,
      Interpro_Pfam_descriptions,
      Interpro_Pfam_Evals,
      InterPro_IDs,
      InterPro_descriptions,
      Interpro_GO_terms,
      Panther_ID, Panther_descriptions, Panther_Evals
    )
  
  # Output file
  out_file <- file.path("02_annotated_OGs", paste0(species, "_annotated_OGs.tsv"))
  write_tsv(annotated, out_file)
  message("[✓] Output written: ", out_file)
}
```
## split the terminal significant gene families using this python script:
```
#!/usr/bin/env python3

from collections import defaultdict
import os

# === Input settings ===
input_file = "Gamma_branch_probabilities.tab"
output_dir = "Species_specific_sig_OGs"
threshold = 0.05

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Dictionary to hold column name -> list of matching FamilyIDs
results = defaultdict(list)

# Read the input file
with open(input_file) as f:
    header = f.readline().strip().split("\t")  # first line = species
    for line in f:
        parts = line.strip().split("\t")
        family_id = parts[0]
        for i in range(1, len(parts)):
            try:
                val = float(parts[i])
                if val <= threshold:
                    results[header[i]].append(family_id)
            except ValueError:
                continue  # skip non-numeric or missing values

# Write each species' significant OGs to its own file
for species, og_list in results.items():
    if og_list:  # only if there are significant OGs
        file_path = os.path.join(output_dir, f"{species}_significant_OGs.txt")
        with open(file_path, "w") as f_out:
            f_out.write("\n".join(og_list))
        print(f"[✓] {species}: {len(og_list)} significant OGs written to {file_path}")
```
## Split these significant gene families into expanding and contracting:

```
# === Libraries ===
library(readr)
library(dplyr)
library(tidyr)
library(stringr)

# === Input files ===
gamma_file <- "Gamma_change.tab"
og_dir <- "01_species_sig_ogs"  # folder containing *_sig_ogs.txt files

# === Read CAFE Gamma_change.tab ===
gamma <- read_tsv(gamma_file)

# === Clean column names ===
colnames(gamma)[1] <- "Orthogroup"

# Clean species names
colnames(gamma) <- str_remove(colnames(gamma), "<.*?>")
colnames(gamma) <- str_replace(colnames(gamma),
                               "05_sericata_protein", "05_sericata")
colnames(gamma) <- str_replace(colnames(gamma),
                               "06_megacephala_protein", "06_megacephala")
colnames(gamma) <- str_replace(colnames(gamma),
                               "01_hilli_protein", "01_hilli")
colnames(gamma) <- str_replace(colnames(gamma),
                               "03_stygia_protein", "03_stygia")
colnames(gamma) <- str_replace(colnames(gamma),
                               "02_quadrimaculata_protein", "02_quadrimaculata")
colnames(gamma) <- str_replace(colnames(gamma),
                               "04_vicina_longest_protein", "04_vicina")

# Remove the model <?> column names
gamma <- gamma[, !colnames(gamma) %in% c("", ".1", ".2", ".3")]


# === Convert to long format ===
gamma_long <- gamma %>%
  pivot_longer(-Orthogroup, names_to = "Species", values_to = "Change") %>%
  mutate(
    Change = as.numeric(Change),
    Species = str_replace(Species, "_protein$", "")
  ) %>%
  filter(Change != 0) %>%
  mutate(Direction = ifelse(Change > 0, "expanded", "contracted"))

# === List species-specific OG files ===
og_files <- list.files(og_dir, pattern = "_sig_ogs.txt", full.names = TRUE)

# === Process each species ===
for (og_file in og_files) {
  
  # Extract species name from file
  species <- str_replace(basename(og_file), "_sig_ogs.txt", "")
  message("Processing species: ", species)
  
  # Read significant OGs
  sig_ogs <- read_lines(og_file)
  
  # Match with gamma data
  changes <- gamma_long %>%
    filter(Species == species, Orthogroup %in% sig_ogs)
  
  # Split into expanded and contracted
  expanded <- changes %>% filter(Direction == "expanded") %>% pull(Orthogroup)
  contracted <- changes %>% filter(Direction == "contracted") %>% pull(Orthogroup)
  
  # Write to files
  write_lines(expanded, file.path(og_dir, paste0(species, "_expanded_OGs.txt")))
  write_lines(contracted, file.path(og_dir, paste0(species, "_contracted_OGs.txt")))
  
  message("[✓] ", species, ": ", length(expanded), " expanded, ", length(contracted), " contracted OGs")
}
```
## Annotate the expanding and contracting gene families based on the Eggnog and Interpro annotations
```
#ANNOTATE the expanded and contracted gene families

library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# === File paths ===
orthogroups_file <- "Orthogroups.tsv"
og_list_dir <- "expanded_OGs"
annotation_dir <- "../01_merge_annotations/02_annotated_OGs"
output_dir <- "/nesi/nobackup/uow03920/06_blowfly_assembly_jan/analysis/05_fetch_species_specific_gene_IDs_and_annotate_GO_terms/03_expanded_contracted_annotated"  

# === Load OrthoFinder orthogroups file ===
orthogroups <- read_tsv(orthogroups_file)

# Create a long-format OG → gene mapping
og_gene_map <- orthogroups %>%
  pivot_longer(-`Orthogroup`, names_to = "Species", values_to = "Genes") %>%
  mutate(Genes = strsplit(Genes, ", ")) %>%
  unnest(Genes)

og_gene_map <- og_gene_map %>%
  mutate(Species = str_replace(Species, "^[0-9]+_", "")) %>%  # remove 01_, 02_, etc
  mutate(Species = str_replace(Species, "_no_dup_longest_isoforms", ""))  # remove suffix


# === Process each species ===
og_files <- list.files(og_list_dir, pattern = "_expanded_OGs.txt", full.names = TRUE)

for (og_file in og_files) {
  species <- str_replace(basename(og_file), "_expanded_OGs.txt", "")
  message("Processing species: ", species)
  
  # Load OG list
  og_ids <- read_lines(og_file)
  og_subset <- og_gene_map %>%
    filter(`Orthogroup` %in% og_ids & str_detect(Genes, species))  # Filter OG + species
  
  
  # Load merged annotation
  annotation_file <- file.path(annotation_dir, paste0(species, "_annotated.tsv"))
  ann <- read_tsv(annotation_file, show_col_types = FALSE)
  
  # Clean GeneID column name
  colnames(ann)[1] <- "GeneID"
  
  
  # Merge with OG data
  annotated <- og_subset %>%
    rename(GeneID = Genes) %>%
    left_join(ann, by = "GeneID") %>%
    select(
      Orthogroup, GeneID,
      EggNOG_PFAMs_domains,
      EggNOG_GO_terms,
      EggNOG_KEGG_KO,
      EggNOG_KEGG_paths,
      EggNOG_BRITE_class,
      Interpro_Pfam_Ids,
      Interpro_Pfam_descriptions,
      InterPro_IDs,
      InterPro_descriptions,
      Interpro_GO_terms
    )
  
  # Output file
  out_file <- file.path(output_dir, paste0(species, "_annotated_expanded_OGs.tsv"))
  write_tsv(annotated, out_file)
}


##################################### Contracted ##################################################################
# === File paths ===
orthogroups_file <- "Orthogroups.tsv"
og_list_dir <- "contracted_OGs"
annotation_dir <- "../01_merge_annotations/02_annotated_OGs"
output_dir <- "/nesi/nobackup/uow03920/06_blowfly_assembly_jan/analysis/05_fetch_species_specific_gene_IDs_and_annotate_GO_terms/03_expanded_contracted_annotated"  

# === Load OrthoFinder orthogroups file ===
orthogroups <- read_tsv(orthogroups_file)

# Create a long-format OG → gene mapping
og_gene_map <- orthogroups %>%
  pivot_longer(-`Orthogroup`, names_to = "Species", values_to = "Genes") %>%
  mutate(Genes = strsplit(Genes, ", ")) %>%
  unnest(Genes)

og_gene_map <- og_gene_map %>%
  mutate(Species = str_replace(Species, "^[0-9]+_", "")) %>%  # remove 01_, 02_, etc
  mutate(Species = str_replace(Species, "_no_dup_longest_isoforms", ""))  # remove suffix


# === Process each species ===
og_files <- list.files(og_list_dir, pattern = "_contracted_OGs.txt", full.names = TRUE)

for (og_file in og_files) {
  species <- str_replace(basename(og_file), "_contracted_OGs.txt", "")
  message("Processing species: ", species)
  
  # Load OG list
  og_ids <- read_lines(og_file)
  og_subset <- og_gene_map %>%
    filter(`Orthogroup` %in% og_ids & str_detect(Genes, species))  # Filter OG + species
  
  
  # Load merged annotation
  annotation_file <- file.path(annotation_dir, paste0(species, "_annotated.tsv"))
  ann <- read_tsv(annotation_file, show_col_types = FALSE)
  
  # Clean GeneID column name
  colnames(ann)[1] <- "GeneID"
  
  
  # Merge with OG data
  annotated <- og_subset %>%
    rename(GeneID = Genes) %>%
    left_join(ann, by = "GeneID") %>%
    select(
      Orthogroup, GeneID,
      EggNOG_PFAMs_domains,
      EggNOG_GO_terms,
      EggNOG_KEGG_KO,
      EggNOG_KEGG_paths,
      EggNOG_BRITE_class,
      Interpro_Pfam_Ids,
      Interpro_Pfam_descriptions,
      InterPro_IDs,
      InterPro_descriptions,
      Interpro_GO_terms
    )
  
  # Output file
  out_file <- file.path(output_dir, paste0(species, "_annotated_contracted_OGs.tsv"))
  write_tsv(annotated, out_file)
}
```

I then manually went through these tsv files and counted how many gene families were associated with the PFAMs that are associated with transposable elements. I created a heat map of these using this script: 

```
library(tidyverse)

df <- tribble(
  ~Pfam,     ~HIL, ~QUA, ~STY, ~VIC, ~SER, ~MEG,
  "PF00075", 1,    0,    0,    0,    5,    7,
  "PF00078", 2,    0,    2,    2,    89,   36,
  "PF00665", 0,    1,    1,    0,    36,   14,
  "PF02925", 0,    0,    0,    0,    0,    0,
  "PF02992", 0,    0,    0,    0,    0,    0,
  "PF03184", 1,    0,    1,    0,    1,    0,
  "PF03221", 1,    0,    0,    0,    0,    0,
  "PF03732", 0,    0,    1,    0,    11,   5,
  "PF04687", 0,    0,    0,    0,    0,    0,
  "PF05699", 1,    0,    1,    0,    0,    1,
  "PF05840", 0,    0,    0,    0,    0,    0,
  "PF05970", 0,    0,    2,    0,    2,    3,
  "PF07727", 0,    0,    1,    0,    7,    7,
  "PF08283", 0,    0,    0,    0,    0,    0,
  "PF08284", 0,    0,    0,    0,    1,    0,
  "PF10551", 0,    0,    0,    0,    0,    0,
  "PF13358", 2,    0,    2,    0,    7,    4,
  "PF13359", 5,    0,    0,    0,    2,    0,
  "PF13456", 0,    0,    0,    0,    0,    0,
  "PF13837", 1,    0,    0,    0,    0,    0,
  "PF13976", 0,    0,    0,    0,    5,    6,
  "PF14214", 0,    0,    1,    0,    4,    3,
  "PF14223", 2,    0,    3,    0,    10,   7,
  "PF14529", 0,    0,    0,    1,    20,   13
)

# Remove rows with all zeros
df_filtered <- df %>%
  filter(rowSums(across(-Pfam)) > 0)

# Convert to long format
df_long <- df_filtered %>%
  pivot_longer(
    cols = -Pfam,
    names_to = "Species",
    values_to = "Count"
  )

df_long$Species <- factor(
  df_long$Species,
  levels = c("QUA", "STY", "HIL", "VIC", "SER", "MEG")
)

# Heatmap
ggplot(df_long, aes(x = Species, y = Pfam, fill = Count)) +
  geom_tile(color = "white") +
  
  scale_fill_gradientn(
    colours = c("#FFC7DE", "#EB71B9", "#3C0D70"),
    trans = "log1p",
    breaks = c(0, 10, 40, 80)
  ) +
  labs(
    x = NULL) +
  
  theme_minimal() +
  theme(
    text = element_text(size = 16, colour = "black"),
    axis.text.x = element_text(size = 16, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    panel.grid = element_blank()
  )

ggsave(
  "TE_gene_family_heatmap.png",
  dpi = 300
)
```
