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

# QC of all assemblies using QUAST, BUSCO, and KAT 

# busco (do a loop for this if multiple assemblies) 
```
#!/bin/bash -e
#SBATCH --account=uow03920
#SBATCH --job-name=compleasm_purged
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paige.matheson14@gmail.com
#SBATCH --output compleasm_%j.out    # save the output into a file
#SBATCH --error compleasm_%j.err     # save the error output into a file

# purge all other modules that may be loaded, and might interfare
module purge

## load tools
module load compleasm/0.2.5-gimkl-2022a

### Compleasm
compleasm.py run -a /nesi/nobackup/uow03920/05_blowfly_assembly_march/19_scaffold/02_quadrimaculata/02_quadrimaculata_scaffold.fasta -o /nesi/nobackup/uow03920/05_blowfly_assembly_march/19_scaffold/00_busco/02_quadrimaculata_compleasm -l diptera_odb10 -t 8;
```

# quast 
```
#!/bin/bash -e

#SBATCH --account=uow03920
#SBATCH --job-name=compleasm_purged
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paige.matheson14@gmail.com
#SBATCH --output compleasm_%j.out    # save the output into a file
#SBATCH --error compleasm_%j.err     # save the error output into a file

quast.py -t 16 -o /nesi/nobackup/uow03920/05_blowfly_assembly_march/19_scaffold/00_quast -l '01_hilli, 02_quadrimaculata, 03_stygia, 04_vicina' /nesi/nobackup/uow03920/05_blowfly_assembly_march/19_scaffold/01_hilli/01_hilli_scaffold.fasta /nesi/nobackup/uow03920/05_blowfly_assembly_march/19_scaffold/02_quadrimaculata/02_quadrimaculata_scaffold.fasta /nesi/nobackup/uow03920/05_blowfly_assembly_march/19_scaffold/03_stygia/03_stygia_scaffold.fasta /nesi/nobackup/uow03920/05_blowfly_assembly_march/19_scaffold/04_vicina/04_vicina_scaffold.fasta
```

# kat

```
#!/bin/bash -e
#SBATCH --account=uow03920
#SBATCH --job-name=kat_analysis
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=2
#SBATCH --mem=50G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paige.matheson14@gmail.com
#SBATCH --output kat_analysis_%j.out    # save the output into a file
#SBATCH --error kat_analysis_%j.err     # save the error output into a file

# purge all other modules that may be loaded, and might interfere
module purge
ml KAT/2.4.2-gimkl-2018b-Python-3.7.3

########## Loop ##############
for i in 01_hilli 02_quadrimaculata 03_stygia 04_vicina;
do
    cd /nesi/nobackup/uow03920/05_blowfly_assembly_march/19_scaffold/00_kat
    #### Link files
    #assembly
    ln -s /nesi/nobackup/uow03920/05_blowfly_assembly_march/19_scaffold/${i}/${i}_scaffold.fasta
    ####kat command
    kat comp -t 16 -o ${i}_kat "/nesi/nobackup/uow03920/05_blowfly_assembly_march/06_concatenate/${i}.fastq" ${i}_scaffold.fasta
    cd ../../../
done
```


# BlobToolKit for QC and detecting contamination

BlobToolKit is used here for quality control and to detect potential contamination in the genome assembly. Several files need to be generated for BlobToolKit to function effectively:

Blast output (to identify taxonomic affiliations)
Aligned Nanopore reads (BAM file with .csi index)
Assembly (FASTA format)
BUSCO report (CSV format)

## BLAST the assembly against the NT database

```
#!/bin/bash -e
#SBATCH --account=uow03744
#SBATCH --job-name=blastn
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=46
#SBATCH --mem=120G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=meeranhussain1996@gmail.com
#SBATCH --output Blastn_%j.out
#SBATCH --error Blastn_%j.err

# Purge any pre-loaded modules
module purge

# Load BLAST modules
module load BLASTDB/2024-07 BLAST/2.13.0-GCC-11.3.0

# Run BLAST for each scaffolded genome
for i in 06 07 08 13 16 19 40; do 
    blastn \
    -query ./Maethio_${i}/ragtag.scaffold.fasta \
    -task megablast \
    -db nt \
    -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
    -culling_limit 10 \
    -num_threads 42 \
    -evalue 1e-3 \
    -out ./Maethio_${i}/Maethio_${i}_megablast.out;
done
```

# blob to check for contamination

Blob uses the older version of busco so need to run the old one. 

`ml busco`
`ml miniprot`

```busco -i /nesi/nobackup/uow03920/05_blowfly_assembly_march/19_scaffold/01_hilli/01_hilli_scaffold.fasta -l diptera_odb10 -o 01_hilli -m genome --out_path ./01_hilli```

Then run blobtools. You need to download the taxdump file from NCBI and put it in the directory with a link to it. 

pip install blobtools

```
#!/bin/bash -e
#SBATCH --account=uow03920
#SBATCH --job-name=blob
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=25G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paige.matheson14@gmail.com
#SBATCH --output blob_%j.out
#SBATCH --error blob_%j.err

blobtools create --busco /nesi/nobackup/uow03920/05_blowfly_assembly_march/23_busco_for_blob/01_hilli/01_hilli/run_diptera_odb10/full_table.tsv --fasta /nesi/nobackup/uow03920/05_blowfly_assembly_march/19_scaffold/01_hilli/01_hilli_scaffold.fasta --cov /nesi/nobackup/uow03920/05_blowfly_assembly_march/20_contamination/01_minimap/01_hilli.bam --hits /nesi/nobackup/uow03920/05_blowfly_assembly_march/20_contamination/01_hilli_megablast.out --replace --taxdump /nesi/nobackup/uow03920/05_blowfly_assembly_march/22_blob/new_taxdump 01_hilli

```

to view these data is ............. boarderline impossible ..... i had the worst time imaginable trying to figure it out. It ended up being kinda easy but idk. it gives you a cool snail plot though.

download the blobtoolkit-api and blobtoolkit-viewer from https://github.com/genomehubs/blobtoolkit/releases. Put them in your working directory. Make them executable (```chmod +x blobtoolkit-api blobtoolkit-viewer```). Use nesi's ondemand virtual desktop. In the virtual desktop, open a terminal and cd to the working directory (with the api and viewer and data needed; `cd /nesi/nobackup/uow03920/05_blowfly_assembly_march/22_blob`). Then add the api to your path: `export=$PATH:$(pwd)`. Then create the viewer `blobtools host --port 8081 --api-port 9001 --hostname localhost /01_hilli`. Then open the link it sends you in the browser within virtual envivronment. 


download the taxonomy table to use to filter out contaminants!!!

# filtration of scaffolded genome assembly using meeran's custom script!

This script removes contaminants and filters out contigs that are less than 1000bp in length (script provided at https://github.com/meeranhussain/Genome_assembly_AND_annotation). I altered the script because i had long read data and meeran's script used qualimap (i used seqkit instead to filter out reads less than 1000bp). 

This is to remove low quality, small and contaminated contigs. It actually didn't really make a huge deal to my assmeblies but good to do anyway.

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

run by using the following slurm script: 

```
#!/bin/bash -e
#SBATCH --account=uow03920
#SBATCH --job-name=filter_contaminants
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paige.matheson14@gmail.com
#SBATCH --output filter_contaminants_%j.out
#SBATCH --error filter_contaminants_%j.err

module purge 
ml SeqKit

for i in 01_hilli 02_quadrimaculata 03_stygia 04_vicina; do 
./filter.sh -i /nesi/nobackup/uow03920/05_blowfly_assembly_march/19_scaffold/${i}/${i}_scaffold.fasta -b /nesi/nobackup/uow03920/05_blowfly_assembly_march/23_busco_for_blob/${i}/${i}/run_diptera_odb10/full_table.tsv -f ${i}.csv -o ${i}_filtered_contaminants;
done
```

then do busco, quast, etc. to check the quality of the assemblies and make sure that they meet the desired assembly standards.... YAY you are done with assemblies, now time to annotate xx

# annotate! 
## rename the contigs with a sequential naming scheme that adheres to EDTA's requirements
First - you need to rename the contig IDs because EDTA requires the input FASTA file to have contig IDs that are no longer than 13 characters. 

use the following code (make a .sh file)

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

the first fasta is the one we one to rename contigs for and the second is the new output. this script will also produce a mapping file in case we need to see what the original names of the contigs were. 

## identify the repeat regions in the genomes

We want to identify repeat regions in our genomes eg TEs, simple microsatellite repeats, low complexity regions, duplications) so that we can mask them in our downstream gene prediction/annotation analyses. 

with anno --1, we are soft masking, so changing repeat regions to lowercase letters. 


```
  GNU nano 5.6.1                                                                        repeat_annotation.sl                                                                                  
#!/bin/bash -e
#SBATCH --account=uow03920
#SBATCH --job-name=EDTA_03
#SBATCH --time=25:00:00
#SBATCH --cpus-per-task=36
#SBATCH --mem=50G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paige.matheson14@gmail.com
#SBATCH --output=EDTA_03_%j.out
#SBATCH --error=EDTA_03_%j.err

# Purge conflicting modules
module purge

# Load EDTA module
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


# repeat masker
This code masks the repeats that we identified from the above code in a fasta genome assembly. This will avoid the prediction of false positive gene structures in repetitive and low complexitiy regions.

```
#!/bin/bash -e
#SBATCH --account=uow03920
#SBATCH --job-name=Repeatmasker
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=50G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paige.matheson14@gmail.com
#SBATCH --output=Repeatmasker_%j.out
#SBATCH --error=Repeatmasker_%j.err

# Purge conflicting modules
module purge

# Load RepeatMasker
ml RepeatMasker/4.1.0-gimkl-2020a

# Loop through selected strains and run RepeatMasker
for i in 01_hilli 03_stygia; do
    cd ${i}
    RepeatMasker -pa 32 -xsmall -lib ${i}_final.fasta.mod.EDTA.TElib.fa ${i}_filtered_contaminants.fasta
    cd ../
done
```

# gene annotation
we will annotate genomes on the repeat-masked genome using homology-based evidence from protein databases, leveraging BRAKER3 pipelines.

A soft link was created to the masked genome file *_scfld_fil_mod.fasta.masked and renamed as *_masked.fasta within the respective species folder

Protein evidence used for annotation was a concatenated database comprising:
Uniprot-SwissProt proteins
OrthoDB diptera proteins (https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/)
Busco proeins


## combine the two proteins into one database
```
cat Arthropoda.fasta Uniprot_sprot.fasta > proteins.fasta
```

I struggled a bit trying to download braker3. I tried to download singularity, genemark, millions of Perl dependencies but all I needed to run was `ml Apptainer` then `singularity build braker3.sif docker://teambraker/braker3:latest
` lol

have to make sure AUGUSTUS will work since there are some annoying stuff to do with permissions. create a folder for each species and run these code in each of them (i found it easier to run parallel rather than in a loop):
```
#Make a directory in your project to hold AUGUSTUS config
mkdir -p /nesi/nobackup/uow03920/05_blowfly_assembly_march/28_annotation/01_hilli/augustus_config

#module load AUGUSTUS and then add it to the right path so that you have a writeable version in your project folder. 
ml AUGUSTUS
cp -r $AUGUSTUS_CONFIG_PATH/* /nesi/nobackup/uow03920/05_blowfly_assembly_march/28_annotation/01_hilli/augustus_config/
```

nowwwww, you can run the code for the annotation (hopefully)

```
#!/bin/bash -e
#SBATCH --account=uow03920
#SBATCH --job-name=braker_03
#SBATCH --time=40:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=90G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paige.matheson14@gmail.com
#SBATCH --output=braker_%j.out
#SBATCH --error=braker_%j.err

#Clean environment
module purge

#Load apptainer module
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
#!/bin/bash -e
#SBATCH --account=uow03920
#SBATCH --job-name=annotation_qc
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=5
#SBATCH --mem=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paige.matheson14@gmail.com
#SBATCH --output=annotation_qc_%j.out
#SBATCH --error=annotation_qc_%j.err

# purge all other modules that may be loaded, and might interfare
module purge

## load tools
ml BUSCO
module load compleasm/0.2.5-gimkl-2022a

LINEAGE=/nesi/nobackup/uow03920/05_blowfly_assembly_march/29_annotation_qc/diptera_odb10

# Run BUSCO on BRAKER3-predicted protein sequences for each strain
compleasm.py protein -p /nesi/nobackup/uow03920/05_blowfly_assembly_march/28_annotation/01_hilli/braker/braker.aa -o /nesi/nobackup/uow03920/05_blowfly_assembly_march/29_annotation_qc/01_hilli -l $LINEAGE -t 12
```

my proteome assemblies had quite a bit of duplications so I ran a bit of extra code to keep only the longest isoforms, and then repeat busco to check that it made duplication go down. 

NOTES (A gene can often produce more than one RNA transcript, which can then be translated into proteins. Each different transcript of the same gene is called an isoform. Why multiple isoforms? Alternative splicing: exons are combined in different ways.Alternative transcription start sites or polyadenylation sites. Effect: Different isoforms can produce proteins with slightly different sequences or functions.). We keep only the longest isoform (i.e. the most complete protein) to avoid redundancy. Then repeat busco to make sure it reduces duplication.

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

`
# gene validator

```
#!/bin/bash -e
#SBATCH --account=uow03920
#SBATCH --job-name=genevalidator_hilli
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=150G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paige.matheson14@gmail.com
#SBATCH --output=genevalidator_%j.out
#SBATCH --error=genevalidator_%j.err

# Clean environment
module purge


# Load BLAST module (required for homology search)
ml BLAST/2.16.0-GCC-12.3.0

# Optional: Create BLAST database (do once before running the loop)
makeblastdb -in ./uniprot_sprot.fasta -dbtype prot -parse_seqids

DB=$(readlink -f ./uniprot_sprot.fasta)

for i in 01_hilli 02_quadrimaculata 03_stygia 04_vicina; do
  /nesi/nobackup/uow03920/05_blowfly_assembly_march/31_genevalidator/genevalidator/bin/genevalidator \
    -d "$DB" \
    --num_threads 2 \
    -m 8 \
    -o ${i}/02_genevalidator_uni_swissprot \
    ../30_clean_duplicate_proteins/${i}/annotation_longest.aa
done
```

# functional annotation using egg nog

```
#!/bin/bash -e
#SBATCH --account=uow03920
#SBATCH --job-name=eggnog
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=15G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --output=eggnog_%j.out
#SBATCH --error=eggnog_%j.err

# Clean environment
module purge

# Load EggNOG-mapper module
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



