# abdoulie_bojang_01
Assemble and profile Staphylococcus aureus samples from clinical samples

<!-- MarkdownTOC autolink="true" levels="1,2,3" -->

- [Set up analysis](#set-up-analysis)
	- [Set working directory](#set-working-directory)
	- [Trim primer sequences from reads using Trimmomatic](#trim-primer-sequences-from-reads-using-trimmomatic)
	- [Combine single end reads into the same FASTQ](#combine-single-end-reads-into-the-same-fastq)
- [Assemble Staphylococcus genomes](#assemble-staphylococcus-genomes)
	- [Assemble reads for each sample using SPAdes](#assemble-reads-for-each-sample-using-spades)
	- [Subset assemblies to keep only contigs >5k bp](#subset-assemblies-to-keep-only-contigs-5k-bp)
	- [Run a BLASTn search on contigs to identify specifically Staphylococcus contigs](#run-a-blastn-search-on-contigs-to-identify-specifically-staphylococcus-contigs)
		- [Conduct BLASTn search against nt](#conduct-blastn-search-against-nt)
		- [Construct a list of accession IDs from BLASTn hits](#construct-a-list-of-accession-ids-from-blastn-hits)
		- [Map accession hits to taxa names](#map-accession-hits-to-taxa-names)
		- [Map contig hits to taxa names](#map-contig-hits-to-taxa-names)
		- [Create subset assembly that only contains Staphlyococcus contigs](#create-subset-assembly-that-only-contains-staphlyococcus-contigs)
- [Identify sequence types for each Staphylococcus assembly](#identify-sequence-types-for-each-staphylococcus-assembly)
	- [Insert sample identifier into FASTA headers and create a combined FASTA](#insert-sample-identifier-into-fasta-headers-and-create-a-combined-fasta)
	- [Split combined FASTA into 10 FASTAs for MLST querying](#split-combined-fasta-into-10-fastas-for-mlst-querying)
	- [Query split FASTAs for MLST profile](#query-split-fastas-for-mlst-profile)

<!-- /MarkdownTOC -->


# Set up analysis

## Set working directory

```{bash, eval = F}
WORKING_DIR=/data/SGSlab/mchung/abdoulie_bojang_01
```

## Trim primer sequences from reads using Trimmomatic

Based off Abdoulies primer sequences, he is using TruSeq3-PE-2 primers:

Adapter: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
Adapter Read 2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

```{bash, eval = F}
module load trimmomatic/0.39
for SAMPLE in $(find $WORKING_DIR/reads | sed "s/_L001.*//g" | sort -n | uniq)
do
sbatch --time=48:00:00 --mem=245g --wrap="trimmomatic PE \
-phred33 \
-threads 8 \
"$SAMPLE"_L001_R1_001.fastq.gz \
"$SAMPLE"_L001_R2_001.fastq.gz \
$SAMPLE.1.trimmed.fastq.gz \
$SAMPLE.1.se.trimmed.fastq.gz \
$SAMPLE.2.trimmed.fastq.gz \
$SAMPLE.2.se.trimmed.fastq.gz \
ILLUMINACLIP:/usr/local/apps/trimmomatic/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads \
LEADING:3 TRAILING:3 MINLEN:36"
done
```

## Combine single end reads into the same FASTQ

```{bash, eval = F}
for SAMPLE in $(find $WORKING_DIR/reads | grep -v trimmed | sed "s/_L001.*//g" | sort -n | uniq)
do
cat "$SAMPLE".1.se.trimmed.fastq.gz "$SAMPLE".2.se.trimmed.fastq.gz > "$SAMPLE".se.trimmed.fastq.gz
done
```

# Assemble Staphylococcus genomes

## Assemble reads for each sample using SPAdes

```{bash, eval = F}
module load spades/3.15.3

for SAMPLE in $(find $WORKING_DIR/reads | grep -v trimmed | sed "s/_L001.*//g" | sort -n | uniq)
do

sbatch --time=48:00:00 --wrap="spades.py \
--isolate \
-1 "$SAMPLE".1.trimmed.fastq.gz \
-2 "$SAMPLE".2.trimmed.fastq.gz \
-s "$SAMPLE".se.trimmed.fastq.gz \
-o "$WORKING_DIR"/spades/"$(basename "$SAMPLE")""
done
```

## Subset assemblies to keep only contigs >5k bp

```{bash, eval = F}
module load seqkit/2.1.0

for FASTA in $(find $WORKING_DIR/spades -name "contigs.fasta")
do
seqkit seq -m 5000 "$FASTA" > "$(echo "$FASTA" | sed "s/[.]fasta/.l5000.fasta/g")"
done
```

## Run a BLASTn search on contigs to identify specifically Staphylococcus contigs

### Conduct BLASTn search against nt
```{bash, eval = F}
module load blast/2.9.0+

for FASTA in $(find $WORKING_DIR/spades -name "contigs.l5000.fasta")
do
sbatch --cpus-per-task=8 --time=48:00:00 --wrap="blastn \
-query "$FASTA" \
-db /fdb/blastdb/nt \
-evalue 1e-6 \
-outfmt 6 \
-max_target_seqs 1 \
-max_hsps 1 \
-out "$(echo "$FASTA" | sed "s/[.]fasta/.blastn.tsv/g")" \
-num_threads 8"
done
 ```

### Construct a list of accession IDs from BLASTn hits
```{bash, eval = F}
rm $WORKING_DIR/spades/accession_id.list
for BLASTN_OUTPUT in $(find "$WORKING_DIR"/spades -name "*blastn.tsv")
do
cut -f2 "$BLASTN_OUTPUT" >> "$WORKING_DIR"/spades/accession_id.list
done
cat $WORKING_DIR/spades/accession_id.list | sort -n | uniq > "$WORKING_DIR"/spades/temp
mv "$WORKING_DIR"/spades/temp "$WORKING_DIR"/spades/accession_id.list
```

### Map accession hits to taxa names

```{bash, eval = F}
module load edirect/8.60
cat "$WORKING_DIR"/spades/accession_id.list | epost -db nuccore | esummary -db nuccore | xtract -pattern DocumentSummary -element Caption,TaxId > "$WORKING_DIR"/spades/accession2taxid.list
cat "$WORKING_DIR"/spades/accession2taxid.list | cut -f2 | epost -db taxonomy | esummary -db taxonomy | grep "<TaxId>" | sed "s/.*<TaxId>//g" | sed "s/<\\/TaxId>//g" > "$WORKING_DIR"/spades/col1
cat "$WORKING_DIR"/spades/accession2taxid.list | cut -f2 | epost -db taxonomy | esummary -db taxonomy | grep "<Genus>" | sed "s/.*<Genus>//g" | sed "s/<\\/Genus>//g" > "$WORKING_DIR"/spades/col2
cat "$WORKING_DIR"/spades/accession2taxid.list | cut -f2 | epost -db taxonomy | esummary -db taxonomy | grep "<ScientificName>" | sed "s/.*<ScientificName>//g" | sed "s/<\\/ScientificName>//g" > "$WORKING_DIR"/spades/col3

paste "$WORKING_DIR"/spades/col1 "$WORKING_DIR"/spades/col2 "$WORKING_DIR"/spades/col3 > "$WORKING_DIR"/spades/taxid2taxa.map
rm col*

rm "$WORKING_DIR"/spades/accession2taxid2taxa.list
IFS=""
cat "$WORKING_DIR"/spades/accession2taxid.list  |while read LINE
do
ACCESSION=$(echo $LINE | awk -F "\t" '{print $1}')
TAXID=$(echo $LINE | awk -F "\t" '{print $2}')
TAXA=$(awk -v TAXID=$TAXID -F "\t" '$1 == TAXID {print $0}' "$WORKING_DIR"/spades/taxid2taxa.map)

echo -e ""$ACCESSION"\t"$TAXA"" >> "$WORKING_DIR"/spades/accession2taxid2taxa.list

done
```

### Map contig hits to taxa names

```{bash, eval = F}
IFS=$' \t\n'
for BLASTN_OUTPUT in $(find "$WORKING_DIR"/spades -name "*blastn.tsv")
do
rm $(dirname $BLASTN_OUTPUT)/accession2taxa.map
cat $BLASTN_OUTPUT | while read LINE
do
CONTIG=$(echo $LINE | cut -d " " -f1)
ACCESSION=$(echo $LINE | cut -d " " -f2 | sed "s/[.].*//g")
TAXA=$(grep $ACCESSION "$WORKING_DIR"/spades/accession2taxid2taxa.list | awk -F "\t" '{print $3"\t"$4}')
echo -e ""$CONTIG"\t"$ACCESSION"\t""$TAXA" >> $(dirname $BLASTN_OUTPUT)/accession2taxa.map
done
done
```

### Create subset assembly that only contains Staphlyococcus contigs

```{bash, eval = F}
module load samtools/1.13

for MAP in $(find "$WORKING_DIR"/spades -name "accession2taxa.map")
do
awk '$3 == "Staphylococcus" {print $1}' $(dirname "$MAP")/accession2taxa.map > $(dirname "$MAP")/staphylococcus_contigs.list

samtools faidx $(dirname "$MAP")/contigs.fasta -r $(dirname "$MAP")/staphylococcus_contigs.list -o $(dirname "$MAP")/contigs.staphylococcus.fasta
done
```

# Identify sequence types for each Staphylococcus assembly

## Insert sample identifier into FASTA headers and create a combined FASTA 
```{bash, eval = F}
rm "$WORKING_DIR"/spades/combined.fna
for FASTA in $(find "$WORKING_DIR"/spades/ -name "contigs.staphylococcus.fasta")
do
sed -i "s/>/>$(dirname $FASTA | sed "s/.*\\///g")\|/g" $FASTA
cat $FASTA >>"$WORKING_DIR"/spades/combined.fna
done
```

## Split combined FASTA into 10 FASTAs for MLST querying

```{bash, eval = F}
module load ucsc/423
faSplit sequence "$WORKING_DIR"/spades/combined.fna 10 "$WORKING_DIR"/spades/split
```

## Query split FASTAs for MLST profile

MLST profiles were queried on https://pubmlst.org/bigsdb?db=pubmlst_saureus_seqdef&page=batchSequenceQuery
