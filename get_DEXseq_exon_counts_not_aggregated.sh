#NOTE THAT THIS ANALYSIS USES BAMS PROVIDED BY GTAC WHICH ARE BASED ON REF GENOME: "Ensembl_R76_primary/STAR"
#"R76" of Ensembl implies use of REF GENOME: GRCm38.p2 
#For this reason, the most recent Ensembl annotations (without producing new BAMs) that can be used are: Ensembl 102: Nov 2020 (GRCm38.p6)

export BASEDIR=/storage1/fs1/alberthkim/Active/data/hao_usp7_project/20220120_dexseq_analysis 
export OUTDIR=$BASEDIR/exon_counts/not_aggregated
export BAMDIR=/storage1/fs1/alberthkim/Active/data/hao_usp7_project/20220120_USP7_culture_RNAseq/files/pd6ZzJdl/Kim_s5703_MGI2530

#1. Create a GFF file for use with DEXSeq
cd $BASEDIR/refs
wget http://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz
gunzip Mus_musculus.GRCm38.102.gtf.gz

isub -m 10 -i 'docker(malachig/htseq-0.12.4)'
python /usr/local/DEXSeq/dexseq_prepare_annotation.py --aggregate no $BASEDIR/refs/Mus_musculus.GRCm38.102.gtf $BASEDIR/exon_counts/not_aggregated/Mus_musculus.GRCm38.102.gff
exit

#2. Create annotations from the GFF to help annotate the practically useless exon count file that will be produced in the next step
cd $BASEDIR/exon_counts/not_aggregated/
grep -v aggregate_gene Mus_musculus.GRCm38.102.gff | perl -ne 'chomp; $gid="ZZZ"; $eid="ZZZ"; if ($_ =~ /exonic\_part\_number\s+\"(\w+)\".*gene\_id\s+\"(\w+)\"/){$eid=$1; $gid=$2;} print "$gid:$eid\t$_\n"' | sort > Mus_musculus.GRCm38.102.annotations.tsv

#NOTE: Initially I was NOT understanding that the DEXseq package in R will reassociate counts with annotations using the GFF file.  So the step above is probably not needed.

#3. Prepare exon counts using DEXseq-count and HTSeq-count
#you need to inform the script about unstranded data by specifying the option -s no. If your library preparation protocol reverses the strand (i.e., reads appear on the strand opposite to their gene of origin), use -s reverse. In case of paired-end data, the default (-s yes) means that the read from the first sequence pass is on the same strand as the gene and the read from the second pass on the opposite strand (“forward-reverse” or “fr” order in the parlance of the Bowtie/TopHat manual) and the options -s reverse specifies the opposite case.

#details from the collaborator on library prep: ...

#contruct exon counting commands, one for each sample BAM file
bsub -G compute-oncology -q oncology -g /mgriffit/small -M 20000000 -R 'select[mem>20000] span[hosts=1] rusage[mem=20000]' -oo $OUTDIR/div6_cd_1.exon_counts.log -a 'docker(malachig/htseq-0.12.4)' 'python /usr/local/DEXSeq/dexseq_count.py --stranded no --paired yes --format bam --order pos $OUTDIR/Mus_musculus.GRCm38.102.gff $BAMDIR/div6_cd_1.AACTCCGATC-ATCCGTTGGC/div6_cd_1.AACTCCGATC-ATCCGTTGGC.genome_accepted_hits.bam $OUTDIR/div6_cd_1.exon_counts.tsv'

bsub -G compute-oncology -q oncology -g /mgriffit/small -M 20000000 -R 'select[mem>20000] span[hosts=1] rusage[mem=20000]' -oo $OUTDIR/div6_cd_2.exon_counts.log -a 'docker(malachig/htseq-0.12.4)' 'python /usr/local/DEXSeq/dexseq_count.py --stranded no --paired yes --format bam --order pos $OUTDIR/Mus_musculus.GRCm38.102.gff $BAMDIR/div6_cd_2.GTCCGAGGAA-CTTCAGGTAA/div6_cd_2.GTCCGAGGAA-CTTCAGGTAA.genome_accepted_hits.bam $OUTDIR/div6_cd_2.exon_counts.tsv'

bsub -G compute-oncology -q oncology -g /mgriffit/small -M 20000000 -R 'select[mem>20000] span[hosts=1] rusage[mem=20000]' -oo $OUTDIR/div6_cd_3.exon_counts.log -a 'docker(malachig/htseq-0.12.4)' 'python /usr/local/DEXSeq/dexseq_count.py --stranded no --paired yes --format bam --order pos $OUTDIR/Mus_musculus.GRCm38.102.gff $BAMDIR/div6_cd_3.CGACGTACAC-TTACGTTCTC/div6_cd_3.CGACGTACAC-TTACGTTCTC.genome_accepted_hits.bam $OUTDIR/div6_cd_3.exon_counts.tsv'

bsub -G compute-oncology -q oncology -g /mgriffit/small -M 20000000 -R 'select[mem>20000] span[hosts=1] rusage[mem=20000]' -oo $OUTDIR/div6_wt_1.exon_counts.log -a 'docker(malachig/htseq-0.12.4)' 'python /usr/local/DEXSeq/dexseq_count.py --stranded no --paired yes --format bam --order pos $OUTDIR/Mus_musculus.GRCm38.102.gff $BAMDIR/div6_wt_1.AAGTATCACG-ATCCTCTTGG/div6_wt_1.AAGTATCACG-ATCCTCTTGG.genome_accepted_hits.bam $OUTDIR/div6_wt_1.exon_counts.tsv'

bsub -G compute-oncology -q oncology -g /mgriffit/small -M 20000000 -R 'select[mem>20000] span[hosts=1] rusage[mem=20000]' -oo $OUTDIR/div6_wt_2.exon_counts.log -a 'docker(malachig/htseq-0.12.4)' 'python /usr/local/DEXSeq/dexseq_count.py --stranded no --paired yes --format bam --order pos $OUTDIR/Mus_musculus.GRCm38.102.gff $BAMDIR/div6_wt_2.AGGTGAGAGC-CTGCACACTG/div6_wt_2.AGGTGAGAGC-CTGCACACTG.genome_accepted_hits.bam $OUTDIR/div6_wt_2.exon_counts.tsv'

bsub -G compute-oncology -q oncology -g /mgriffit/small -M 20000000 -R 'select[mem>20000] span[hosts=1] rusage[mem=20000]' -oo $OUTDIR/div6_wt_3.exon_counts.log -a 'docker(malachig/htseq-0.12.4)' 'python /usr/local/DEXSeq/dexseq_count.py --stranded no --paired yes --format bam --order pos $OUTDIR/Mus_musculus.GRCm38.102.gff $BAMDIR/div6_wt_3.AATTGGCACA-CTTGGAGCGA/div6_wt_3.AATTGGCACA-CTTGGAGCGA.genome_accepted_hits.bam $OUTDIR/div6_wt_3.exon_counts.tsv'

bsub -G compute-oncology -q oncology -g /mgriffit/small -M 20000000 -R 'select[mem>20000] span[hosts=1] rusage[mem=20000]' -oo $OUTDIR/div7_cd_1.exon_counts.log -a 'docker(malachig/htseq-0.12.4)' 'python /usr/local/DEXSeq/dexseq_count.py --stranded no --paired yes --format bam --order pos $OUTDIR/Mus_musculus.GRCm38.102.gff $BAMDIR/div7_cd_1.AGACCGAATG-CTGGTCGTTC/div7_cd_1.AGACCGAATG-CTGGTCGTTC.genome_accepted_hits.bam $OUTDIR/div7_cd_1.exon_counts.tsv'

bsub -G compute-oncology -q oncology -g /mgriffit/small -M 20000000 -R 'select[mem>20000] span[hosts=1] rusage[mem=20000]' -oo $OUTDIR/div7_cd_2.exon_counts.log -a 'docker(malachig/htseq-0.12.4)' 'python /usr/local/DEXSeq/dexseq_count.py --stranded no --paired yes --format bam --order pos $OUTDIR/Mus_musculus.GRCm38.102.gff $BAMDIR/div7_cd_2.AGGTCAGCGA-TACACGAAGA/div7_cd_2.AGGTCAGCGA-TACACGAAGA.genome_accepted_hits.bam $OUTDIR/div7_cd_2.exon_counts.tsv'

bsub -G compute-oncology -q oncology -g /mgriffit/small -M 20000000 -R 'select[mem>20000] span[hosts=1] rusage[mem=20000]' -oo $OUTDIR/div7_cd_3.exon_counts.log -a 'docker(malachig/htseq-0.12.4)' 'python /usr/local/DEXSeq/dexseq_count.py --stranded no --paired yes --format bam --order pos $OUTDIR/Mus_musculus.GRCm38.102.gff $BAMDIR/div7_cd_3.CAGTTCATGG-CTTGGACCTG/div7_cd_3.CAGTTCATGG-CTTGGACCTG.genome_accepted_hits.bam $OUTDIR/div7_cd_3.exon_counts.tsv'

bsub -G compute-oncology -q oncology -g /mgriffit/small -M 20000000 -R 'select[mem>20000] span[hosts=1] rusage[mem=20000]' -oo $OUTDIR/div7_wt_1.exon_counts.log -a 'docker(malachig/htseq-0.12.4)' 'python /usr/local/DEXSeq/dexseq_count.py --stranded no --paired yes --format bam --order pos $OUTDIR/Mus_musculus.GRCm38.102.gff $BAMDIR/div7_wt_1.CCGACACATC-GACATGTAGA/div7_wt_1.CCGACACATC-GACATGTAGA.genome_accepted_hits.bam $OUTDIR/div7_wt_1.exon_counts.tsv'

bsub -G compute-oncology -q oncology -g /mgriffit/small -M 20000000 -R 'select[mem>20000] span[hosts=1] rusage[mem=20000]' -oo $OUTDIR/div7_wt_2.exon_counts.log -a 'docker(malachig/htseq-0.12.4)' 'python /usr/local/DEXSeq/dexseq_count.py --stranded no --paired yes --format bam --order pos $OUTDIR/Mus_musculus.GRCm38.102.gff $BAMDIR/div7_wt_2.CGATCAGAAG-AGTAATCGAC/div7_wt_2.CGATCAGAAG-AGTAATCGAC.genome_accepted_hits.bam $OUTDIR/div7_wt_2.exon_counts.tsv'

bsub -G compute-oncology -q oncology -g /mgriffit/small -M 20000000 -R 'select[mem>20000] span[hosts=1] rusage[mem=20000]' -oo $OUTDIR/div7_wt_3.exon_counts.log -a 'docker(malachig/htseq-0.12.4)' 'python /usr/local/DEXSeq/dexseq_count.py --stranded no --paired yes --format bam --order pos $OUTDIR/Mus_musculus.GRCm38.102.gff $BAMDIR/div7_wt_3.GTATGGCTAA-GCCAGCTATC/div7_wt_3.GTATGGCTAA-GCCAGCTATC.genome_accepted_hits.bam $OUTDIR/div7_wt_3.exon_counts.tsv'

bsub -G compute-oncology -q oncology -g /mgriffit/small -M 20000000 -R 'select[mem>20000] span[hosts=1] rusage[mem=20000]' -oo $OUTDIR/div8_cd_1.exon_counts.log -a 'docker(malachig/htseq-0.12.4)' 'python /usr/local/DEXSeq/dexseq_count.py --stranded no --paired yes --format bam --order pos $OUTDIR/Mus_musculus.GRCm38.102.gff $BAMDIR/div8_cd_1.CCTAAGCTGA-ACCAGATTCA/div8_cd_1.CCTAAGCTGA-ACCAGATTCA.genome_accepted_hits.bam $OUTDIR/div8_cd_1.exon_counts.tsv'

bsub -G compute-oncology -q oncology -g /mgriffit/small -M 20000000 -R 'select[mem>20000] span[hosts=1] rusage[mem=20000]' -oo $OUTDIR/div8_cd_2.exon_counts.log -a 'docker(malachig/htseq-0.12.4)' 'python /usr/local/DEXSeq/dexseq_count.py --stranded no --paired yes --format bam --order pos $OUTDIR/Mus_musculus.GRCm38.102.gff $BAMDIR/div8_cd_2.TCAGTCGGCT-GCAGATGTAA/div8_cd_2.TCAGTCGGCT-GCAGATGTAA.genome_accepted_hits.bam $OUTDIR/div8_cd_2.exon_counts.tsv'

bsub -G compute-oncology -q oncology -g /mgriffit/small -M 20000000 -R 'select[mem>20000] span[hosts=1] rusage[mem=20000]' -oo $OUTDIR/div8_cd_3.exon_counts.log -a 'docker(malachig/htseq-0.12.4)' 'python /usr/local/DEXSeq/dexseq_count.py --stranded no --paired yes --format bam --order pos $OUTDIR/Mus_musculus.GRCm38.102.gff $BAMDIR/div8_cd_3.AATCCGAGTA-TTAGAATCGG/div8_cd_3.AATCCGAGTA-TTAGAATCGG.genome_accepted_hits.bam $OUTDIR/div8_cd_3.exon_counts.tsv'

bsub -G compute-oncology -q oncology -g /mgriffit/small -M 20000000 -R 'select[mem>20000] span[hosts=1] rusage[mem=20000]' -oo $OUTDIR/div8_wt_1.exon_counts.log -a 'docker(malachig/htseq-0.12.4)' 'python /usr/local/DEXSeq/dexseq_count.py --stranded no --paired yes --format bam --order pos $OUTDIR/Mus_musculus.GRCm38.102.gff $BAMDIR/div8_wt_1.CCATCTGTCG-CTCTCATGCG/div8_wt_1.CCATCTGTCG-CTCTCATGCG.genome_accepted_hits.bam $OUTDIR/div8_wt_1.exon_counts.tsv'

bsub -G compute-oncology -q oncology -g /mgriffit/small -M 20000000 -R 'select[mem>20000] span[hosts=1] rusage[mem=20000]' -oo $OUTDIR/div8_wt_2.exon_counts.log -a 'docker(malachig/htseq-0.12.4)' 'python /usr/local/DEXSeq/dexseq_count.py --stranded no --paired yes --format bam --order pos $OUTDIR/Mus_musculus.GRCm38.102.gff $BAMDIR/div8_wt_2.AGATCGCTCG-GACGCATTGA/div8_wt_2.AGATCGCTCG-GACGCATTGA.genome_accepted_hits.bam $OUTDIR/div8_wt_2.exon_counts.tsv'

bsub -G compute-oncology -q oncology -g /mgriffit/small -M 20000000 -R 'select[mem>20000] span[hosts=1] rusage[mem=20000]' -oo $OUTDIR/div8_wt_3.exon_counts.log -a 'docker(malachig/htseq-0.12.4)' 'python /usr/local/DEXSeq/dexseq_count.py --stranded no --paired yes --format bam --order pos $OUTDIR/Mus_musculus.GRCm38.102.gff $BAMDIR/div8_wt_3.CCGTGAACTG-TATGCAACTC/div8_wt_3.CCGTGAACTG-TATGCAACTC.genome_accepted_hits.bam $OUTDIR/div8_wt_3.exon_counts.tsv'

#Create clean versions of the DEXseq count files without the statistics lines at the bottom
cd $OUTDIR
grep -v "^_" div6_cd_1.exon_counts.tsv > div6_cd_1.exon_counts.clean.tsv
grep -v "^_" div6_cd_2.exon_counts.tsv > div6_cd_2.exon_counts.clean.tsv
grep -v "^_" div6_cd_3.exon_counts.tsv > div6_cd_3.exon_counts.clean.tsv
grep -v "^_" div6_wt_1.exon_counts.tsv > div6_wt_1.exon_counts.clean.tsv
grep -v "^_" div6_wt_2.exon_counts.tsv > div6_wt_2.exon_counts.clean.tsv
grep -v "^_" div6_wt_3.exon_counts.tsv > div6_wt_3.exon_counts.clean.tsv
grep -v "^_" div7_cd_1.exon_counts.tsv > div7_cd_1.exon_counts.clean.tsv
grep -v "^_" div7_cd_2.exon_counts.tsv > div7_cd_2.exon_counts.clean.tsv
grep -v "^_" div7_cd_3.exon_counts.tsv > div7_cd_3.exon_counts.clean.tsv
grep -v "^_" div7_wt_1.exon_counts.tsv > div7_wt_1.exon_counts.clean.tsv
grep -v "^_" div7_wt_2.exon_counts.tsv > div7_wt_2.exon_counts.clean.tsv
grep -v "^_" div7_wt_3.exon_counts.tsv > div7_wt_3.exon_counts.clean.tsv
grep -v "^_" div8_cd_1.exon_counts.tsv > div8_cd_1.exon_counts.clean.tsv
grep -v "^_" div8_cd_2.exon_counts.tsv > div8_cd_2.exon_counts.clean.tsv
grep -v "^_" div8_cd_3.exon_counts.tsv > div8_cd_3.exon_counts.clean.tsv
grep -v "^_" div8_wt_1.exon_counts.tsv > div8_wt_1.exon_counts.clean.tsv
grep -v "^_" div8_wt_2.exon_counts.tsv > div8_wt_2.exon_counts.clean.tsv
grep -v "^_" div8_wt_3.exon_counts.tsv > div8_wt_3.exon_counts.clean.tsv
wc -l *.clean.tsv *annotations.tsv

