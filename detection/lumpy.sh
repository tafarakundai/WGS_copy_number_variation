#First conda deactivate so that the python modules do not conflict

module add chpc/BIOMODULES lumpy
module add samtools
module add sambamba/0.7.1
module add samblaster 
module add python/2.7.15 

module add python/3.9.6

# NAVIGATE TO WERE YOUR BAM FILES ARE 

# Step 1: Step 1: Extract Discordant Paired-End Alignments

for bam in *_dedup.bam; do
    base=$(basename "$bam" _dedup.bam)
    echo "Extracting discordant reads for $base"
    samtools view -b -F 1294 "$bam" > "${base}_discordants.unsorted.bam"
done

# Use this if you have some files already processed

for bam in *_dedup.bam; do
    base=$(basename "$bam" _dedup.bam)
    discordant_file="${base}_discordants.unsorted.bam"
    
    if [[ -f "$discordant_file" ]]; then
        echo "Skipping $base - discordant file already exists: $discordant_file"
        continue
    fi
    
    echo "Extracting discordant reads for $base"
    samtools view -b -F 1294 "$bam" > "$discordant_file"
done

# Step 2: Extract Split-Read Alignments

for bam in *_dedup.bam; do
    base=$(basename "$bam" _dedup.bam)
    splitter_file="${base}.splitters.unsorted.bam"
    
    if [[ -f "$splitter_file" ]]; then
        echo "Skipping $base - splitter file already exists: $splitter_file"
        continue
    fi
    
    echo "Extracting split reads for $base"
    samtools view -h "$bam" | \
    /home/apps/chpc/bio/lumpy-0.2.13/scripts/extractSplitReads_BwaMem -i stdin | \
    samtools view -Sb - > "$splitter_file"
done

# use this to specify breeds that might have been left out 

for bam in P*_dedup.bam; do
    base=$(basename "$bam" _dedup.bam)
    splitter_file="${base}.splitters.unsorted.bam"
    
    echo "Extracting split reads for $base"
    samtools view -h "$bam" | \
    /home/apps/chpc/bio/lumpy-0.2.13/scripts/extractSplitReads_BwaMem -i stdin | \
    samtools view -Sb - > "$splitter_file"
done

# Step 3: Sort both bam files 

# Sort discordant BAMs
for file in *_discordants.unsorted.bam; do
    base=$(basename "$file" _discordants.unsorted.bam)
    echo "Sorting discordant BAM for $base"
    samtools sort "$file" -o "${base}.discordants.bam"
done

# Sort splitter BAMs
for file in *.splitters.unsorted.bam; do
    base=$(basename "$file" .splitters.unsorted.bam)
    echo "Sorting splitters BAM for $base"
    samtools sort "$file" -o "${base}.splitters.bam"
done

# Step 4: Indexing both files

for bam in *.discordants.bam *.splitters.bam; do
    echo "Indexing $bam"
    samtools index "$bam"
done

#Separately 
# Index discordant BAM files
for bam in *.discordants.bam; do
    echo "Indexing discordant BAM: $bam"
    samtools index "$bam"
done

# Index splitter BAM files
for bam in *.splitters.bam; do
    echo "Indexing splitter BAM: $bam"
    samtools index "$bam"
done


# RUN LUMPY EXPRESS

#Create your own config file because if you are working on a cluster, the installation might not match the default 
# LUMP configuration 

# Copy the original config
cp /home/apps/chpc/bio/lumpy-0.2.13/scripts/lumpyexpress.config ./my_lumpyexpress.config

# Edit your copy
nano ./my_lumpyexpress.config

OR create it this way

# Create the corrected config file directly
cat > ./my_lumpyexpress.config << 'EOF'
#!/bin/bash -e
# general
LUMPY_HOME=/home/apps/chpc/bio/lumpy-0.2.13
# HEXDUMP is used to determine if a file is a CRAM
HEXDUMP=`which hexdump || true`
LUMPY=/home/apps/chpc/bio/lumpy-0.2.13/bin/lumpy
SAMBLASTER=`which samblaster || true`
# either sambamba or samtools is required
SAMBAMBA=`which sambamba || true`
SAMTOOLS=`which samtools || true`
# python 2.7 or newer, must have pysam, numpy installed
PYTHON=`which python || true`
# python scripts
PAIREND_DISTRO=$LUMPY_HOME/scripts/pairend_distro.py
BAMGROUPREADS=$LUMPY_HOME/scripts/bamkit/bamgroupreads.py
BAMFILTERRG=$LUMPY_HOME/scripts/bamkit/bamfilterrg.py
BAMLIBS=$LUMPY_HOME/scripts/bamkit/bamlibs.py
EOF

Then use your config file with the -K option:

# Run with your custom config
/home/apps/chpc/bio/lumpy-0.2.13/scripts/lumpyexpress \
    -B "IKA1_dedup.bam" \
    -S "IKA1.splitters.bam" \
    -D "IKA1.discordants.bam" \
    -K ./my_lumpyexpress.config \
    -o "IKA1.vcf"

# RUN THE LUMPY EXPRESS
for file in *_dedup.bam; do
    # Extract the first 4 letters as sample name
    sample=$(echo "$file" | cut -c1-4)
    
    echo "Processing sample: $sample"
    echo "Files: ${sample}*_dedup.bam, ${sample}*.splitters.bam, ${sample}*.discordants.bam"
    
    /home/apps/chpc/bio/lumpy-0.2.13/scripts/lumpyexpress \
        -B "${file}" \
        -S "${sample}"*.splitters.bam \
        -D "${sample}"*.discordants.bam \
        -K ./my_lumpyexpress.config \
        -o "${sample}.vcf"
        
    echo "Completed $sample at $(date)"
    echo "---"
done

# GENOTYPING

module add svtyper
# example for genotyping 1 sample
svtyper \
       -i IBR1x.vcf \
       -B IBR1_dedup.bam \
       -o IBR1.gt.vcf \
       --verbose
                  
                  
 # Example for genotyping multiple samples
 
 # Loop through all IHL*.vcf files
 for vcf in IHL*.vcf; do
     # Extract the base name 
     base=$(basename "$vcf" .vcf)
     
     svtyper \
         -i "$vcf" \
         -B /mnt/lustre/users/mmalima/All_projects/AsianZebu/Final_Data2/CNV/CNVcaller/${base}_dedup.bam \
         -o ${base}.gt.vcf \
         --verbose
done

# POST PROCESSING 

# Process all SV types
python svtype.py IBR1.gt.vcf IBR1.txt

# Process only deletions and duplications
python svtype.py  IBR1.gt.vcf IBR1.txt --sv-types DEL DUP

# Process specific SV types
python svtype.py  IBR1.gt.vcf IBR1.txt --sv-types DEL DUP INV BND
                             
python svtype.py  IDN1.gt.vcf IDN1.txt --sv-types DEL DUP INV 
                                                   
# Multiple files   

for file in *.gt.vcf; do
    # Extract the identifier (4th character and beyond)
    base=$(basename "$file" .gt.vcf)
    python3 svtype.py "$file" "${base}.txt" --sv-types DEL DUP 
done                            
                             
# If you have samples that gave you problems and you want to run them seperately run this!!!
echo "Extracting discordant reads for IKA1" && samtools view -b -F 1294 IKA1_dedup.bam > IKA1_discordants.unsorted.bam

samtools sort IKA1_discordants.unsorted.bam -o IKA1.discordants.bam

samtools index IKA1.discordants.bam
