#!/bin/bash

# Set variables
ANNOVAR_DIR="/apps/chpc/bio/annovar"
WORK_DIR="/mnt/lustre/users/mmalima/All_projects/AsianZebu/Final_Data2/CNV/ANNOVAR/refgenome"
DB_DIR="$WORK_DIR/annovar_db"
SPECIES="custom"
BUILD="GCF_029378745"

# Change to working directory
cd $WORK_DIR

# Create database directory
mkdir -p $DB_DIR/$SPECIES

# Copy your bosindicus_refGene.txt to the database directory with proper naming
cp $WORK_DIR/bosindicus_refGene.txt $DB_DIR/${SPECIES}/${BUILD}_refGene.txt

# Run ANNOVAR annotation
perl $ANNOVAR_DIR/annotate_variation.pl \
    -geneanno \
    -dbtype refGene \
    -buildver $BUILD \
    -out annotated_CNVs \
    $WORK_DIR/Filtered_CNVRs.avinput \
    $DB_DIR/$SPECIES/

# The output files will be:
# annotated_CNVs.variant_function - contains gene-based annotations
# annotated_CNVs.exonic_variant_function - contains detailed exonic annotations (if any)

echo "Annotation complete!"
echo "Check the following output files:"
echo "- annotated_CNVs.variant_function"
echo "- annotated_CNVs.exonic_variant_function"