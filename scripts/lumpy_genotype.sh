# Process all PRS samples simultaneously
for sample in PRS1 PRS2 PRS3 PRS4 PRS5; do  # Replace with your actual sample names
    svtyper \
        -i ${sample}.vcf \
        -B /mnt/lustre/users/mmalima/All_projects/AsianZebu/Final_Data2/CNV/CNVcaller/${sample}_dedup.bam \
        -o ${sample}.gt.vcf \
        --verbose &
done

# Wait for all background jobs to complete
wait
echo "All svtyper jobs completed!"