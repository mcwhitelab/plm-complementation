# After adding a new species
# Only new proteomes and orthologies will be downloaded and formatted


awk -F',' '{print $2}' speckey.txt | tail -n +2 > speclist.txt
bash download_proteomes.sh

cat proteomes/inp*fasta > all_proteomes.fasta

bash download_orthologies.sh
bash format_orthologies.sh



