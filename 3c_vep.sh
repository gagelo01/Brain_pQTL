cd /mnt/sda/boujer01/Transcriptomique/ensembl-vep
export PERL5LIB=/home/boujer01/perl5/lib/perl5/x86_64-linux-gnu-thread-multi/
output_name=/home/gagelo01/workspace/Projects/Brain_pQTL/Data/Modified/VEPoutput.txt
input_vcf=/home/gagelo01/workspace/Projects/Brain_pQTL/Data/Modified/VEPinput.vcf
./vep -i $input_vcf --o $output_name --cache --dir_cache ./ --offline --force_overwrite
