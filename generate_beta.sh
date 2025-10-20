#! usr/bin/bash

for file in ../raw_data/*.beta ; do
    base=$(basename "$file" .beta)
    echo "Procesando archivo: $base.beta"
    wgbstools view "$file" -r chrX:1-156040895 |
    awk 'BEGIN {OFS="\t"} 
     { 
       if (($4+$5) > 0) 
         beta = $4/($4+$5); 
       else 
         beta = "NA"; 
       print $1, $2, $3, $4, $5, beta 
     }' "$file" > "../processed_data/${base}_processed.bed"

     > "processed_data/${base}_processed.bed"
    
    echo "Archivo procesado guardado como: ../processed_data/${base}_processed.bed"
    echo "-----------------------------------"
done