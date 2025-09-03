while read s; do
        filename=proteomes/inp9_${s}.fasta
        echo $filename
        if [[  ! -e "$filename"  ]]; then
            wget -O "$filename" "https://inparanoidb.sbc.su.se/download/proteome/${s}&prot"
        fi
done < speclist.txt

