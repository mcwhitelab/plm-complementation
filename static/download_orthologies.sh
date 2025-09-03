while read f; do
    while read g; do
        
        if [[ "$f" < "$g" ]]; then
            
            filename=orthologies/inp9_${f}_${g}_prot.sqltable
            echo $filename
            if [[  ! -e "$filename"  ]]; then
                wget -O "$filename" "https://inparanoidb.sbc.su.se/download/sqltable/${f}&${g}&prot"
            fi

        fi
    done < speclist.txt
done < speclist.txt

