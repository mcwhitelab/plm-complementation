for f in orthologies/*sqltable; do

       outfile=${f%.sqltable}.pairs
       echo $outfile
       if [[ ! -e $outfile ]]; then

          python ../scripts/format_orthology.py -f $f -o $outfile
       fi
    done
