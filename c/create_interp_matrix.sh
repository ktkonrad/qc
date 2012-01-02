#!/bin/bash
if [ $# -lt 5 ]
then
    echo 'usage: create_interp_matrix k dx M upsample outfile'
else
    matlab -nodisplay -nodesktop -nojvm -r "create_interp_matrix($1, $2, $3, $4, '$5')" >/dev/null 2>&1
fi
