#!/bin/bash
if [ $# -lt 5 ]
then
    echo 'not enough args'
else
    matlab -nodisplay -nodesktop -nojvm -r "create_interp_matrix($1, $2, $3, $4, '$5')"
fi