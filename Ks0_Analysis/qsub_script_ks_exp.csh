#!/bin/tcsh
if ( -e $GROUP_DIR/group_env.csh ) then
        source $GROUP_DIR/group_env.csh
endif

cd /global/homes/l/lwen1990/pwg/embedding/Ks0/analysis
root4star -b -q cuts_exp_ks.C\(0,700\)
