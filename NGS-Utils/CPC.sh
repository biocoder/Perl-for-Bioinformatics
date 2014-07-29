#!/bin/bash

# Kranti Konganti
# 10/30/2013

# $LastChangedRevision: 63539213c8805a94fb7070c08b80c941d02c69aa $
# $LastChangedDate: Tue Jul 29 07:56:22 2014 -0500 $
# $LastChangedBy: github.com/biocoder $
# $ChangeLog: $
usage() {
cat <<EOF

Usage: ${0} "CPC_ARG_1 CPC_ARG_2 ..." <path to blastall binary>

Version: 64

Author: Kranti Konganti

EOF

exit;
}

construct_cmd() {
    args='';
    for arg in $1
    do
	args="$args $arg";
    done
}

if [ $# -lt 2 ]; then
    usage;
    exit;
fi


IFS='|';

# Construct CPC call and execute
construct_cmd "$1";
export PATH=$2:$PATH;

cd $CPC_HOME;
cpc_cmd_call="./bin/run_predict.sh$args";

echo $cpc_cmd_call;
eval "$cpc_cmd_call";

