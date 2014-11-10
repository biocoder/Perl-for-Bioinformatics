#!/bin/bash

# Kranti Konganti
# 10/30/2013

# $LastChangedBy: konganti $ =~ m/.+?\:(.+)/;
# $LastChangedDate: 2014-11-10 10:05:27 -0500 (Mon, 10 Nov 2014)  $ =~ m/.+?\:(.+)/;
# $LastChangedRevision: 0600 $ =~ m/.+?(\d+)/;
# $AUTHORFULLNAME = 'Kranti Konganti';

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

