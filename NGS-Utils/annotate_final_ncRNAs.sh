#!/bin/bash

# (C) Kranti Konganti
# This program is distributed as Artistic License 2.0
# 06/18/2014
# Coordinate with lncRNApipe output to parse out and add Infernal annotation to final ncRNA transcripts.

# $LastChangedBy: konganti $ =~ m/.+?\:(.+)/;
# $LastChangedDate: 2015-23-02 12:30:27 -0500 (Mon, 23 Feb 2015)  $ =~ m/.+?\:(.+)/;
# $LastChangedRevision: 0.7.2 $ =~ m/.+?\:\s*(.*)\s*.*/;
# $AUTHORFULLNAME = 'Kranti Konganti';

if [  -z "$FINAL_GTF" ]  ||
    [ -z "$CM_TXT_OUT"  ] ||
    [ -z "$CPC_TXT_OUT" ] ||
    [ -z "$COV" ] ||
    [ -z "$INF_GTF" ]; then
    echo -e "\nERROR!\n------\nUsage: $0 <FINAL_GTF> <CM_TXT_OUT> <CPC_TXT_OUT> <COV> <INF_GTF>\n";
    exit -1;
fi

grep noncoding $CPC_TXT_OUT | cut -f 1 | sort -n | uniq | while read unetrid; do 
    trid=${unetrid//\./\\.};
    contig_id=`grep -P "\ttranscript\t.+?\"$trid\".+" $FINAL_GTF | awk '{print $1}'`;
    contig_st=`grep -P "\ttranscript\t.+?\"$trid\".+" $FINAL_GTF | awk '{print $4}'`;
    contig_en=`grep -P "\ttranscript\t.+?\"$trid\".+" $FINAL_GTF | awk '{print $5}'`;
    trlen=`grep -P "\"$trid\"" $FINAL_GTF | grep -oP 'transcript_length \"\d+\"' | head -n 1 | perl -e '\$line = <>; if (\$line =~ m/.+?(\d+)/) {print \$1;}'`;
    hitlen=`grep -P "\s+$trid\s+" $CM_TXT_OUT | sort -k15,15nr | uniq | head -n 1 | awk '{if(\$10=="+") print \$9-\$8; else print \$8-\$9;}'`;
    annot=`grep -P "\s+$trid\s+" $CM_TXT_OUT | sort -k15,15nr | uniq | awk '{$1=""; sc=$15; for(i=3;i<=17;i++) $i=""; if (length($0) != 0) print $0, " | BitScore: ",sc}' | cut -d " " -f 2- | sed -e 's/\s\+/ /g' | head -n 1`;

    if [ -z "$trlen" ] || [ -z "$hitlen" ]; then
	calc_cov=0.0;
    else
	calc_cov=$(echo "scale=2; $hitlen * 100 / $trlen" | bc);
    fi
    
    per_symb='%';

    if [ -z "$annot" ] || [[ $(echo "$calc_cov < $COV" | bc) -eq 1 ]]; then 
	annot='No Annotation'; 
    else
	annot="$annot | SequenceCoverage: $calc_cov$per_symb";
    fi 

    grep_trid="\\\""$trid"\\\"";
    grep -P $grep_trid $FINAL_GTF | sed 's/Infernal.*//' | sed "s/$/Infernal_prediction \"$annot\"\;/";
    
    inf_ann_id=0;
    grep -P "\s+$trid\s+" $CM_TXT_OUT |  while read gene_name gene_id query q_acc mdl mdl_from mdl_to seq_from seq_to strand trunc pass gc bias score e_value inc desc; do
	inf_ann_id=$(($inf_ann_id + 1));
	
	if [ "$strand" == "-" ]; then
	    inf_match_en=$(($seq_from + $contig_st));
	    inf_match_st=$(($seq_to + $contig_st));
	else
	    inf_match_en=$(($seq_to + $contig_st));
	    inf_match_st=$(($seq_from + $contig_st));
	fi
	
	if [ "$inc" == "!" ]; then
	    signi="yes";
	else
	    signi="no";
	fi
	
	echo -e "$contig_id\tlncRNApipe-Infernal\texon\t$inf_match_st\t$inf_match_en\t$score\t$strand\t.\tgene_id \"$gene_id\"; gene_name \"$gene_name\"; transcript_id \"${unetrid}.${inf_ann_id}\"; exon_number \"1\" e_value \"$e_value\"; significant_match \"$signi\"; description \"$desc\";" >> $INF_GTF;
    done;
done;
