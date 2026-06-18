#!/bin/bash
# Extract DUP info from merged VCF (SURVIVOR version with multi-source)

if [ $# -ne 2 ]; then
    echo "Usage: $0 <merge_vcf>"
    exit 1
fi

in_vcf=$1

awk -F'\t' '
BEGIN {
    OFS = "\t";
    print "CHR", "POS", "ID", "END", "SVLEN", "TYPE_HINT", "SOURCE";
}
!/^#/ {
    svtype = ""; svlen = "NA"; end = "NA"; type_hint = "DUP"; supp_vec = "";
    split($8, info, ";");

    for (i in info) {
        if (info[i] ~ /^SVTYPE=/) {
            split(info[i], t, "=");
            if (t[2] ~ /DUP/) { svtype = "DUP"; }
        }
        if (info[i] ~ /^END=/)       { split(info[i], a, "="); end = a[2]; }
        if (info[i] ~ /^SUPP_VEC=/) { split(info[i], e, "="); supp_vec = e[2]; }
        if (info[i] ~ /^SVLEN=/)     {                                                                                             
            split(info[i], b, "=");                                                                                                
            svlen = b[2];                                                                                                          
            gsub(/[^0-9.-]/, "", svlen);                                                                                           
            svlen = (svlen < 0) ? -svlen : svlen;                                                                                  
        }                                                                                                                          
    }                                                                                                                              
                                                                                                                                   
    if (svtype == "DUP") {                                                                                                         
        source = "";
        if (supp_vec != "") {
            if (substr(supp_vec, 1, 1) == "1") source = (source == "" ? "SyRI" : source ",SyRI");
            if (substr(supp_vec, 2, 1) == "1") source = (source == "" ? "svim-asm" : source ",svim-asm");
            if (substr(supp_vec, 3, 1) == "1") source = (source == "" ? "pbsv" : source ",pbsv");
            if (substr(supp_vec, 4, 1) == "1") source = (source == "" ? "svim" : source ",svim");
            if (substr(supp_vec, 5, 1) == "1") source = (source == "" ? "cuteSV" : source ",cuteSV");
            if (substr(supp_vec, 6, 1) == "1") source = (source == "" ? "Sniffles2" : source ",Sniffles2");
        }

        if (source == "") {
            if ($3 ~ /^pbsv/) source = "pbsv";
            else if ($3 ~ /^Sniffles2/) source = "Sniffles2";
            else if ($3 ~ /^svim/) source = "svim";
            else source = "unknown";
        }

        if ($3 ~ /DUP_TANDEM/ || $8 ~ /TANDEM/) type_hint = "DUP_TANDEM";
        else if ($3 ~ /INS\.DUP/ || $3 ~ /DUP_INT/) type_hint = "INS+DUP";                                                                           
                                                                                                                                   
        print $1, $2, $3, end, svlen, type_hint, source;                                                                           
    }                                                                                                                              
}' "$in_vcf" > DUP.location.txt
