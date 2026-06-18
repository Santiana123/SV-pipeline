#!/bin/bash
# Extract and filter INV breakpoints from merged VCF (SURVIVOR version)

if [ $# -ne 2 ]; then
    echo "Usage: $0 <merge_vcf>"
    exit 1
fi

merge_vcf=$1

awk -F'\t' '
BEGIN {
    OFS = "\t";
    print "CHR", "POS", "ID", "CHR2", "END", "SVLEN", "SOURCE", "Counted";
}
!/^#/ {
    svtype = ""; chr2 = "NA"; end = "NA"; svlen = "NA"; counted = "YES"; supp_vec = "";
    split($8, info, ";");

    for (i in info) {
        if (info[i] ~ /^SVTYPE=/) { split(info[i], a, "="); svtype = a[2]; }
        if (info[i] ~ /^CHR2=/)   { split(info[i], b, "="); chr2 = b[2]; }
        if (info[i] ~ /^END=/)    { split(info[i], c, "="); end = c[2]; }
        if (info[i] ~ /^SUPP_VEC=/) { split(info[i], e, "="); supp_vec = e[2]; }
        if (info[i] ~ /^SVLEN=/)  {
            split(info[i], d, "=");
            svlen = d[2];
            gsub(/[^0-9.-]/, "", svlen);
            svlen = (svlen < 0) ? -svlen : svlen;
        }
    }

    if ((svlen == "NA" || svlen == 0) && end != "NA") {
        diff = (end + 0) - ($2 + 0);
        svlen = (diff < 0) ? -diff : diff;
    }

    if (svtype == "INV") {
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
            if ($3 ~ /^svim_asm/) source = "svim";
            else if ($3 ~ /^pbsv/) source = "pbsv";
            else if ($3 ~ /^Sniffles2/) source = "Sniffles2";
            else source = "unknown";
        }

        if ($3 ~ /BND/ && $1 == chr2 && end != "NA" && ($2 + 0) > (end + 0)) {
            counted = "NO";
        }

        print $1, $2, $3, chr2, end, svlen, source, counted;
    }
}' "$merge_vcf" > INV.breakpoint.txt
