#!/usr/bin/gawk -f
BEGIN {FS="\t"; OFS="\t"}

/^@/ {print; next} # header lines

{ 
    if ( and($2,4) == 0 ) # is aligned properly
    {
    $1 = "*"
    $7 = "="
    $8 = 0
    $9 = 0
    $11 = "*"   # remove quality values
    NF=11
    print
    }
}

