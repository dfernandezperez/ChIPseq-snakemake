mm9=$1
dm6=$2

comm -23 <(samtools view $mm9 | cut -f1 | sort) <(samtools view $dm6 | cut -f1 | sort) > $mm9.uniqID
comm -13 <(samtools view $mm9 | cut -f1 | sort) <(samtools view $dm6 | cut -f1 | sort) > $dm6.uniqID

python scripts/ExtractReads_fromBam.py --bam $mm9 --names $mm9.uniqID --out $mm9.tmp
python scripts/ExtractReads_fromBam.py --bam $dm6 --names $dm6.uniqID --out $dm6.tmp

mv $mm9.tmp $mm9
mv $dm6.tmp $dm6

rm $mm9.uniqID
rm $dm6.uniqID
