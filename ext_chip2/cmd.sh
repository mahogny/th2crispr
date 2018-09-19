for f in *.bed
do
        echo $f
	annotatePeaks.pl $f mm10 > $f.annot
done


