rm OriginalSifts.txt
for i in *.jpg; do
    echo -n $i >> OriginalSifts.txt
    echo -n "," >> OriginalSifts.txt
    ./siftFromImage $i >> OriginalSifts.txt
done