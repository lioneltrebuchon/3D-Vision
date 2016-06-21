rm VIPSifts.txt
for i in *VIP.jpg; do
    echo -n $i >> VIPSifts.txt
    echo -n "," >> VIPSifts.txt
    ./siftFromImage $i >> VIPSifts.txt
done