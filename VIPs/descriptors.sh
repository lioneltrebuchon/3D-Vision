rm VIPSifts.txt
for i in *VIP.jpg; do
    echo $i
    ./siftFromImage $i >> VIPSifts.txt
done