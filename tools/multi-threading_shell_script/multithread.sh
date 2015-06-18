#!/bin/bash
# Yonggang Li, ygli@theory.issp.ac.cn, 12/26/2014
echo "#=======================================================================#"
echo "#              Multi-threading run of IM3D based on shell               #"
echo "#                               *  *  *                                 #"
echo "#                  About 35% acceleration efficiency.                   #"
echo "#           Yonggang Li, ygli@theory.issp.ac.cn, 12/26/2014             #"
echo "#=======================================================================#"

start1=$(date +%s)

echo
sysctl -n machdep.cpu.brand_string
SEND_THREAD_NUM=`sysctl -n machdep.cpu.thread_count`  # mac_osx
echo "Number of threads for this sytem: $SEND_THREAD_NUM"
tmp_fifofile="/tmp/$$.fifo"
mkfifo "$tmp_fifofile"
exec 6<>"$tmp_fifofile"

for ((i=0; i<$SEND_THREAD_NUM; i++)); do
    echo
done >&6

num_ions=`awk -F= 'NR==5 {print $2}' ./Config.in`
echo "Totoal number of ions: $num_ions"
echo
num_ions=$[$num_ions/$SEND_THREAD_NUM]
seed1=`awk -F= 'NR==41 {print $2}' ./Config.in`
seed2=`awk -F= 'NR==42 {print $2}' ./Config.in`

for i in `seq 1 $SEND_THREAD_NUM`; do
    mkdir ./temp$i
    cp im3d ./temp$i/im3d
    cp file1.ply2 ./temp$i/file1.ply2
    cp -r data ./temp$i/data
done

start2=$(date +%s)
for i in `seq 1 $SEND_THREAD_NUM`; do
    read -u6
    {
        echo "Thread $i submitted!"

        seed1=$[$seed1+$i]
        seed2=$[$seed2+$i]
        echo "Seeds for thread $i: $seed1, $seed2"
        sed '5s/.*/max_no_ions='$num_ions'/' ./Config.in >./temp$i/Config1.in
        sed '41s/.*/seed1='$seed1'/' ./temp$i/Config1.in >./temp$i/Config2.in
        sed '42s/.*/seed2='$seed2'/' ./temp$i/Config2.in >./temp$i/Config.in

        mkdir ./temp$i/output
        cd temp$i
        ./im3d >out
        cd ..

        echo "-> Thread $i finished!"

        #sleep 3
        echo >&6
    } &

    #pid=$!
    #echo $pid
done

wait

exec 6>&-
end2=$(date +%s)

# xhzheng, 12/26/2014
./sumfile.sh

echo
echo "Done!"

end1=$(date +%s)
echo "Total consumed time: $[$end1 - $start1] s; Calculation time: $[$end2 - $start2] s."

exit 0
