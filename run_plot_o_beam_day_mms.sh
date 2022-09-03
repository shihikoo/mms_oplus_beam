#!/bin/sh
# This is a comment!
echo Program starts        # This is a comment, too!

start_date="2019-09-17"

for i in {0..364}
do
    idate=`date +"%Y-%m-%d" -d "${start_date} + ${i} + day"`
    echo $idate
    idl -e "plot_o_beam_day_mms, ps=1, store_tplot =1, save_data =1, time_start='"+$idate+"', time_duration = 1"
    cp ~/data/misc/CDFLeapSeconds.txt ~/data/
done

echo Program ends
