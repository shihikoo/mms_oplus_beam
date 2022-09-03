#!/bin/sh
# This is a comment!

echo start

#idl -e "run_sort_o_beam_map_mms, time_start = '2016-01-01', time_end ='2020-12-31/23:59:59'"  &
#idl -e "run_sort_o_beam_map_mms, time_start = '2016-01-01', time_end ='2020-12-31/23:59:59', reduce = 1" &
#idl -e "run_sort_o_beam_map_mms, time_start = '2019-01-01', time_end ='2019-12-31/23:59:59'"  &
#idl -e "run_sort_o_beam_map_mms, time_start = '2019-01-01', time_end ='2019-12-31/23:59:59', reduce = 1" &
#idl -e "run_sort_o_beam_map_mms, time_start = '2019-04-17', time_end ='2019-08-17'"  &
#idl -e "run_sort_o_beam_map_mms, time_start = '2019-04-17', time_end ='2019-08-17', reduce = 1" &

#idl -e "run_sort_o_beam_map_mms, time_start = '2016-01-01', time_end ='2016-12-31/23:59:59'"  &
#idl -e "run_sort_o_beam_map_mms, time_start = '2017-01-01', time_end ='2017-12-31/23:59:59'"  &
#idl -e "run_sort_o_beam_map_mms, time_start = '2018-01-01', time_end ='2018-12-31/23:59:59'" &
#idl -e "run_sort_o_beam_map_mms, time_start = '2020-01-01', time_end ='2020-12-31/23:59:59'" &

idl -e "run_sort_o_beam_map_mms, time_start = '2016-01-01', time_end ='2020-12-31/23:59:59', avoid_2019 = 1"  &
