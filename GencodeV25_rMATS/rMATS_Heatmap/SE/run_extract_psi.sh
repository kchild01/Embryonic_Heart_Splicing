sc_PSI_count=../scripts/extract_PSI_count.py
declare -a event_array=("SE")
declare -a counttype_array=("JC")
for event in "${event_array[@]}"; do for counttype in "${counttype_array[@]}"; do python $sc_PSI_count ./filtered2_${event}.MATS.${counttype}.txt ./PSI; done
> done
