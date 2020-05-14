SCRIPTSFOLDER=$1
find . -name "*.frip" -not -name "*init*" -exec cat {} \; | python ${SCRIPTSFOLDER}/qc_get_frip_stats_table.py
