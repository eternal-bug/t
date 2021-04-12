# =========== configuration ==========
# WARN: don't add space between value and "="
# sequence similarity
export similarity_p=0.70
# sequence continuous align length
export length_p=0.60
# if a sequence in the database contain more than one region alignment,
#   keep_1=1 will only keep the longest alignment;
#   keep_1=0 will accpet all pass above limit.
export keep_1=1
# the length of gene upstream will be extract.
export upstream=200
# discard limit
# if the sequence only near several sequence will be discard
export isolate=3
