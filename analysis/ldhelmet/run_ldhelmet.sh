block=$1
outdir=block_$1

for fname in data/references/fasta/filtered/*fasta; do
    base=$(basename $fname .fasta)
    echo "Currently on ${base}"

    time ./bin/ldhelmet find_confs \
    --num_threads 10 \
    --window_size 50 \
    --output_file data/ldhelmet/${outdir}/${base}.conf ${fname}

    sleep 1

    echo "table_gen for ${base}"

    time ./bin/ldhelmet table_gen \
    --num_threads 10 \
    --conf_file data/ldhelmet/${outdir}/${base}.conf \
    --theta 0.03 \
    --rhos 0.0 0.1 10.0 1.0 100.0 \
    --output_file data/ldhelmet/${outdir}/${base}.lk > table_gen 2> table_gen2
    
    sleep 1

    rm table_gen*

    sleep 1

    time ./bin/ldhelmet pade \
    --num_threads 10 \
    --conf_file data/ldhelmet/${outdir}/${base}.conf \
    --theta 0.03 \
    --output_file data/ldhelmet/${outdir}/${base}.pade

    sleep 1

    time ./bin/ldhelmet rjmcmc \
    --num_threads 30 \
    --window_size 50 \
    --seq_file ${fname} \
    --lk_file data/ldhelmet/${outdir}/${base}.lk \
    --pade_file data/ldhelmet/${outdir}/${base}.pade \
    --num_iter 1000000 \
    --burn_in 100000 \
    --block_penalty ${block} \
    --mut_mat_file data/ldhelmet/mut_mat \
    --output_file data/ldhelmet/${outdir}/${base}.post

    sleep 1

    time ./bin/ldhelmet post_to_text \
    --mean \
    --perc 0.025 \
    --perc 0.50 \
    --perc 0.975 \
    --output_file data/ldhelmet/${outdir}/${base}.txt \
    data/ldhelmet/${outdir}/${base}.post

    sleep 3

    echo "Removing temp files..."
    rm data/ldhelmet/${outdir}/${base}.conf
    rm data/ldhelmet/${outdir}/${base}.lk
    rm data/ldhelmet/${outdir}/${base}.pade
    rm data/ldhelmet/${outdir}/${base}.post

done 

