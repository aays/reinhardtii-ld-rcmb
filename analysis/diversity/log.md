
## 12/7/2019

today: 
1. port CO density analysis from old folder
2. redo LDhelmet rho ~ CO density correlation
3. clean up CO density ~ diversity code
4. convert CO density to cM/Mb in both of the above cases

first - need to create a new version of the Liu breakpoints file 
(`data/diversity/all_break_points.txt`) that has mean rho for each CO

```bash
time python3.5 analysis/diversity/get_COs_rho.py \
--filename data/diversity/all_break_points.txt \
--table data/correlates/annotation_table_rho.txt.gz \
--out data/diversity/all_break_points_rho.txt
```

then, use annotation table w/ rho values to get distribution of rho:

```bash
zcat data/correlates/annotation_table_rho.txt.gz \
| grep -v '#' \
| awk '{print $(NF-0)}' \
| sort -n - \
| uniq -c > data/diversity/rho_counts.txt
```

while this and the above script run, can get cracking on porting the CO density ~ diversity
analysis since it doesn't involve the rho data
