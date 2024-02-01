for i in 1.30 1.35 1.40 1.45 1.50 1.55 1.60 1.65 1.70 1.80 1.90 2.00
do
  mcl  $1 --abc -I $i -t 8 -o "clr_mcl_clustering.$i.labels"
done