filename2="$(basename $0 .sh)"

col=3
#col=2

if [ $col == "2" ]; then
  ylabel="ln{/Symbol P}"
fi
if [ $col == "3" ]; then
  ylabel="U/N"
  col="(\$3/\$1)"
fi

gnuplot << EOF
set terminal postscript eps enhanced color "Arial"
set encoding iso_8859_1
set size 0.65,0.55
set output "${filename2}.eps"
set xlabel "N"
set ylabel "$ylabel"
set key top left
p 'colMat' u 1:$col w lines lt -1 lc -1 noti
EOF

  #produce jpegs
  gs -sDEVICE=jpeg -dJPEGQ=100 -dNOPAUSE -dBATCH -dSAFER -r300 -sOutputFile=${filename2}.jpg ${filename2}.eps
  mogrify -trim -resize 1800x1600 ${filename2}.jpg
  display -resize 600x600 ${filename2}.jpg &

