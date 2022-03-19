reset
set terminal gif animate delay 1
set output 'kapitza.gif'
set xrange [-10:10]
set yrange [-10:10]
set size square 1,1
set pointsize 2

set style line 2 lc rgb '#0060ad' pt 7

do for [i=0:4000] {
    plot 'kapitza.dat' u (0):(0):3:4 every ::2*i::2*i w vectors lc "grey" nohead notitle, \
         x2=y2=NaN '' u (x1=x2,x2=$3):(y1=y2,y2=$4):(x1-x2):(int($0)%2==0 ? NaN: y1-y2) every ::2*i::2*i+1 w vectors lc "grey" nohead notitle, \
         ''    u 3:4 every ::2*i::2*i     w p pt 7 ps 3 lc "red" title "1", \
         ''    u 3:4 every ::2*i+1::2*i+1 w p pt 7 ps 2 lc "blue" title "2", \
}



