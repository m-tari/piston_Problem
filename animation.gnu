	
set xrange [-10:110]
set yrange [-10:15.5]
set nokey	
set grid

#set terminal pdf
#set output "animation.pdf"

do for [in=1:1600] {plot 'flow.dat' i in u 1:2 w lp lc rgb "black" t columnheader(1) ; pause 0.01}


pause -1
