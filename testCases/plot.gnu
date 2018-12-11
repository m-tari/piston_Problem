#set xrange [-0.2:1.1]
#set yrange [-0.1:1.4]
#set nokey	
set grid
set key

do for[in=0:100] {plot 'flow.dat' i in u 1:5 w lp t columnheader(1); pause 1}

plot 'flow.dat' i 0 u 1:5 w lp t columnheader(1) ,\
	 'flow.dat' i 1 u 1:5 w lp t columnheader(1)

#unset colorbox

pause -1