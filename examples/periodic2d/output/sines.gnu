# mesh PLOT
reset
set palette defined (-1 "dark-blue", 0 "white", 1 "dark-red")
set cbrange [  -1.0 : 1.0] 
set notics
unset border
set nokey
unset colorbox
set size ratio -1
name = "output/mesh"

set term post eps color solid
set output name.".ps"
set xrange [-1.0 : 1.00]
set yrange [-1.0 : 1.00]
# insert label

plot \
   name."_P00" w image, \
   name."_P01" w image, \
   name."_P02" w image, \
   name."_P03" w image
