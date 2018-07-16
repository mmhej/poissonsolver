# mesh PLOT
reset
set palette defined (-1 "dark-blue", 0 "white", 1 "dark-red")

c = 10.0
#set cbrange [ -4.0*c*exp(-c) : 4.0*c*exp(-c) ] 
#set cbrange [ -exp(-c) : exp(-c) ] 
set cbrange [ -0.3*c*exp(-c) : 0.3*c*exp(-c) ]

set notics
unset border
set nokey
unset colorbox
set size ratio -1
name = "output/mesh"

set term post eps color solid
set output name.".ps"
set xrange [-1.0 : 1.0]
set yrange [-1.0 : 1.0]
# insert label

plot \
   name."_P00" w image, \
   name."_P01" w image, \
   name."_P02" w image, \
   name."_P03" w image
