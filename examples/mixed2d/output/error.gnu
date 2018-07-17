# mesh PLOT
reset
set palette defined (-1 "dark-blue", 0 "white", 1 "dark-red")

#set cbrange [ -4.0*10.0*exp(-10.0) : 4.0*10.0*exp(-10.0) ] 
#set cbrange [ -exp(-10.0) : exp(-10.0) ] 
#set cbrange [ -0.3*c*exp(-10.0) : 0.3*c*exp(-10.0) ]

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
