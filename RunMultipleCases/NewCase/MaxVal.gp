#!/usr/bin/gnuplot
p 'time-n199' binary format='%6lf' u 4:5 w d
min_y = GPVAL_DATA_Y_MIN
max_y = GPVAL_DATA_Y_MAX
min_x = GPVAL_DATA_X_MIN
max_x = GPVAL_DATA_X_MAX
Th = 0.25*(abs(min_y)+abs(max_y)+abs(min_x)+abs(max_x))
print "theta = ", Th
