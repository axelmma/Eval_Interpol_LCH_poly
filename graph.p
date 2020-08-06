      set terminal png

      set xtic 300 font ",20" scale 1.5                      	 
      set ytic 100  font ",20"  
      set ytics out scale 1.0 nomirror

    	set xtics out nomirror offset 0,-1,0                      
      set xlabel "$p$" offset 0,-2,0 font ",20"
      set ylabel "Seconds" offset -5,0,0 font ",20"


      set tmargin 0
      set rmargin 0

      set xr [2:1600]
      set yr [0:500]

      set key left top Left width -2 box lw 0.3 height 1.2 spacing 1  font ",15"
      set output "Graphetikz.png"
      set title "multipoint evaluation" font ",30"

      #Change Multi-point_LCH_r according to the chosen extension of the field
      plot "Multi-point_LCH_r_2.dat" using 1:2 smooth acsplines dt 1 lc 2 lw 2.5 title "Algorithm 1",\
           "Multi-point_LCH_r_2.dat" using 1:3 smooth acsplines dt 4 lc 6 lw 2.5 title "Alg. 1 suboptimal.",\
           "Multi-point_LCH_r_2.dat" using 1:4 smooth acsplines dt 2 lc 7 lw 2.5 title "State of art"

