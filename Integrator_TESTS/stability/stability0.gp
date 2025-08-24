set terminal epslatex size 5,4 standalone color colortext 10 font ",12" header "\\usepackage{amsmath} \\usepackage{amsfonts} \\usepackage[T1]{fontenc} \\usepackage[utf8]{inputenc}"
set output "stability0.tex"

# Define the grey levels
set style line 1 lc rgb "#90EE90" pt 7 ps 1.3  # Lightest green
set style line 2 lc rgb "#000000" pt 7 ps 1.3  # Steel blue
set style line 3 lc rgb "#6495ED" pt 7 ps 1.3  # Cornflower blue
set style line 4 lc rgb "#8B4513" pt 7 ps 1.3  # Saddle brown


# Set the titles and labels for the first plot
#set title "Log-log plot for nu = 0.0"
set xlabel "$\\Delta t / (N_\\text{acc}\\Delta t_c)$"
set ylabel "$E_\\text{stability}$"
set logscale y
#set format x "$10^{%L}$"
set format y "$10^{%L}$"

# Configure ticks
unset mxtics
unset mytics

# Set the key (legend) properties
set key top left
#unset key 

# Plot the data from the files for the first model (nu = 0.0)
#set label "Euler" at graph 0.4,0.72 center rotate by 16
#set label "velocity-Verlet" at graph 0.45,0.44 center rotate by 26
#set label "RKN4" at graph 0.2,0.15 center rotate by 36
#set label "Beeman" at graph 0.4,0.2 center rotate by 45
plot [0.002:1][0.1:1e10] "Euler0.txt" u 1:4 title "Euler" with linespoints linestyle 1, \
                  "VelocityVerlet0.txt" u 1:4 title "Velocity Verlet" with linespoints linestyle 2, \
                  "Beeman0.txt" u ($1/2):4 title "Beeman" with linespoints linestyle 3, \
                  "RK40.txt" u ($1/4):4 title "RKN4" with linespoints linestyle 4

