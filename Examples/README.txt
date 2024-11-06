tree: paml-bin/evolver-2(get random rooted tree)-nomber of species 5- no. of trees 1-random number seeds 123-with branch lengths from birth death model- birth rate 0.1- death rate 0.2- sampling frac 0.3- mutation rate 0.0015- theta 0.6-


EX_1: tree stationary_dist -v 0.0015                                                   (change mutation rate)
EX_2: tree stationary_dist -v 0.0015 -fd fitness_dist                                  (input fitness distributions)
EX_3: tree stationary_dist -v 0.0015 -fd fitness_dist -g 3333                          (change genome length)
EX_4: tree stationary_dist -v 0.0015 -fd fitness_dist -r 5e07                          (change recombination rate)
EX_5: tree stationary_dist -v 0.0015 -fd fitness_dist -B                               (with backup)
EX_6: tree stationary_dist -v 0.0015 -fd fitness_dist -n 500                           (change population size)
EX_7: tree stationary_dist -v 0.0015 -fd fitness_dist --nonWF                          (using non-WF)
EX_8: tree none -v 0.0015 -N                                                           (Neutral option)
EX_9: tree stationary_dist -v 0.0015 -fd fitness_dist -S                               (dn/ds)
EX_10: tree stationary_dist -v 0.0015 -fd fitness_dist -m mutation_matrix.csv          (mutation rate matrix)
EX_11: tree stationary_dist -v 0.0015 -fd fitness_dist -d change_parameters.yaml       (varying population parameters)
EX_12: tree stationary_dist -v 0.0015 -fd fitness_dist -k consensus                    (sampling method)
EX_13: tree stationary_dist -v 0.0015 -fd fitness_dist -c                              (counting substitutions) 
