# 
# Keywords:
# 
opti conp c6  
# 
# Options:
# 
title
ASE calculation                                                                 
end
cell
   4.510281   4.510281   3.067830  90.000000  90.000000  90.000000
fractional  1   
O     core 0.3057000 0.3057000 0.0000000 -2.0000000 1.00000 0.00000             
O     core 0.6940319 0.6940319 0.0000000 -2.0000000 1.00000 0.00000             
O     core 0.1940319 0.8056999 0.5000000 -2.0000000 1.00000 0.00000             
O     core 0.8056999 0.1940319 0.5000000 -2.0000000 1.00000 0.00000             
Ti    core 0.9998659 0.9998659 0.0000000 4.00000000 1.00000 0.00000             
Ti    core 0.4998659 0.4998659 0.5000000 4.00000000 1.00000 0.00000             
totalenergy          -247.2832267005 eV
species   2
O      core   -2.000000                  
Ti     core    4.000000                  
buck     
O     core Ti    core  4590.72790     0.261000  0.000000      1.20 12.00
buck     
O     core O     core  1388.77000     0.362620  175.0000      1.20 12.00
lbfgs_order 5000
time        300.0
maxcyc opt     1500
maxcyc fit     1500
dump temp.res                                                    
