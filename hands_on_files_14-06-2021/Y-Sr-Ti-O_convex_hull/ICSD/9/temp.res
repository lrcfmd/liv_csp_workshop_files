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
   3.955158   3.955158   3.955158  90.000000  90.000000  90.000000
fractional  1   
O     core 0.0000000 0.5000000 0.5000000 -2.0000000 1.00000 0.00000             
O     core 0.5000000 1.0000000 0.5000000 -2.0000000 1.00000 0.00000             
O     core 0.5000000 0.5000000 1.0000000 -2.0000000 1.00000 0.00000             
Sr    core 0.0000000 1.0000000 1.0000000 2.00000000 1.00000 0.00000             
Ti    core 0.5000000 0.5000000 0.5000000 4.00000000 1.00000 0.00000             
totalenergy          -158.5835462699 eV
species   3
O      core   -2.000000                  
Sr     core    2.000000                  
Ti     core    4.000000                  
buck     
O     core Sr    core  1952.39000     0.336850  19.22000      1.20 12.00
buck     
O     core Ti    core  4590.72790     0.261000  0.000000      1.20 12.00
buck     
O     core O     core  1388.77000     0.362620  175.0000      1.20 12.00
lbfgs_order 5000
time        300.0
maxcyc opt     1500
maxcyc fit     1500
dump temp.res                                                    
