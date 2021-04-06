# CPMOD

Main class to run experiments is mtree.tests.OD_Test.java with the arguments like:

--algorithm cpmod_rd --R 1.9 --W 10000 --k 50 --slide 500 --datafile tao.txt --numberWindow 100 --samplingTime 100 

where cpmod_rd is the algorithm name, R is default radius threshold, k is default neighbor count threshold, W is the default window size, and slide is the default slide sizes. All these parameters can be changed in main file OD_Test.java for multiple queries. numberWindow is the number of micro-slide, and samplingTime is frequency of CPU and memory sampling. 


