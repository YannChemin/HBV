# HBV
HBV hydrological model[^1], in parallelized C language with OpenMP (very very FAST !)

## Description

The HBV hydrological model[^1], or Hydrologiska Byråns Vattenbalansavdelning model, is a computer simulation used to analyze river discharge and water pollution. Developed originally for use in Scandinavia, this hydrological transport model has also been applied in a large number of catchments on most continents.

This is a translation of the HBV hydrology model in parallelized C language, the speedup is enormous as it uses all of the computer cores ! From hours to seconds...

## Run times

**i7 4.9GHz - 8 cores = 8 threads - HP Omen (20220730, Ubuntu)**

```
real 0m39.214s
user 5m2.379s
sys 0m0.228s
```

**i5 2.5GHz - 4 cores = 8 threads - MacBookPro 2014 (20220730, Ubuntu)**

```
real	1m40,002s
user	9m53,303s
sys	0m1,027s
```

## Example dataset is Karkeh Basin (Iran)[^2][^3] 

```
#********************************************************************************
#*     HBV-7 Karkheh								*
#********************************************************************************
#*1=Doab									*
#*2=Pole_Char									*
#*3=Dobe_Merek									*
#*4=Ghor_B									*
#*5=Darttot_TS									*
#*6=Holilan									*
#*7=Kaka_R_Pole_d_Cham_A							*
#*8=Joligir									*
#*9=Pole_z_Paye_P								*
#*										*
#*- Version 1.2 (3-2-2005): incorporation of sub-basin dependent initial	*
#*	conditions,correction coefficients,elevation differences (between	*
#*	subbasin and station) and surface area fractions forests and fields 	*
#*- Version 1.3.1 (14-2-2005): incorporation of additional model 		*
#*	performance criteria							*
#*- Version 1.3.2 (5-9-2005):	incorporation of output infordmation for all	*
#*	sub-basins and additional performance criterion				*
#*- Version 1.4 (30-9-2005): incorporation of uncertainty analysis for		*
#*	climate change purposes and corresponding output information (see	*
#*	file 10)								*
#********************************************************************************
```

## Compilation & Run

```
make
> gcc -o hbv main.c arrays.c hbv_model.c hbv_report.c hbv_performance.c readcsv.c -lm -fopenmp -Wall
./hbv
```



[^1]: https://en.wikipedia.org/wiki/HBV_hydrology_model
[^2]: Mutuwatte, L., 2005. Calibration of a semi distributed hydrological model using discharge and remote sensing data. PhD Thesis, ITC, The Netherlands.
[^3]: https://github.com/YannChemin/HBV/blob/main/Calibration_of_a_semi_distributed_hydrol.pdf
