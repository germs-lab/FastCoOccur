# co-occurrence
This repository helps you to run co-occurrence anlaysis

### Step 1 : Prepare table
Your input file needs to be in one table in CSV format with treatment info in front.

See convert_data_format.R


### Step 2: run co-occurrence 
compile
```
c++ get-co-occurrence-table.cpp -o FastCoOccur -fopenmp
```

run
```
./FastCoOccur -f data_table -t col#_of_treatment -s col#_of_starting_data -p 0.05> result.tsv
```
* -f : input file name
* -t : number of comlumn that contain treatment
* -s : number of comlumn start actual data
* -p : p-value cut off

example
```
./FastCoOccur -f sample_input_file.csv -t 2 -s 3 -p 0.05 > sample_result_p0.05.tsv
```

### Step 3: filter result
filter by rho
```
python filter_co_result.py sample_result_p0.05.tsv 0.6 > sample_result_p0.05_r0.6.tsv
```
### Step 4: Plot or additional analysis in R


## Run on MSU's HPCC
Here example is provided to run co-occurence in MSU's HPCC. Copy example qsub file into your working directory
```
cp /mnt/research/germs/softwares/co-occurrence/example.qsub .
```
Change filename, then submit job
```
qsub example.qsub
```



### Tips for MSU's HPCC
####to install R library 

Make directory into your home directory such as;
```
~/R/library
```

Then add export into .bashrc file. 
```
echo "export R_LIBS_USER=$/mnt/home/YOUR_ID/R/library" >> ~/.bashrc
```

Load module
```
module load R/3.2.0
```
#### Use multiple core
```
export OMP_NUM_THREADS=10
```
#### compile c++ using multi-core
```
module load GNU/4.9
module load OpenMPI/1.10.0
```
