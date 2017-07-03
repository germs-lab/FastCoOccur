# co-occurrence
This repository helps you to run co-occurrence anlaysis

### Step 1 : Prepare table
Your input file needs to be in one table in CSV format with treatment info in front.

### Step 2: run co-occurrence 
compile
```
c++ get-co-occurrence-table.cpp -o co -fopenmp
```

run
```
./co -f data_table -t col#_of_treatment -s col#_of_starting_data -p 0.05> result.tsv
```
-f : input file name
-t : number of comlumn that contain treatment
-s : number of comlumn start actual data
-p : p-value cut off

example
```
./co -f test.unix.csv -t 2 -s 3 -p 0.05
```

### Step 3: filter result
filter by rho
```
python filter_co_result.py result_file.tsv 0.6 > filtered.result.tsv
```
### Step 4: Plot or additional analysis in R
