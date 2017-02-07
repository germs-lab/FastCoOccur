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
./co data_table col#_of_treatment col#_of_starting_data > result.tsv
```
example
```
./co table_for_co_core.csv 18 19 0.05 > co_core_p0.05.tsv
```

### Step 3: filter result
filter by rho
```
python filter_co_result.py result_file.tsv 0.6 > filtered.result.tsv
```
### Step 4: Plot or additional analysis in R
