# co-occurrence
co-occurrence

### run co-occurrence 
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

### filter result
filter by rho
```
python filter_co_result.py result_file.tsv 0.6 > filtered.result.tsv
```
