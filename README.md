## Credits/Citations
This project is forked from Ishan Goel's repository: https://github.com/quackduck/annealer
Then it's built onto it

**Single threaded**
1) main_job.sh runs Makefile that runs annealer.cpp, detanneal.cpp, vertexing.cpp for single threaded SA and DA runs
   In SA, try Souvik's sweep based SA and Ishan's SA algorithms by appropriately renaming destination folders

2) binning_job_py.sh runs AtMakefile_py that runs binning_the_efficiency_results_approach2.py
   This is for binning the convergence efficiencies and finding average anneal time

3) LinearRegressionOnAnnealTimeVsSweeps.ipynb in Jupyter Lab for Time/Anneal vs sweeps using .txt output of binning_the_efficiency_results_approach2.py

**Multi threaded** (Yet to be finalized)
1) main_job_parallel.sh runs Makefile_parallel that runs annealer.cpp or annealer_parallelized_sweeps.cpp, detanneal.cpp, vertexing.cpp for single
   threaded SA and DA runs
   In SA, try Souvik's sweep based SA and Ishan's SA algorithms by appropriately renaming destination folders

2) binning_job_py.sh runs AtMakefile_parallel_py that runs binning_the_efficiency_results_approach2.py
   This is for binning the convergence efficiencies and finding average anneal time

3) LinearRegressionOnAnnealTimeVsSweeps.ipynb in Jupyter Lab for Time/Anneal vs sweeps using .txt output of binning_the_efficiency_results_approach2.py
