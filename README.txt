This repository contains code used in the manuscript:
Mönttinen H.A.M and Löytynoja A. "Template switching in DNA replication can create and maintain RNA hairpins"


Prerequests:

	The code is designed to run in conda, to install conda, please see:
	https://conda.io/projects/conda/en/latest/user-guide/install/index.html


	Also you will need MongoDB installed:

	https://docs.mongodb.com/guides/server/install/


	Create conda environment with the required python packages:

	conda create --name py3 --file req.txt

	
	Activate the environment:
	
	conda activate py3
	
	install ete3:
	
	conda install -c etetoolkit ete3 

	
	Compile the fpa2.cpp file:

	c++ -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` fpa2.cpp -o fpa_ext4`python3-config --extension-suffix`


Running scripts:

	The main scripts can be found from main/ directory. To run them please see the doctring in the
	beginning of each script.

