Annotated PPI Network
===================
OVERVIEW
--------

HARDWARE REQUIREMENTS
---------------------
The hardware requirements to run the program are:

1. At least 8GB of RAM

2. At least 4 cores in the CPU are recommended

SOFTWARE REQUIREMENTS
---------------------
The software requirements to run the program are:

1. Ubuntu (>=16.04). You can try other OS, but they are not supported

2. Python (>=3.7)

3. All Python libraries specified under requirements.txt

INPUT FILES
-----------
There are 4 files which are necessary to run this program:

1. A protein-protein interaction network in the format of:

    **EntrezID** **EntrezID** **CONFIDENCE**
    
    
    Such as:
    
    216 216 0.76
    
    3679 1134 0.73
    
    55607 71 0.65
    
    5552 960 0.63
    
    2886 2064 0.9
    
    5058 2064 0.73
    
    1742 2064 0.87
    
    An example of such a network can be found [here](http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/download.php)


2. An experiment file, in Excel format. Rows show individual entries, columns are the values of the following properties:
	- **Base mean**
	- **Log fold change**
	- **Statistic**
	- **Adjusted p value**
	- **Ensembl id**
	- **Gene symbol**
	- **Entrez id**

    The file may contain other columns too, but the indices and names of the above columns must be entered to the configuration file.


3. A configuration file. The path to the file is input at the beginning of the program. Details
can be found under README_CONFIGURATION.md

HOW TO RUN
----------
