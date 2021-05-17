Example source code to automatically process data from the FDC.  

Usage: ./download_treatment.py test/result.example.txt  
As the result, in the test folder, you will see two new folders:  
1, raw: raw data download from FDC server, including meta data annotation and processed expression matrices (FDC can automatically generate it from public repositories).  
2, diff: differential expression profiles. Each file will also include a count map file, for the number of replicates per condition. If a dataset has Sub Condition annotated, a *.sep.gz will also be generated through merging replicates at Sub Condition levels.  

This example code only works for treatment experiment with annotation fields, Treatment, Condition, Sub Condition, Dose, and Duration. For other curation scheme, you need to design your own program.
