Additional information on input data:

- assemblies (reference). 
    
    Should be determined in config/default.yaml. Section "sample" (noted as "parameters").
    
    Extension is determined in config/default.yaml. Section "assemblies_ext" (noted as "Extensions").

- reads
    
    Extension is determined in config/default.yaml. Section "reads_ext" (noted as "Extensions").
    
    Required structure example:

    └── reads

        ├── sample1

            ├── sample1_1.fastq.gz

            └── sample1_2.fastq.gz

        ├── sample2

            ├── sample2_1.fastq.gz

            └── sample2_2.fastq.gz
