### What is CitFinder?

Protein citrullination (or deimination) is an irreversible post-translational modification implicated in several physiological and pathological processes, including gene expression regulation, apoptosis, Rheumatoid arthritis, and Alzheimer's disease. However, challenges in sample preparation and data analysis make it difficult to confidently identify and validate Citrullinated proteins (also known as Citrullinome) and its sites. To overcome these limitations, we generated an algorithm called ‘CitFinder’ to analyze mass spectrometry based spectral libraries to confidently identify and validate citrullinated sites. 

### Dependencies
* python2.7
* Numpy

### Installation
```
git clone https://github.com/Citrullinome/CitFinder.git
```

### Parameters

Usage: -i -o -f -m -g -s -r
Use –h or --help for detailed help for parameter

There are 7 input parameters for CitFinder.py

commond line input	| description
--------------------| ------------
-i, --in	| SpectraST non_consensus_library.splib in txt format
-o, --out	|  Output file of modified peptides pairs with RT information, neutral loss and skyline validation results in csv format
-f, --fasta		| Fasta file required for modification site and 10 amino acid information
-m, --modification	| Please specify one targed mod at a time. For example: R[157] OR R
-g, --grouping	| (Optional) Specify the grouping information and comma seprate them. For example: Heart,Lung,Liver,Muscle,Kidney,Brain. Default will be no grouping
-s, --skyline	| (Optional) Skyline report for validation. Note: the file name and modification mass should be consistent with splib
-r, --rtShift	| (Optional) If rt shift is True, it will only provide the modifed peptide pairs with >= 5 mins rt shift.Otherwise, it will provide all the modified peptide pairs. Default: True

### Running
* Step1: prepare your input files:

  * ```Cit_Mouse_Organs_SpecLib.splib```: Splib file generated from spectrast tool
  * ```Cit_Mouse_Organs_SpecLib_CitFinder_Skyline.csv```: Skyline reports for validation

* Step2: ```python CitFinder.py -i Cit_Mouse_Organs_SpecLib.splib -o Cit_Mouse_Organs_SpecLib_CitFinder_Skyline.csv -g Heart,Lung,Liver,Muscle,Kidney,Brain -f UP_Mouse_Rev_Canonical_20180228_DECOY_iRT.fasta -m R[157] -s Cit_Mouse_Organs_SpecLib_Skyline.csv```

Upon completion, ```Cit_Mouse_Organs_SpecLib_CitFinder_Skyline.csv``` will contain the modified peptides pairs with RT information, neutral loss and skyline validation results.

### Benchmark Datasets
Please refer to the data in the Example folder. We have seperated the data into good, okay and bad categories. And we have put the skyline manual validation spectrums for comparison. For example, to run the good data, simply running the following:
```
python CitFinder.py -i Cit_Good_Examples.splib -o Good_Examples_CitFinder_skyline.csv -f UP_Mouse_Rev_Canonical_20180228_DECOY_iRT.fasta -m R[157] -g Heart,Lung,Liver,Muscle,Kidney,Brain -s skyline_report.csv
```

### Support
If you have any questions about CitFinder, please contact Justyna Fert-Bober (Justyna.Fertbober@cshs.org) or Vidya Venkatraman (vidya.venkatraman@cshs.org)

### Citation
Justyna Fert-Bober, Vidya Venkatraman, Christie Hunter, Ruining Liu, Erin L. Crowgey, Rakhi Pandey, Ronald Holewinski, Alexander Scotland, Ben Berman, Jennifer E. Van Eyk, “Enriched Ion Library for Mouse Citrullinome across multiple organ systems”, Manuscript Submitted (2019)

### Licence
See the [LICENSE](./LICENSE) file for license rights and limitations (Apache2.0).

