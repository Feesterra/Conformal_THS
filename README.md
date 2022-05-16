# Conformal prediction models for MIEs in the thyroid hormone system

The models are part of the [RiskMix project](https://www.aces.su.se/riskmix/) and described in Dracheva et al.[^1].

#### Dependencies

- Python (version 3.5)
- numpy (version 1.16.6)
- scipy (version 0.18.1)
- pandas (version 0.21.1)
- scikit-learn (version 0.18.2)
- nonconformist

#### Usage

Models are stored in model_files folder and are available via the link:
        https://drive.google.com/file/d/1yTfS9bLOEY3jUdDRTsyc_dSzKsXBb9aE/view?usp=sharing  
Archive password: conformal

1. Thyroid-specific models:
    - TRHR antagonist;
    - TSHR agonist/antagonist;
    - NIS inhibition;
    - DIO1,2,3 inhibition;
    - TPO inhibition;
    - TR&#946; agonist/antagonist;
    - TTR binding.
2. General toxicity models:
    - AhR;
    - CAR agonists/antagonist;
    - PPARD agonist/antagonist;
    - PPARG agonist/antagonist;
    - PXR agonist.

The predictions for all the abovementioned MIEs simultaneously are obtained using the CP_predict script.

Help on CP_predict:

```
Predicts compound activities with CP models

optional arguments:
  -h, --help            show this help message and exit
  -f FILE, --file FILE  File with chemicals to predict (should contain column
                        with SMILES structures named "SMILES")
  -dir DIRECTORY, --directory DIRECTORY
                        Directory for saving the predictions
  -si SIGNIFICANCE, --significance SIGNIFICANCE
                        Sets significance level (1-confidence). If set,
                        results will be converted into the prediction regions.
                        If not provided, the results will be returned in the
                        form of p-values.
  -n NAME, --name NAME  Name of the file with results. If not provided, the
                        datestamp in the form of YYYYMMDD_hm will be used as a
                        filename
```

File with structures to be predicted has to be in the form of .csv document with two columns: chemical IDs and SMILES.
Chemical IDs column can contain any unique identifiers of the molecules; the second column, however, has to be named
        "SMILES" and contain curated chemical structures in SMILES format. To obtain accurate predictions, users are
        strongly recommended to follow the standardization procedure described in Dracheva et al.[^1]

[^1]: Dracheva, E.; Norinder, U.; Rydén, P.; Engelhardt, J.; Weiss, J. M.; Andersson, P. L. _In Silico_ Identification of Potential Thyroid Hormone System Disruptors among Chemicals in Human Serum and Chemicals with a High Exposure Index. _Environ. Sci. Technol._ 2022
DOI: 10.1021/acs.est.1c07762
