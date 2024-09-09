# MoCoProject
In this repo we host code for the work on the ongoing MoCo project at Neurobiology Research Unit, Copenhagen University Hospital (Denmark). Check out the project webpage [here](https://sites.google.com/view/melanieganz/research-projects/imaging-children-without-anesthesia)

Overall project abstract:

It is standard procedure in most Danish hospitals to employ general anesthesia when small children are in need of medical imaging procedures that demand anesthesia. Anesthesia is used to prevent the child from leaving the scanner, to reduce motion artefacts and for dealing with children that are too anxious to enter the scanning environment. Movement during the image acquisition can cause serious distortions of the imaging data, which can invalidate its quality and potentially lead to an erroneous conclusion. Even though complications directly related to general anesthesia are rare, there is an increasing concern about the potential neurotoxic effects of general anesthesia. Additionally, there are logistic and financial challenges associated with the use of general anesthesia.

The aim of our project is to utilize recent developments in medical imaging hardware technology in order to show that most young children can undergo medical imaging procedures without anesthesia and that it is still possible to obtain a high-quality diagnostic scan. We want to test the clinical efficacy of a marker-less motion tracking device that can register the child’s movements while scanning which allows for motion correction of the acquired images.

Our goal is to demonstrate the clinical utility of a novel approach of imaging including training, preparation and the use of the tracking device in children with cerebral palsy that undergo MRI scans for diagnostic purposes. Results from this project can be used when the new dedicated childhood brain damage MRI system funded by the Elsass Foundation is installed at BørneRiget, the new childrens’s hospital to open in Copenhagen in 2024.


## Projects
### ChildrenHeadMotionMRI:
Characterisation of children's head motion for Magnetic Resonance Imaging with and without general anaesthesia
Paper published in Frontiers in Radiology 2021 and accessible [here](https://www.frontiersin.org/articles/10.3389/fradi.2021.789632/full)

### ImageQualityMetrics:
Evaluating the match of image quality metrics with radiological assessment in a dataset with and without motion artifacts
[Abstract presented at ISMRM 2022](https://archive.ismrm.org/2022/2061.html). The data used in thsi work is shared here [OpenNeuro](https://openneuro.org/datasets/ds004332).

### MotionCorrectedClinicalMRProtocol:
Evaluating the performance of markerless prospective motion correction and selective reacquisition in a clinical protocol for brain MRI
Code for the work with a dataset of 22 healthy volunteers to validate the quality of preopsective motion correction techniques prior to the ongoing clinical study at Neurobiology Research Unit, Copenhagen University Hospital (Denmark). Fully anonymized data (BIDSified image data and motion tracking files) are available on [OpenNeuro](https://openneuro.org/datasets/ds004332). Raw data in ISMRM-RD format can be shared upon signing the Open Brain consent and upon request - please e-mail mganz@nru.dk.


### RealNoiseMRI challenge:
MRI reconstruction challenge with realistic noise, first round associated with MICCAI 2021, second round asociated with MedNeurIPS
You can read more about the challenge [here](https://realnoisemri.grand-challenge.org/). Fully anonymized data (BIDSified image data and motion tracking files) are available on [OpenNeuro](https://openneuro.org/datasets/ds004332). Raw data in ISMRM-RD format can be shared upon signing teh Open Brain consent and upon request - please e-mail mganz@nru.dk. 

### Kom med i scanneren - mobile phone app
In order to get access to or be able to re-build the freely available code for the "Kom med i scanneren" - mobile phone app, please have a look at the ReadMe [here](./KomMedIScanneren/README.md).

## Data
### Datasets with and without deliberate head movements for evaluating the performance of markerless prospective motion correction and selective reacquisition in a general clinical protocol for brain MRI
We are sharing 22 datasets with and without deliberate head movements that were used in the projects on ImageQualityMetrics, MotionCorrectedClinicalMRProtocol and the RealNoiseMRI challenge publically on [OpenNeuro](https://openneuro.org/datasets/ds004332). The dataset consists of fully anonymized data (BIDSified image data and motion tracking files). Raw data in ISMRM-RD format can be shared upon signing the Open Brain consent and upon request - please e-mail mganz@nru.dk.

