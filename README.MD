# AREDCI: Assessing REproducibility and Differential Chromatin Interactions

## Project Overview
AREDCI is a novel approach designed to address the challenges in analyzing ChIA-PET sequencing data. It integrates data preprocessing, normalization, reproducibility assessment, and differential chromatin interaction (DCI) identification to provide robust evaluation of sample consistency and enhance the accuracy and reliability of DCI identification.

## Technologies Used
- **Programming Language**: Python
- **Dependencies**: NumPy, flask, matplotlib, statsmodels, scipy, xgboost, sklearn, networkx, idr

## Installation Guide
### System Requirements
- Python 3.8 or higher

### Installation Steps
1. Clone the repository to your local machine
`git clone https://github.com/lyotvincent/AREDCI.git`

2. Install dependencies according to requirements.txt

## Usage
### Using the Web Interface
1. `python web.py`;
2. Open your browser and go to `http://127.0.0.1:1080/`;
3. Upload path of your ChIA-PET data and customize the parameters;
4. Click the "Run DCI analysis" button to start the analysis process;

### Web Interface
![Web Interface Screenshot](/static/images//web_interface.jpeg)

## Core Features
### Data Preprocessing
AREDCI has a comprehensive preprocessing workflow to enhance data quality.

### Normalization Techniques
AREDCI incorporates multiple normalization techniques including TMM, MA, MD, and SSRD.

### Reproducibility Assessment
AREDCI develops a self-similarity-based algorithm for evaluating inter-sample reproducibility.

### Differential Chromatin Interaction (DCI) Identification
AREDCI designs algorithms based on KDE, local contextual information, MD, and SSRD to enhance the accuracy and reliability of DCI identification.

## Algorithm Details
AREDCI employs several algorithms to improve the accuracy of DCI identification:
- **Kernel Density Estimation (KDE)**: Used to estimate the distribution function of the error ratio.
- **Local Contextual Information**: A network-based approach to analyze DCI by integrating local context features of interactions.
- **Minus, or difference vs. Distance (MD)**: A method based on the ratio of fitted values between sample groups.
- **Signed Squared Relative Difference (SSRD)**: An improved method for identifying DCI between sample groups.

## Experimental Results
AREDCI demonstrates commendable precision, recall, and F1-score in both simulated and real-world datasets, showcasing its effectiveness in analyzing chromatin interaction data.

## Experimental Data
The accession numbers for the data sets used in this study are as follows:  
* GM12878:
  - Replicate1: ENCFF409UJV and ENCFF810FRZ;
  - Replicate2: ENCFF048QWJ and ENCFF062ISF.  
* K562:
  - Replicate1: ENCFF485GHA and ENCFF116MYA;
  - Replicate2: ENCFF915YCD and ENCFF245TLN.  
* HCT116:
  - Replicate1: ENCFF695OOH and ENCFF424OGG;
  - Replicate2: ENCFF869HLZ and ENCFF933EOS.  
* MCF-7:
  - Replicate1: ENCFF391NJE and ENCFF465LQF;
  - Replicate2: ENCFF592KNN and ENCFF342KQM.  
* Mus musculus strain B6CASTF1/J embryonic stem cell male embryo:
  - Replicate1: ENCFF533ZMZ and ENCFF776ZNO;
  - Replicate2: ENCFF399JRD and ENCFF538TLL  
* The Cell Differentiation dataset is available under GEO accession number GSE44067.  