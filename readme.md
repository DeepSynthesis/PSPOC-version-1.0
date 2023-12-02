## About The Project

Protein pKa is a fundamental physicochemical parameter that dictates protein structure and function. However, accurately determining protein site-pKa values remains a substantial challenge, both experimentally and theoretically. In this study, we introduce a physical organic approach, leveraging a protein structural and physical-organic-parameter-based representation (P-SPOC), to develop a rapid and intuitive model for protein pKa prediction. Our P-SPOC model achieves state-of-the-art predictive accuracy, with a minimal absolute error (MAE) of only 0.33 pKa units. Furthermore, we have incorporated advanced protein structure prediction models, like AlphaFold2, to approximate structures for proteins lacking three-dimensional representations, which enhances the applicability of our model in the context of structure-undetermined protein research. To promote broader accessibility within the research community, an online prediction interface was also established at isyn.luoszgroup.com.


## Getting Started

### Installation

1. Clone the repo
   ```sh
   git clone https://github.com/DeepSynthesis/PSPOC-version-1.0
   ```
2. Create new python environment
   ```sh
   conda create -n HPSPOC-env python=3.11
   conda activate HPSPOC-env
   pip install rdkit-pypi
   conda install pandas=2.0
   conda install numpy=1.24
   conda install scikit-learn=1.2.2
   conda install xgboost=1.7
   conda install matplotlib=3.7
   ```

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- USAGE EXAMPLES -->
## Usage

1. Descriptor generation for training and test set: 
   ```sh
   cd script
   python PSPOC_web_v3.py
   ```
2. Model training and model test: 
   ```sh
   python xgboostMethod.py
   ```
3. Predict pKa for your own data:
  
   open `./script/pspoc_v1.py` ,then you can find codes as below around line 53.
   ```sh
   if __name__ == '__main__':    

    PATH='./Pred/'
    datafilename='your data file'
   ```
   change line 53 `datafilename='your data file'` to your datafile name. The data (.csv) must include column named `"PDB_ID"` `"Chain_ID"` `"Residue_ID"` `"Residue"` .

   then run
   ```sh
   python pspoc_v1.0.py
   ```
   After prediction, the descriptors of your data will be saved in PredicPSPOC.csv, and the result of your data will be saved in `After******.csv` .

<p align="right">(<a href="#readme-top">back to top</a>)</p>



