# Mesoricin Design
This repository contains the code for designing and analyzing mesoricin, a short antifungal peptide sourced from Mesorhizobium sp., to enhance its antifungal activity and reduce toxicity. The design involves three computational steps using the Antifungal Index (AFI), a comprehensive metric that integrates machine learning-predicted minimum inhibitory concentrations (MICs) for multiple representative fungal species.

## Project Overview

### Purpose and Key Steps
This project focuses on optimizing the antifungal peptide Mesoricin through three computational steps using the Antifungal Index (AFI):

1. **Segmentation**: The initial peptide (Mesoricin1) is segmented into smaller peptides, and the variant with the lowest AFI is selected (Mesoricin2).
2. **Single-Point Mutation**: Mesoricin2 is subjected to single-point mutations to identify the variant with the lowest AFI (Mesoricin3).
3. **Global Optimization**: Positions in Mesoricin3 are optimized simultaneously to yield the final variant (Mesoricin4). The search space includes 19 positions, each having 19 possible amino acid substitutions, resulting in a total of 19^19 potential combinations.

### Highlights
- **AFI Framework**: AFI serves as a quantitative metric for guiding peptide optimization by considering antifungal activity while reducing hemolytic and cytotoxic effects. This comprehensive in silico assessment allows for the overall optimization of peptide antifungal activity[1,2].
- **Potential Benefits**: The Antifungal Index (AFI) is developed using publicly available experimental data on antimicrobial peptides. Peptides designed using AFI may possess not only lower MICs but also additional benefits such as high fungicidal activity, potent microbial biofilm inhibition and eradication, stability under extreme conditions, low likelihood of inducing microbial resistance, and minimal hemolytic and cytotoxic effects.
- **Cross-Activity**: Due to the correlation between antimicrobial activities, our extensive previous experiments have shown that AFI is also highly effective in screening antibacterial peptides[3,4,5].

### Usage Instructions
To use this project:

1. **Set up the environment**:
    ```bash
    git clone https://www.github.com/JinZhangLab/mesoricin_design.git
    cd mesoricin_design # Replace with the actual directory path
    conda env create -f environment.yml
    conda activate mesoricin_env
    ```
2. **Run the design script**:
    ```bash
    python ./mesoricin_design.py
    ```

## References

1. J. Zhang, et al. In Silico Design and Synthesis of Antifungal Peptides Guided by Quantitative Antifungal Activity. J. Chem. Inf. Model., 2024, 64 (10): 4277-4285. DOI: 10.1021/acs.jcim.4c00142.
2. J. Zhang, et al. Large-Scale Screening of Antifungal Peptides Based on Quantitative Structure-Activity Relationship. ACS Med. Chem. Lett., 2022, 13(1): 99-104. DOI: 10.1021/acsmedchemlett.1c00556.
3. Yang et al., Antimicrobial peptide DvAMP combats carbapenem-resistant Acinetobacter baumannii infection. International Journal of Antimicrobial Agents, 2024. DOI: 10.1016/j.ijantimicag.2024.107106
4. Yang et al., Novel antimicrobial peptide DvAMP serves as a promising antifungal agent. Bioorganic Chemistry, 2023. DOI: 10.1016/j.bioorg.2023.106679
5. Tian et al., The antibacterial activity and mechanism of a novel peptide MR-22 against multidrug-resistant Escherichia coli. Journal of Antimicrobial Chemotherapy, 2024. DOI: 10.3389/fcimb.2024.1334378

## License

This repository is distributed under the MIT License.

## Contact

For further inquiries, please contact:
- **Email:** jzhang@gmc.edu.cn
