# Data and Code for "Observation of tandem running behavior in dealates of Asian dampwood termite, _Hodotermopsis sjostedti_"
 
## Article information

This repository provides access to the data and source code used for the manuscript  
**Observation of tandem running behavior in dealates of Asian dampwood termite, _Hodotermopsis sjostedti_**

**Authors:**  
**Nobuaki Mizumoto**<sup>1,2</sup>, **William Chambliss**<sup>1</sup>, **Carroll P Elijah**<sup>1</sup>, **Tomohiro Nakazono**<sup>3</sup>, **Taisuke Kanao**<sup>4</sup>  

<sup>1</sup> Department of Entomology & Plant Pathology, Auburn University, Auburn, AL, 36849, USA<br>
<sup>2</sup> Okinawa Institute of Science and Technology, Onna-son, Okinawa, 904-0495, Japan<br>
<sup>3</sup> Laboratory of Insect Ecology, Graduate School of Agriculture, Kyoto University, Kyoto, Japan 606-8502<br>
<sup>4</sup> Faculty of Science, Yamagata University, Yamagata 990-8560, Japan<br>

**Paper DOI:**   
Paper accepted at Journal of Ethology [10.1007/s10164-025-00857-5](https://doi.org/10.1007/s10164-025-00857-5)  
Preprint at bioRxiv [10.1101/2025.05.02.651862](https://doi.org/10.1101/2025.05.02.651862)  

This study describes the tandem running behavior in the termite *Hodotermopsis sjostedti*, using high-resolution position data extracted from video tracking. The analysis focuses on identifying tandem runs, measuring their durations, detecting leader roles (male or female), and analyzing how behaviors vary with arena size.  
This repository includes raw tracking data and the Python and R scripts to analyze them.  

## Repository Structure

- **`/analysis/codes/`**: Contains scripts for analysis.
  - **`analysis.R`**: R script for conducting statistical analysis and generating figures.
  - **`sleap_processing.py`**: Python script to process and clean data from SLEAP tracking outputs.
  
- **`/data_raw/`**: Contains raw `.h5` data produced by SLEAP.
- **`/data_fmt/`**: Contains formatted datasets generated during analysis.
- **`/output/`**: Contains output files, including figures and analysis results.

## Setup & Dependencies


This project is written in R. You’ll need the following packages:

```r
install.packages(c("stringr", "data.table", "arrow", "dplyr", "MASS", "ggplot2",
                   "patchwork", "knitr", "survival", "survminer", "zoo",
                   "cowplot", "coxme", "tidyr"))
```
This project also uses Python. You’ll need the following Python packages:

```ini
pandas==1.3.5
h5py==3.1.0
numpy==1.19.5
scipy==1.7.3
```

## Citation
@article{mizumoto2025, title={Observation of tandem running behavior in dealates of Asian dampwood termite, Hodotermopsis sjostedti}, author={Mizumoto, Nobuaki and Chambliss, William and Carroll, Elijah P and Nakazono, Tomohiro and Kanao, Taisuke}, journal={Journal of Ethology}, year={2025}, doi={10.1007/s10164-025-00857-5} }

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact
William Chambliss: wlc0018@auburn.edu  
Nobuaki Mizumoto: nzm0095@auburn.edu
