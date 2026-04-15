# CIRSense
CIRSense utilizes Channel State Information (CSI) data from Wi-Fi signals for non-contact sensing applications, including distance estimation and respiration rate monitoring. 

This project implements the algorithms described in the paper:  **CIRSense: Rethinking WiFi Sensing with Channel Impulse Response**

### BibTeX Citation

```bibtex
@ARTICLE{11457720,
  author={Kong, Ruiqi and Chen, He},
  journal={IEEE Transactions on Mobile Computing}, 
  title={CIRSense: Rethinking WiFi Sensing With Channel Impulse Response}, 
  year={2026},
  volume={},
  number={},
  pages={1-15},
  doi={10.1109/TMC.2026.3678326}}
```

## Files Overview

### Core Scripts

- **`cirsense_distance.m`**: Single-target distance estimation script. Processes CSI data to estimate the distance to a single target and compares it with ground truth.
  
- **`cirsense_bpm.m`**: Single-target respiration rate monitoring script. Estimates breathing rate from CSI data and compares with ground truth respiration data.

- **`cirsense_distance_multitarget.m`**: Multi-target distance estimation script. 

- **`cirsense_bpm_multitarget.m`**: Multi-target respiration rate monitoring script.
  
### Supporting Functions

- **`build_default_params`**: Defines default parameters for the algorithms.
- **`build_default_config`**: Sets default configuration paths and settings.
- **`get_window_length`**: Calculates smoothing window lengths.
- **`cirsense`**: Core CIRSense algorithm for distance estimation and respiration trace extraction.
- **`initialization`**: Initializes matrices and indices for signal processing.
- **`Domino`**: Preprocessing algorithm to correct hardware distortions (e.g., STO, CFO, ...) in CSI data.
- **`estimate_tau1_freq`**: Estimates time delays in the frequency domain.

### Data
- The input data can be downloaded from https://ieee-dataport.org/documents/cirsense-dataset-real-world-80211ax-csi-measurements-wireless-sensing

## Contact

ruiqikong@cuhk.edu.hk

## License

This project is for research purposes. Please cite appropriately if used in publications.
