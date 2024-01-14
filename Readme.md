A 4D Flow Processing Module to Compute the Pulsatility Index Transmission Coefficient (PITC) (Based on [Quantitative Velocity Tool (QVT)](https://github.com/uwmri/QVT))
=========
## Introduction ##
This library was developed to measure Pulsatility Transmission (PITC), some extra functionality has been added to extract time continuous boundary conditions which can be used to drive flow in cerebral hemodynamic models as well. 
An example extrapolation of PITC is shown below from the submitted manuscript describing the biomarker.

![TestSizeCoarse400](https://github.com/ABI-Animus-Laboratory/QVTplus/assets/108192400/0f48fa5a-f434-4e9f-852f-3f3d02526581)


Major changes from QVT include the development of an improved flow  quality metric (Q) which allows the automatic extraction of the highest quality flow from each cerebral vessel. This was necessary as there are several "good" flow positions which increase the interobserver error when using this flow profile as input into patient specific modelling or when weighting the measurement of PITC.

## Installation ##
Requires MATLAB version > 2018
### Dependencies ###
**Required Matlab Add-Ons** \
Image Processing Toolbox (for medfilt3) \
Curve Fitting Toolbox (for csaps) \
Statistics and Machine Learning Toolbox (for kmeans)
### Compatibility ###
Currently, .hdf5 and .dat formatting is not supported with the PITC measurement and only loading DICOMS allow this functionality. Resolutions can now be anisotropic.
## Usage
Run paramMap as per QVT, once the file has been processed, fill in the table with appropriate branch labels (P1) and sampling locations if you intend to extract local flow as well (P5). If the artery is disconnected into several segments (common in ICA) the label will change, so the connection matrix should be supplied (P2). Simply comma between numbers 67,112,etc. You must also include labels to remove posterior communicants (P3), and you can specify left and right ACA's as well in (P4). Once the small major artery Label section is full, click Done. 

![Figure2_Interaction](https://github.com/ABI-Animus-Laboratory/QVTplus/assets/108192400/c52a5d6a-3cee-49e1-9559-9ca4bcfbbbd4)

Run the separate function "UPDATE" which will ask you to point to where the subject data is saved. It will take the label table and measure PITC using Q weights by first connecting the branches and then assigning pulsatility distance (see manuscript for details). If flow locations were also specified, run flow_pipeline for flow and error extraction. The average flow (based on quality cutoff) is used for boundary conditions. Plots of flow and error are also included.

If at any point you want to manually save your own best chosen flow, follow standard QVT instruction and run flow_pipeline_manual, filling in necessary paths at top of script.

### Acknowledging ### 
Standing on the shoulders of giants, this module would not be possible without the original creators work. Please properly cite and acknowledge the following publications:

- [Roberts GS, Hoffman CA, Rivera-Rivera LA, Berman SE, Eisenmenger LB, Wieben O. Automated hemodynamic assessment for cranial 4D flow MRI. Magn Reson Imaging. 2022 Dec 26:S0730-725X(22)00231-4. doi: 10.1016/j.mri.2022.12.016. Epub ahead of print. PMID: 36581214](https://pubmed.ncbi.nlm.nih.gov/36581214/)
- [Schrauben E, Wahlin A, Ambarki K, Spaak E, Malm J, Wieben O, Eklund A. (2015). Fast 4D flow MRI intracranial segmentation and quantification in tortuous arteries. J Magn Reson Imaging, 42(5), 1458-1464. doi:10.1002/jmri.24900](https://pubmed.ncbi.nlm.nih.gov/25847621/)

And our own
- (Fill when accepted)

### Contact ### 
Please contact Sergio Dempsey if you have a dataset that is not loading correctly and that functionality will be added within a reasonable timeframe.

**Sergio Dempsey**: sdem348@aucklanduni.ac.nz
