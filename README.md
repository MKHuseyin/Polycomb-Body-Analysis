# Polycomb-Body-Analysis
Macros and R script for segmenting and analysing Polycomb bodies from spinning disk microscopy z-stacks of live cells expressing endogenous HaloTag Polycomb fusion proteins

ImageJ macros are used to process files once z-stacks are deconvolved and nuclei segmented, according to methods described in Huseyin and Klose 2020 (https://doi.org/10.1101/2020.04.25.061358). Both are designed to run with specific directory structures (described in the respective files), and will need to be modified to run correctly otherwise.

Similarly, R script is to be run on the resulting directory following both ImageJ macros and addition of tables of measurements of background pixel intensity for each z-stack in txt format. The included script is a generic example, including the extras for analysing multiple conditions.
