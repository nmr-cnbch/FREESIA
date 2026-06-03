# FREESIA
A Python script for analysis of spectra obtained from FREESIA experiments.

## Running the script
To run the script first make sure the following packages are installed:
* [uv](https://docs.astral.sh/uv/getting-started/installation/)

Then open the terminal in the directory with the script files and run it using the following command in the terminal:
```
uv run freesia.py accordion_dim input_dir output_dir
```
For the list of optional arguments, use the command `uv run freesia.py --help` or `uv run freesia.py -h`, or see the [Optional arguments section](#optional-arguments).

Before running the script, `uv` checks if the appropriate virtual environment already exists. If not, it will create it and download suitable Python version and required packages.

### Required arguments
* `accordion_dim` - number of the indirect dimension whose evolution time was co-incremented with the irradiation frequency.
* `input_dir` - path of the directory which contains the input files.
* `output_dir` - path of the directory where the output files will be saved. The directory will be created if it does not exist. If it already exists, the files in the directory may be overwritten!

### Optional arguments
* `-h, --help` - show the help message and exit
* `-w, --weighting [WEIGHTING ...]` - weighting function by which the cross sections get multiplied, written as `function[parameter]`, for example `c[4]`. Four types of functions are supported: c, s, e, g. They are defined in Section 3 of the Supporting Information of the publication (doi: ). Multiple functions can be used (e.g. `-w c[4] s[6]`), each function will be used on the data separately and more subfolders will be created with output for every used weighting function. (default: `g[5]`)
* `-spectrum, --spectrum_files [spectrum_file_name ...]` - name or path of the spectrum file. Use this argument if you wish to use a different file name or directory than default (or want to use the sum of multiple spectra). Spectrum used should have been created from data which was submitted to apodization in every dimension except the accordion dimension. Apodization in accordion dimension will lead to worse or incorrect results. Multiple spectra can be provided, they will be added together, but they all must have exactly the same parameters and numbers of points in each dimension - otherwise the data in the output might be incorrect! (default: ucsf)
* `-fws, --fully_weighted_spectrum_files [fully_weighted_spectrum_real_file_name ...]` - name or path of the fully weighted spectrum file. Use this argument if you wish to use a different file name or directory than default (or want to use the sum of multiple spectra). This spectrum is used to determine overlap within cross sections of the spectrum. If overlap is detected, the neighbouring cross sections (along all dimensions except the accordion dimension) are checked and the one with the lowest level of overlap is selected for profile calculation. This spectrum must have been created from data which was submitted to apodization in every dimension, including the accordion dimension (where the apodization must remove the sinc disturbance of the lineshape - the recommended apodization function is squared cosine). Multiple spectra can be provided, they will be added together, but they all must have exactly the same parameters and numbers of points in each dimension. (default: weighted_spectrum.ucsf)
* `-plf, --peak_list_file peak_list_file_name` - name or path of the peak list file. Use this argument if you wish to use a different file name or directory than default. (default: peak.list)
* ` -flf, --freq_list_file frequency_list_file_name` - name or path of the frequency list file. Use this argument if you wish to use a different file name or directory than default. (default: fq1list)
* `-tol, --tolerance TOLERANCE` - acceptable tolerance (chemical shift difference in ppm) between the detected chemical shift and the reference chemical shift. Requires the `--reference` argument with the appropriate file. (default: 0.25)
<!--* `--latex` - create a text file that will include a latex code for a table containing the main results of the script: peak names, their positions in the spectrum (chemical shifts), chemical shifts determined from the amplitude profiles for the accordion dimension, reference chemical shifts (if -ref argument was provided) and the assigned profile group (default: False)-->
* `--figures {individual,combined,both}` - available options: individual - a plot of every cross section and profile is saved into a separate png file; combined - plots for profiles are combined into fewer figures, one figure having 35 profile plots; both - both individual and combined options are used (default: None)
* `--text {cs,profiles,both}` - available options: cs - intensities of consecutive  cross section points are saved as text files; profiles - intensities of consecutive profile points are saved as text files; both - values for both cross sections and profiles are saved as text files (default: None)
* `-ref, --reference_file reference_file_name` - name or path of the file with reference chemical shifts which will be used to compare the detected chemical shifts. Detected shifts are interpreted as correct if they match the reference values within the tolerance passed in the `--tolerance` argument (0.25 ppm by default). It should be a two-column text file with the name of the peak in the first column and its reference chemical shift in the second (columns separated by white spaces). (default: None)
* `-v, --verbosity` - Increase verbosity level for the logging (default: 0)

### Input files
* Text file with a list of irradiation frequencies used in the experiment. The format is identical as the format required by TopSpin program for Bruker spectrometers - the first line should contain `bf ppm` and then every subsequent line should contain a floating point value of an irradiation frequency in ppm. By default, the name of the file should be `fq1list`, otherwise it must be specified in the `-flf, --freq_list_file` argument.
* Text file with a list of peaks in the spectrum, in the [Sparky](https://nmrfam.wisc.edu/nmrfam-sparky-distribution/) format. By default, the name of the file should be `peak.list`, otherwise it must be specified in the `-plf, --peak_list_file` argument.
* Spectrum obtained by Fourier transformation (in each evolved dimension) of data measured in a FREESIA experiment, in ucsf format. Spectrum used should have been created from data which was submitted to apodization in every dimension except the accordion dimension. By default, the name of the file should be `ucsf`, otherwise it must be specified in the `-spectrum, --spectrum_files` argument.
* Spectrum measured in a FREESIA experiment, in ucsf format, but the data should have been submitted to apodization in every dimension including the accordion dimension. By default, the name of the file should be `weighted_spectrum.ucsf`, otherwise it must be specified in the `-fws, --fully_weighted_spectrum_files` argument.
* Optional: text file with reference chemical shifts for each peak. It must be specified in the `-ref, --reference_file` argument. Every line of the file should contain the name of the peak (the same as in the list of peaks in the spectrum), at least one white space and then the value of the reference chemical shift (in ppm).

### Output
The output directory contains several files:
* `output_table.csv` with peak names, chemical shifts in every spectral dimension, reference chemical shifts (if provided via the `-ref` argument), chemical shifts determined from obtained amplitude profiles for every used weighting function, and validity of obtained profiles.
* `FREESIA_output_summary.txt` contains numbers and percentage of valid results
* Log file with different level of information, according to the use of `--verbosity` argument. If the log file already exists in the output directory, new logs are appended to it.
<!--* Optional: text file `latex_table.txt` with output data in a table written in LaTeX. Created if `--latex` option is used.-->

The output directory contains one subdirectory for cross-sections before any processing and one subdirectory with output for each weighting function used.
* Optional: a set of .png files with the plots of cross-sections and profiles, either as one file per cross-section/profile, or with fewer files with multiple plots in each image. Depends on the use of the `--figures` argument.
* Optional: text files with intensities of each point in the cross-sections and/or profiles, according to the use of `--text` argument.

## Test data
Directory `MBP_input` contains the input files for the script from the H<ins>N</ins>CA(COirr) experiment performed on MBP sample, as reported in the original publication.<!--(add reference to the published paper here)--> Due to GitHub file size limitations, the required ucsf spectra have to be downloaded from the repository [Dane Badawcze UW](https://doi.org/10.58132/5X0W5D) (directory: `MBP_HNCA_COirr_FREESIA/ucsf_spectra`) and added to the `MBP_input`. Once all files are in place, pass `MBP_input` as the input directory for the script. Use `2` as the number of accordion dimension (corresponds to 15N dimension). You can use the following command after replacing `input_dir` and `output_dir` with paths to the input and output directories on your computer: `uv run freesia.py 2 input_dir output_dir -spectrum ucsf_not_weighted_along_Ndim -fws ucsf_weighted_along_Ndim`.

## Citation
Please cite: ...

## Problems or questions?
If you encounter issues or have questions, we recommend writing an email: anzaw@chem.uw.edu.pl

<!--\## References
[1\] Insert our article (and supporting information? here --> 
