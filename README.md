# absolute-abundance-16s

## Absolute quantification of prokaryotes in the microbiome by 16S rRNA qPCR or ddPCR

## citations

## updates

# Table of contents
1. [Description](#description)
2. [Using the scripts and notebooks](#using-the-scripts-and-notebooks)
      * [Setup](#setup)
      * [Examples](#examples)
          * [qPCR](#qpcr)
          * [ddPCR](#ddpcr)
          * [Notebooks](#notebooks)
      * [Prepare input files](#prepare-input-files)
      * [Run scripts and view output](#run-scripts-and-view-output)
      * [User adjustable parameters](#user-adjustable-parameters)

# Description
This repo contains the following directories:
* source code in ```src```
* scripts in ```scripts```
    * qPCR-specific analysis, following the steps described in the protocol to perform quality control and calculate 16S rRNA copies per reaction from qPCR data
    * ddPCR-specific analysis, following the steps described in the protocol to perform quality control and calculate 16S rRNA copies per reaction from ddPCR data
    * universal analysis, following the steps described in the protocol to assess controls and calculate 16S rRNA copies per dry gram of input stool
* key general files that the scripts rely on in ```artifacts```
* examples in ```examples```
    * qPCR example, including toy data input files, bash files to run the scripts, and toy data example output
    * ddPCR example, including toy data input files, bash files to run the scripts, and toy data example output
* data and code specific to manuscripts in ```manuscripts```

# Using the scripts and notebooks
## Setup

**Clone the repository**

	git clone https://github.com/bhattlab/absolute-abundance-16s.git

**Install dependencies**

Install conda (https://developers.google.com/earth-engine/guides/python_install-conda/), if not already installed. This will provide Python (tested version 3.12.1).

Then execute the following commands.

    pip install click==8.1.7
    pip install matplotlib==3.8.3
    pip install pandas==2.2.1
    pip install scipy==1.12.0
    pip install seaborn==0.13.2
    pip install openpyxl==3.1.2
    pip install jupyterlab

Note: These are the tested package versions, but the scripts most likely work with other versions as well.

## Examples
### qPCR

Navigate to

    cd examples/qpcr_v1

Perform qPCR-specific analysis

    bash qpcr_specific_analysis.sh

And then use the output from qPCR-specific analysis for universal analysis

    bash universal_analysis_from_qpcr.sh

As each analysis script proceeds, it will output key information to the command line. The output files will appear in ```output``` and can be compared to the ```output_expected```.

### ddPCR

Navigate to

    cd examples/ddpcr_v1

Perform ddPCR-specific analysis

    bash ddpcr_specific_analysis.sh

And then use the output from ddPCR-specific analysis for universal analysis

    bash universal_analysis_from_ddpcr.sh

As each analysis script proceeds, it will output key information to the command line. The output files will appear in ```output``` and can be compared to the ```output_expected```.

### Notebooks

We highly recommend analyzing data interactively. One way to do this is with a notebook in Jupyter Lab.

To see a notebook version of either the qPCR or ddPCR example analysis, first navigate to ```absolute-abundance-16s```.

Then launch Jupyter Lab:

    jupyter lab

And navigate to ```http://localhost:8888/lab``` in your browser window.

Open the file ```examples/qpcr_v1/qpcr_notebook_template.ipynb``` for qPCR (including universal analysis) or ```examples/ddpcr_v1/ddpcr_notebook_template.ipynb``` for ddPCR (including universal analysis).

These are setup to view the example data, but you may duplicate the notebooks and edit filepaths for your own data. The input files (e.g. [Prepare input files](#prepare-input-files)) and parameters (e.g. [User adjustable parameters](#user-adjustable-parameters)) are the same as in the scripts.
The notebooks do not automatically save any outputs, but they can be edited to do so.

## Prepare input files
### qPCR-specific analysis
**1. qPCR data Excel spreadsheet with two columns**
<table border="0">
 <tr>
    <td><b style="font-size:18px">column name</b></td>
    <td><b style="font-size:18px">description</b></td>
 </tr>
 <tr>
    <td>Well</td>
    <td>Well number on the 384-well plate, e.g. A1 through P24</td>
 </tr>
 <tr>
    <td>Cq</td>
    <td>Numerical value or Undetermined if a given well was empty or was a failed technical replicate</td>
 </tr>
</table>

**2. qPCR layout Excel spreadsheet with seven columns**
<table border="0">
 <tr>
    <td><b style="font-size:18px">column name</b></td>
    <td><b style="font-size:18px">description</b></td>
 </tr>
 <tr>
    <td>Well96</td>
    <td>Well number on the 96-well plate, e.g. A1 through H12</td>
 </tr>
 <tr>
    <td>Name</td>
    <td>Must be unique for each sample, control, and standard dilution. Should be NIST_A for component A of NIST etc. and NIST_mix_A_R for the NIST mixture from Reagent Setup. Other specific names do not matter.</td>
 </tr>
 <tr>
    <td>LiquidHandlerDilution</td>
    <td>Numerical value that corresponds to the liquid handler or 96-well format dilution, e.g. 1000 for 1:1000 dilution. Should be 1 if no dilution is performed.</td>
 </tr>
 <tr>
    <td>SinglePipettorDilution</td>
    <td>Numerical value that corresponds to the single pipettor dilution, e.g. 1000 for 1:1000 dilution. Should be 1 if no dilution is performed.</td>
 </tr>
 <tr>
    <td>Type</td>
    <td>Options in qPCR are "PCRPos" (e.g. for NIST), "PCRNeg" (e.g. for no template control), "DNAPos" (e.g. for Zymo mock), "DNANeg" (e.g. for extraction from water), "Pvul" (e.g. for P. vulgatus standard curve points), "Fpra" (e.g. for F. prausnitzii standard curve points), and "Sample" (e.g. for all samples).</td>
 </tr>
 <tr>
    <td>uLAdded</td>
    <td>The number of diluted uL of sample, control, or standard added to the reaction. This value is 6 uL in our protocol.</td>
 </tr>
 <tr>
    <td>ElutionVolume</td>
    <td>The elution volume (uL) in the final step of DNA extraction. Must be provided for Type=DNANeg, DNAPos, or Sample. This value is 100 uL in our protocol.</td>
 </tr>
</table>
  
### ddPCR-specific analysis
**1. ddPCR data Excel spreadsheet with four columns**
<table border="0">
 <tr>
    <td><b style="font-size:18px">column name</b></td>
    <td><b style="font-size:18px">description</b></td>
 </tr>
 <tr>
    <td>Well</td>
    <td>Well number on the 96-well plate, e.g. A01 through H12</td>
 </tr>
 <tr>
    <td>Accepted Droplets</td>
    <td>The number of accepted droplets, e.g. the sum of positive droplets and negative droplets</td>
 </tr>
 <tr>
    <td>Positives</td>
    <td>The number of positive droplets, e.g. droplets above the threshold set in QX Manager</td>
 </tr>
 <tr>
    <td>Negatives</td>
    <td>The number of negative droplets, e.g. droplets below the threshold set in QX Manager</td>
 </tr>
</table>

**2. ddPCR layout Excel spreadsheet with seven columns**
<table border="0">
 <tr>
    <td><b style="font-size:18px">column name</b></td>
    <td><b style="font-size:18px">description</b></td>
 </tr>
 <tr>
    <td>Well96</td>
    <td>Well number on the 96-well plate, e.g. A1 through H12</td>
 </tr>
 <tr>
    <td>Name</td>
    <td>Must be unique for each sample and control. Should be NIST_A for component A of NIST etc. and NIST_mix_A_R for the NIST mixture from Reagent Setup. Other specific names do not matter.</td>
 </tr>
 <tr>
    <td>LiquidHandlerDilution</td>
    <td>Numerical value that corresponds to the liquid handler or 96-well format dilution, e.g. 1000 for 1:1000 dilution. Should be 1 if no dilution is performed.</td>
 </tr>
 <tr>
    <td>SinglePipettorDilution</td>
    <td>Numerical value that corresponds to the single pipettor dilution, e.g. 1000 for 1:1000 dilution. Should be 1 if no dilution is performed.</td>
 </tr>
 <tr>
    <td>Type</td>
    <td>Options in ddPCR are "PCRPos" (e.g. for NIST), "PCRNeg" (e.g. for no template control), "DNAPos" (e.g. for Zymo mock), "DNANeg" (e.g. for extraction from water), and "Sample" (e.g. for all samples).</td>
 </tr>
 <tr>
    <td>uLAdded</td>
    <td>The number of diluted uL of sample, control, or standard added to the reaction. This value is 6 uL in our protocol.</td>
 </tr>
 <tr>
    <td>ElutionVolume</td>
    <td>The elution volume (uL) in the final step of DNA extraction. Must be provided for Type=DNANeg, DNAPos, or Sample. This value is 100 uL in our protocol.</td>
 </tr>
</table>

### Universal analysis
**1. The output from qPCR-specific or ddPCR-specific analysis, specifically the "for_universal_analysis" sheet**

**2. Weights Excel spreadsheet with six columns**
<table border="0">
 <tr>
    <td><b style="font-size:18px">column name</b></td>
    <td><b style="font-size:18px">description</b></td>
 </tr>
 <tr>
    <td>Name</td>
    <td>Must correspond to the Name from the qPCR or ddPCR layout file</td>
 </tr>
 <tr>
    <td>empty_wt</td>
    <td>The mass (g) of the empty tube used for drying to measure moisture content</td>
 </tr>
 <tr>
    <td>filled_wt</td>
    <td>The mass (g) of the tube used for drying with the wet stool for drying in it</td>
 </tr>
 <tr>
    <td>dry_wt</td>
    <td>The mass (g) of the tube used for drying with the now dry stool from drying in it</td>
 </tr>
 <tr>
    <td>empty_PB</td>
    <td>The mass (g) of the empty PowerBead tube for DNA extraction</td>
 </tr>
 <tr>
    <td>filled_PB</td>
    <td>The mass (g) of the PowerBead tube with the wet stool for DNA extraction in it</td>
 </tr>
</table>

## Run scripts and view output

Create a new directory for the prepared input files from the previous step and create another directory for the script output.

### qPCR

Edit and execute the following command for qPCR-specific analysis to perform quality control and calculate 16S rRNA copies per reaction.

    python full_path_to/scripts/qpcr_specific_analysis.py -q full_path_to/qpcr_data.xlsx -qs Sheet1 -l full_path_to/qpcr_layout.xlsx -ls Sheet1 -f full_path_to/artifacts/format_conversion.tsv -o outputfolder --pcopy insertpcopy --fcopy insertfcopy

The script will output key information to the command line and files to the specified output directory. The qPCR-specific analysis script saves plots of the standard curve after each step:

<table border="0">
 <tr>
    <td><b style="font-size:18px">filename</b></td>
    <td><b style="font-size:18px">description</b></td>
 </tr>
 <tr>
    <td>standard_curve_visualization_post_step_90.pdf</td>
    <td>Initial visualization</td>
 </tr>
 <tr>
    <td>standard_curve_visualization_post_step_91.pdf</td>
    <td>After removing standard curve dilution points with a large technical replicate span</td>
 </tr>
 <tr>
    <td>standard_curve_visualization_post_step_92.pdf</td>
    <td>After removing concentrated dilution points that plateau</td>
 </tr>
 <tr>
    <td>standard_curve_visualization_post_step_93.pdf</td>
    <td>After removing dilute dilution points too near the limit of blank</td>
 </tr>
 <tr>
    <td>standard_curve_visualization_final_post_step_95.pdf</td>
    <td>The resulting dilution points and final regression line</td>
 </tr>
</table>

At the end, the qPCR-specific analysis script outputs an Excel spreadsheet ```qpcr_specific_output.xlsx``` with the following sheets:

<table border="0">
 <tr>
    <td><b style="font-size:18px">sheet name</b></td>
    <td><b style="font-size:18px">description</b></td>
 </tr>
 <tr>
    <td>for_universal_analysis</td>
    <td>All samples and controls that passed quality control</td>
 </tr>
 <tr>
    <td>removed_wells_tech_reps_fails</td>
    <td>All wells removed because they do not have at least two successful technical replicates</td>
 </tr>
 <tr>
    <td>removed_standards</td>
    <td>Standard curve dilution points removed due to a large technical replicate span, being too concentrated, or being too dilute</td>
 </tr>
 <tr>
    <td>removed_high_variation_samples</td>
    <td>Samples and controls removed due to technical replicate variation</td>
 </tr>
 <tr>
    <td>removed_samples</td>
    <td>Samples and controls removed due to being too concentrated, too dilute, or low confidence but not undiluted</td>
 </tr>
</table>

After reviewing these, edit and execute the following command for universal analysis to assess controls and calculate 16S rRNA copies per dry gram of input stool.

    python full_path_to/scripts/universal_analysis.py -d outputfolder/qpcr_specific_output.xlsx -ds for_universal_analysis -w full_path_to/weights.xlsx -ws Sheet1 -n full_path_to/artifacts/nist_expected_values_03262024.xlsx -ns Sheet1 -o outputfolder

The universal analysis script outputs an Excel spreadsheet ```univ_analysis_output.xlsx``` with the following sheets:

<table border="0">
 <tr>
    <td><b style="font-size:18px">sheet name</b></td>
    <td><b style="font-size:18px">description</b></td>
 </tr>
 <tr>
    <td>from_universal_analysis</td>
    <td>All samples</td>
 </tr>
 <tr>
    <td>nist_measured_expected</td>
    <td>NIST positive PCR controls measured to expected</td>
 </tr>
 <tr>
    <td>neg_dna_extract_controls</td>
    <td>Negative DNA extraction controls</td>
 </tr>
 <tr>
    <td>pos_dna_extract_controls</td>
    <td>Zymo mock positive DNA extraction controls</td>
 </tr>
 <tr>
    <td>extract_input_outside_range</td>
    <td>List of samples with DNA extraction input outside of desired range.</td>
 </tr>
 <tr>
    <td>drying_input_outside_range</td>
    <td>List of samples with drying aliquot (e.g. for stool moisture content) outside of desired range.</td>
 </tr>
 <tr>
    <td>dry_stool_amount_low</td>
    <td>List of samples with a small mass of dried stool after drying, which increases the error</td>
 </tr>
 <tr>
    <td>water_fraction_over_cutoff</td>
    <td>List of samples with a water fraction over the cutoff, which similarly indicates increased error</td>
 </tr>
</table>

Note: unlike in qPCR- and ddPCR-specific analysis, the samples identified in the latter four sheets of universal analysis are not removed. It is up to the user to assess the situation on a case-by-case basis.

### ddPCR

Edit and execute the following command for ddPCR-specific analysis to perform quality control and calculate 16S rRNA copies per reaction.

    python full_path_to/scripts/ddpcr_specific_analysis.py -d full_path_to/ddpcr_data.xlsx -ds Sheet1 -l full_path_to/ddpcr_layout.xlsx -ls Sheet1 -o outputfolder

The script will output key information to the command line and files to the specified output directory. At the end, the ddPCR-specific analysis script outputs an Excel spreadsheet ```ddpcr_specific_output.xlsx``` with the following sheets:

<table border="0">
 <tr>
    <td><b style="font-size:18px">sheet name</b></td>
    <td><b style="font-size:18px">description</b></td>
 </tr>
 <tr>
    <td>for_universal_analysis</td>
    <td>All samples and controls that passed quality control</td>
 </tr>
 <tr>
    <td>too_few_droplets</td>
    <td>Wells with a low number of droplets</td>
 </tr>
 <tr>
    <td>too_concentrated</td>
    <td>Samples and controls removed due to being too concentrated (e.g. need more dilution)</td>
 </tr>
 <tr>
    <td>too_dilute</td>
    <td>Samples and controls removed due to being too dilute (e.g. need less dilution)</td>
 </tr>
</table>

After reviewing these, edit and execute the following command for universal analysis to assess controls and calculate 16S rRNA copies per dry gram of input stool.

    python full_path_to/scripts/universal_analysis.py -d outputfolder/ddpcr_specific_output.xlsx -ds for_universal_analysis -w full_path_to/weights.xlsx -ws Sheet1 -n full_path_to/artifacts/nist_expected_values_03262024.xlsx -ns Sheet1 -o outputfolder

See the qPCR section for information on the universal analysis script output.

## User adjustable parameters

The scripts provide default values for all parameters, which align with the recommendations provided in the protocol. However, if the user would like to modify parameters in a particular case, a full list of parameters that can be passed on the command line and their descriptions can be found in the help text of each script. Use the following commands to see the help text (also displayed below) in the terminal window.

    python scripts/qpcr_specific_analysis.py --help
    python scripts/ddpcr_specific_analysis.py --help
    python scripts/universal_analysis.py --help

### qPCR-specific analysis

```
Usage: qpcr_specific_analysis.py [OPTIONS]

Options:
  -q, --qpcr-path                 qPCR data path to Excel file.  [required]
  -qs, --qpcr-sheet               qPCR data Excel sheet name.  [required]
  -l, --layout-path               96-well layout path to Excel file.  [required]
  -ls, --layout-sheet             96-well layout Excel sheet name.  [required]
  -f, --format-conversion-path    format conversion file path.  [required]
  -o, --output-path               output folder path.  [required]
  --pcopy                         16S rRNA copies in uLAdded (e.g. 6 uL) of stock P. vulgatus standard plasmid.
                                  [100000000<=x<=1000000000000; required]
  --fcopy                         16S rRNA copies in uLAdded (e.g. 6 uL) of stock F. prausnitzii standard plasmid.
                                  [100000000<=x<=1000000000000; required]
  --num-of-tech-reps              number of qPCR technical replicates.  [default: 3; 2<=x<=3]
  --max-cq-span-ntc               maximum span of median Cq values of all no template controls.  [default: 2;
                                  0.1<=x<=3.3]
  --max-cq-span-standard-dil-pt   maximum span of Cq values of each standard dilution point.  [default: 2;
                                  0.1<=x<=3.3]
  --max-standard-dil-pts-removed-tech-var
                                  maximum number of standard dilution points removed for technical replicate
                                  variation.  [default: 1; 0<=x<=4]
  --min-cq-gap-conc-standards     minimum Cq gap between concentrated points of the standard curve (e.g. no plateau).
                                  [default: 3.11; 2.81<=x<=3.41]
  --cq-cutoff-conc-standards      maximum Cq value to consider gap between concentrated points of the standard curve.
                                  [default: 15; 8<=x<=18]
  --cq-standards-sep-lob          minimum Cq separation between most dilute standard dilution point and limit of
                                  blank.  [default: 2; 1<=x<=6.6]
  --max-fold-change-pvul-fpra     maximum fold change between P. vulgatus and F. prausnitzii standard curves.
                                  [default: 2; 0.1<=x<=4]
  --most-steep-slope-allowed      most steep slope allowed for the standard curve.  [default: -3.58; -3.98<=x<=-3.31]
  --least-steep-slope-allowed     least steep slope allowed for the standard curve.  [default: -3.11; -3.29<=x<=-2.71]
  --min-r-squared                 minimum R squared value for the standard curve.  [default: 0.98; 0.93<=x<=0.999999]
  --cq-non-standards-sep-lob      minimum Cq separation between samples or controls and limit of blank.  [default: 2;
                                  1<=x<=6.6]
  --overhang-allowed              determines whether samples are allowed to overhang the dilute end of the standard
                                  curve, while remaining the minimum distance from the limit of blank.  [default:
                                  False]
  --max-copies-rxn-lob            maximum acceptable apparent 16S rRNA copies per reaction of the no template
                                  controls.  [default: 500; 5<=x<=3000]
  --max-cq-diff-sample-closest-two-reps
                                  maximum difference between the Cq values of the two closest technical replicates for
                                  samples or controls.  [default: 2; 0.1<=x<=3.3]
  --cq-low-conf-sep-lob           Cq separation between samples or controls and limit of blank to decide a measurement
                                  is low confidence.  [default: 3.3; 1<=x<=6.6]
  --help                          Show this message and exit.
```

### ddPCR-specific analysis

```
Usage: ddpcr_specific_analysis.py [OPTIONS]

Options:
  -d, --ddpcr-path            ddPCR data path to Excel file.  [required]
  -ds, --ddpcr-sheet          ddPCR data Excel sheet name.  [required]
  -l, --layout-path           96-well layout path to Excel file.  [required]
  -ls, --layout-sheet         96-well layout Excel sheet name.  [required]
  -o, --output-path           output folder path.  [required]
  --droplet-volume            volume (nL) of each droplet.  [default: 0.795; 0.75<=x<=0.92]
  --rxn-volume                volume (uL) of ddPCR reaction as setup on the pre-droplet plate.  [default: 22;
                              20<=x<=30]
  --min-accepted-droplets     minimum number of accepted droplets per well.  [default: 10000; 5000<=x<=20000]
  --max-copies-rxn-lob        maximum 16S rRNA copies per reaction of the no template controls.  [default: 25;
                              0<=x<=100]
  --max-copies-rxn-span-ntc   maximum span (e.g. fold change) of 16S rRNA copies per reaction of all no template
                              controls.  [default: 4; 1<=x<=10]
  --min-negative-droplets     minimum number of negative droplets per well.  [default: 10; 0<=x<=100]
  --copies-rxn-loq-mult       multiplied by the limit of blank to define the limit of quantification, under which
                              samples or controls are removed.  [default: 4; 2<=x<=100]
  --help                      Show this message and exit.
```

### Universal analysis

```
Usage: universal_analysis.py [OPTIONS]

Options:
  -d, --data-path                 qPCR- or ddPCR-specific analysis output Excel file path.  [required]
  -ds, --data-sheet               qPCR- or ddPCR-specific analysis output Excel sheet name (e.g.
                                  'for_universal_analysis').  [required]
  -w, --weights-path              sample weights Excel file path.  [required]
  -ws, --weights-sheet            sample weights Excel sheet name.  [required]
  -n, --nist-expected-path        NIST expected values Excel file path. Must contain 16S rRNA copies per undiluted uL
                                  in a column titled 'copies_uL_expected' and a name that matches the name from the
                                  layout for qPCR- or ddPCR-specific analysis in a column titled 'Name'.  [required]
  -ns, --nist-expected-sheet      NIST expected values Excel sheet name.  [required]
  -o, --output-path               output folder path.  [required]
  --nist-max-fold-diff            maximum fold difference between measured and expected NIST control 16S rRNA copies
                                  per undiluted uL.  [default: 5.0; 1.01<=x<=10]
  --neg-extract-ctrl-max-copies   maximum 16S rRNA copies per DNA extraction for negative DNA extraction controls.
                                  [default: 5000.0; 500<=x<=30000]
  --pos-extract-ctrl-max-span     maximum span (e.g. fold change) of 16S rRNA copies per DNA extraction for positive
                                  DNA extraction controls.  [default: 2; 1.1<=x<=10]
  --extract-max-input             maximum amount of stool (g) from which to extract DNA.  [default: 0.25;
                                  0.025<=x<=0.5]
  --extract-min-input             minimum amount of stool (g) from which to extract DNA.  [default: 0.15;
                                  0.025<=x<=0.5]
  --drying-max-input              maximum amount of stool (g) to dry for stool moisture content.  [default: 0.125;
                                  0.025<=x<=0.5]
  --drying-min-input              minimum amount of stool (g) to dry for stool moisture content.  [default: 0.075;
                                  0.025<=x<=0.5]
  --min-dried-dry-mass            minimum dried amount of stool (g) from drying for stool moisture content.  [default:
                                  0.008; 0.002<=x<=0.05]
  --water-fraction-cutoff         cutoff for water fraction, given the error in 16S rRNA copies per dry gram as the
                                  water fraction increases.  [default: 0.9; 0.8<=x<=0.99]
  --help                          Show this message and exit.
```
