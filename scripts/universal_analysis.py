import os
import sys
import click
import matplotlib
import pandas as pd

matplotlib.use("Agg")

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../src')))
import absolute_abundance_functions as aa


@click.command(context_settings=dict(max_content_width=120, show_default=True))
@click.option(
    "-d",
    "--data-path",
    "data_path",
    required=True,
    type=click.Path(exists=True),
    help="qPCR- or ddPCR-specific analysis output Excel file path.",
    metavar="",
)
@click.option(
    "-ds",
    "--data-sheet",
    "data_sheet_name",
    required=True,
    type=str,
    help="qPCR- or ddPCR-specific analysis output Excel sheet name (e.g. 'for_universal_analysis').",
    metavar="",
)
@click.option(
    "-w",
    "--weights-path",
    "weights_path",
    required=True,
    type=click.Path(exists=True),
    help="sample weights Excel file path.",
    metavar="",
)
@click.option(
    "-ws",
    "--weights-sheet",
    "weights_sheet_name",
    required=True,
    type=str,
    help="sample weights Excel sheet name.",
    metavar="",
)
@click.option(
    "-n",
    "--nist-expected-path",
    "nist_expected_path",
    required=True,
    type=click.Path(exists=True),
    help="NIST expected values Excel file path. Must contain 16S rRNA copies per undiluted uL in a column titled 'copies_uL_expected' and a name that matches the name from the layout for qPCR- or ddPCR-specific analysis in a column titled 'Name'.",
    metavar="",
)
@click.option(
    "-ns",
    "--nist-expected-sheet",
    "nist_expected_sheet_name",
    required=True,
    type=str,
    help="NIST expected values Excel sheet name.",
    metavar="",
)
@click.option(
    "-o",
    "--output-path",
    "output_folder_path",
    required=True,
    type=click.Path(exists=True),
    help="output folder path.",
    metavar="",
)
@click.option(
    "--nist-max-fold-diff",
    "NIST_MAX_FOLD_DIFF",
    default=5.0,
    type=click.FloatRange(min=1.01, max=10),
    help="maximum fold difference between measured and expected NIST control 16S rRNA copies per undiluted uL.",
    metavar="",
)
@click.option(
    "--neg-extract-ctrl-max-copies",
    "NEG_EXTRACT_CTRL_MAX_COPIES",
    default=5000.0,
    type=click.FloatRange(min=500, max=30000),
    help="maximum 16S rRNA copies per DNA extraction for negative DNA extraction controls.",
    metavar="",
)
@click.option(
    "--pos-extract-ctrl-max-span",
    "POS_EXTRACT_CTRL_MAX_SPAN",
    default=2,
    type=click.FloatRange(min=1.1, max=10),
    help="maximum span (e.g. fold change) of 16S rRNA copies per DNA extraction for positive DNA extraction controls.",
    metavar="",
)
@click.option(
    "--extract-max-input",
    "EXTRACT_MAX_INPUT",
    default=0.250,
    type=click.FloatRange(min=0.025, max=0.500),
    help="maximum amount of stool (g) from which to extract DNA.",
    metavar="",
)
@click.option(
    "--extract-min-input",
    "EXTRACT_MIN_INPUT",
    default=0.150,
    type=click.FloatRange(min=0.025, max=0.500),
    help="minimum amount of stool (g) from which to extract DNA.",
    metavar="",
)
@click.option(
    "--drying-max-input",
    "DRYING_MAX_INPUT",
    default=0.125,
    type=click.FloatRange(min=0.025, max=0.500),
    help="maximum amount of stool (g) to dry for stool moisture content.",
    metavar="",
)
@click.option(
    "--drying-min-input",
    "DRYING_MIN_INPUT",
    default=0.075,
    type=click.FloatRange(min=0.025, max=0.500),
    help="minimum amount of stool (g) to dry for stool moisture content.",
    metavar="",
)
@click.option(
    "--min-dried-dry-mass",
    "MIN_DRIED_DRY_MASS",
    default=0.008,
    type=click.FloatRange(min=0.002, max=0.050),
    help="minimum dried amount of stool (g) from drying for stool moisture content.",
    metavar="",
)
@click.option(
    "--water-fraction-cutoff",
    "WATER_FRACTION_CUTOFF",
    default=0.9,
    type=click.FloatRange(min=0.8, max=0.99),
    help="cutoff for water fraction, given the error in 16S rRNA copies per dry gram as the water fraction increases.",
    metavar="",
)
def main(
    data_path,
    data_sheet_name,
    weights_path,
    weights_sheet_name,
    nist_expected_path,
    nist_expected_sheet_name,
    output_folder_path,
    **param_dict,
):

    print(f"\nThe files being loaded into the script are")
    print(f"data from qPCR- or ddPCR-specific analysis from {data_path} on sheet {data_sheet_name}")
    print(f"sample weights from {weights_path} on sheet {weights_sheet_name}")
    print(
        f"NIST (or other positive PCR control) expected values from {nist_expected_path} on sheet {nist_expected_sheet_name}"
    )
    print(f"output folder path is {output_folder_path}")

    print(f"\nThe parameters of the script are")
    for key, value in param_dict.items():
        click.echo(f"{key}: {value}")

    # Prepare matplotlib
    aa.prepare_matplotlib()

    # Load data
    data = pd.read_excel(data_path, sheet_name=data_sheet_name)
    weights = pd.read_excel(weights_path, sheet_name=weights_sheet_name)
    nist_expected = pd.read_excel(nist_expected_path, sheet_name=nist_expected_sheet_name)

    # Formatting
    nist_expected2 = aa.format_nist_expected(nist_expected)

    # Step 109
    print(f"\nStep 109")
    nist_measured = aa.calculate_nist_copies_uL(data[data["Type"] == "PCRPos"])
    print(f"16S rRNA copies per undiluted uL for NIST controls were calculated.")

    # Step 110
    print(f"\nStep 110")
    nist_measured_expected = aa.compare_nist_to_expected(nist_measured, nist_expected2, param_dict)
    print(
        f"{nist_measured_expected[nist_measured_expected['within_desired_range'].str.contains('no')].copy().shape[0]:.3g}"
        f" NIST controls were outside of the expected range."
    )

    # Step 111
    print(f"\nStep 111")
    print(f"This step is currently outside of the scope of this script but is an important step in many use cases.")
    print(f"Do not forget to assess data points when a given sample was assayed at multiple dilutions.")

    # Step 112
    print(f"\nStep 112")
    print(f"This step is currently outside of the scope of this script but is an important step in many use cases.")
    print(f"Do not forget to select a dilution when a given sample was assayed at multiple dilutions.")

    # Step 113
    print(f"\nStep 113")
    samples_controls1 = data[data["Type"].isin(["DNAPos", "DNANeg", "Sample"])]
    samples_controls2 = aa.calculate_copies_per_dna_extraction(samples_controls1)
    print(f"16S rRNA copies per DNA extraction were calculated.")

    # Step 114
    print(f"\nStep 114")
    negative_DNA_extraction_controls = samples_controls2[samples_controls2["Type"] == "DNANeg"]
    if negative_DNA_extraction_controls.shape[0] > 0:
        print(
            f"The maximum 16S rRNA copies per DNA extraction of the negative DNA extraction controls on this plate is "
            f"{negative_DNA_extraction_controls['copies_dna_extraction'].max():.3g}."
        )
        if (
            negative_DNA_extraction_controls["copies_dna_extraction"].max()
            > param_dict["NEG_EXTRACT_CTRL_MAX_COPIES"]
        ):
            print(
                f"This is greater than the maximum allowed value of {param_dict['NEG_EXTRACT_CTRL_MAX_COPIES']:.3g}."
            )
            raise RuntimeError("The script aborted. You may rerun it with different parameters, if appropriate.")
        else:
            print(f"This is less than the maximum allowed value of {param_dict['NEG_EXTRACT_CTRL_MAX_COPIES']:.3g}.")
    else:
        print(f"There were no negative DNA extraction controls above the limit of quantification on this plate.")

    # Step 115
    print(f"\nStep 115")
    print(f"This step is currently outside of the scope of this script but is an important step in many use cases.")
    print(
        f"Do not forget to remove samples if the 16S rRNA copies per DNA extraction is less than 4-fold above "
        f"that of the negative DNA extraction control by batch."
    )

    # Step 116
    print(f"\nStep 116")
    positive_DNA_extraction_controls, pos_DNA_extraction_fold_diff = aa.assess_positive_DNA_extraction_controls(
        samples_controls2
    )
    if positive_DNA_extraction_controls.shape[0] > 0:
        print(
            f"The ratio between the 16S rRNA copies per DNA extraction of the highest positive DNA extraction control and "
            f"of the lowest positive DNA extraction control on this plate is {pos_DNA_extraction_fold_diff:.3g}"
        )
        if pos_DNA_extraction_fold_diff > param_dict["POS_EXTRACT_CTRL_MAX_SPAN"]:
            print(
                f"This is greater than the maximum allowed value of {param_dict['POS_EXTRACT_CTRL_MAX_SPAN']:.3g}."
            )
            raise RuntimeError("The script aborted. You may rerun it with different parameters, if appropriate.")
        else:
            print(f"This is less than the maximum allowed value of {param_dict['POS_EXTRACT_CTRL_MAX_SPAN']:.3g}.")
    else:
        print(f"There were no positive DNA extraction controls on this plate.")
    print(f"Do not forget to compare positive DNA extraction controls across all DNA extraction batches.")

    # Step 117
    print(f"\nStep 117")
    weights2, samples_with_DNA_extraction_input_outside_range = aa.calculate_wet_mass_extracted_from(
        weights, param_dict
    )
    print(
        f"{samples_with_DNA_extraction_input_outside_range.shape[0]:.3g} samples were extracted from too much or too little input stool."
    )
    print(f"They are still included in downstream analysis.")

    # Step 118
    print(f"\nStep 118")
    weights3, samples_with_drying_input_outside_range, samples_with_dry_amount_of_stool_low = (
        aa.calculate_wet_and_dry_drying_mass(weights2, param_dict)
    )
    print(f"{samples_with_drying_input_outside_range.shape[0]:.3g} samples had too much or too little stool dried.")
    print(
        f"{samples_with_dry_amount_of_stool_low.shape[0]:.3g} samples have a small absolute dry weight of the drying aliquot, "
        f"which increases the error."
    )
    print(f"They are still included in downstream analysis.")

    # Step 119
    print(f"\nStep 119")
    weights4, water_fraction_over_cutoff = aa.calculate_water_fraction(weights3, param_dict)
    print(
        f"{water_fraction_over_cutoff.shape[0]:.3g} samples had a water fraction over the cutoff "
        f"of {param_dict['WATER_FRACTION_CUTOFF']:.3g}."
    )

    # Step 120
    print(f"\nStep 120")
    weights5 = aa.cutoff_water_fraction(weights4, param_dict)
    print(
        f"The water fraction was set to the cutoff of {param_dict['WATER_FRACTION_CUTOFF']:.3g} for those samples in downstream analysis."
    )

    # Step 121 and 122
    print(f"\nStep 121 and 122")
    weights6 = aa.calculate_effective_dry_stool_extracted_from(weights5)
    copies_weights = samples_controls2[samples_controls2["Type"] == "Sample"].merge(weights6, on="Name", how="left")
    copies_weights2 = aa.calculate_copies_per_wet_stool_g(copies_weights)
    copies_weights3 = aa.calculate_copies_per_dry_stool_g(copies_weights2)
    print(
        f"The 16S rRNA copies per wet gram, 16S rRNA copies per dry gram, and the effective dry stool extracted from (g) were calculated."
    )

    # Outputs
    with pd.ExcelWriter(os.path.join(output_folder_path, "univ_analysis_output.xlsx")) as writer:
        copies_weights3.to_excel(writer, sheet_name="from_universal_analysis", index=False)
        nist_measured_expected.to_excel(writer, sheet_name="nist_measured_expected", index=False)
        negative_DNA_extraction_controls.to_excel(writer, sheet_name="neg_dna_extract_controls", index=False)
        positive_DNA_extraction_controls.to_excel(writer, sheet_name="pos_dna_extract_controls", index=False)
        samples_with_DNA_extraction_input_outside_range.to_excel(
            writer, sheet_name="extract_input_outside_range", index=False
        )
        samples_with_drying_input_outside_range.to_excel(writer, sheet_name="drying_input_outside_range", index=False)
        samples_with_dry_amount_of_stool_low.to_excel(writer, sheet_name="dry_stool_amount_low", index=False)
        water_fraction_over_cutoff.to_excel(writer, sheet_name="water_fraction_over_cutoff", index=False)
    print(f"\nOutputs Complete\n")


if __name__ == "__main__":
    main()
