import os
import sys
import click
import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../src')))
import absolute_abundance_functions as aa


@click.command(context_settings=dict(max_content_width=120, show_default=True))
@click.option(
    "-q",
    "--qpcr-path",
    "qpcr_path",
    required=True,
    type=click.Path(exists=True),
    help="qPCR data path to Excel file.",
    metavar="",
)
@click.option(
    "-qs",
    "--qpcr-sheet",
    "qpcr_sheet_name",
    required=True,
    type=str,
    help="qPCR data Excel sheet name.",
    metavar="",
)
@click.option(
    "-l",
    "--layout-path",
    "layout_path",
    required=True,
    type=click.Path(exists=True),
    help="96-well layout path to Excel file.",
    metavar="",
)
@click.option(
    "-ls",
    "--layout-sheet",
    "layout_sheet_name",
    required=True,
    type=str,
    help="96-well layout Excel sheet name.",
    metavar="",
)
@click.option(
    "-f",
    "--format-conversion-path",
    "format_conversion_path",
    required=True,
    type=click.Path(exists=True),
    help="format conversion file path.",
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
    "--pcopy",
    "PCOPY",
    required=True,
    type=click.FloatRange(min=10**8, max=10**12),
    help="16S rRNA copies in uLAdded (e.g. 6 uL) of stock P. vulgatus standard plasmid.",
    metavar="",
)
@click.option(
    "--fcopy",
    "FCOPY",
    required=True,
    type=click.FloatRange(min=10**8, max=10**12),
    help="16S rRNA copies in uLAdded (e.g. 6 uL) of stock F. prausnitzii standard plasmid.",
    metavar="",
)
@click.option(
    "--num-of-tech-reps",
    "NUM_OF_TECH_REPS",
    default=3,
    type=click.IntRange(min=2, max=3),
    help="number of qPCR technical replicates. 2 replicates will select Rep1 and Rep2 (e.g. A1, A2) from the format conversion file and 3 replicates will select Rep1, Rep2, and Rep3 (e.g. A1, A2, B1) from the format conversion file.",
    metavar="",
)
@click.option(
    "--max-cq-span-ntc",
    "MAX_CQ_SPAN_NTC",
    default=2,
    type=click.FloatRange(min=0.1, max=3.3),
    help="maximum span of median Cq values of all no template controls.",
    metavar="",
)
@click.option(
    "--max-cq-span-standard-dil-pt",
    "MAX_CQ_SPAN_STANDARD_DIL_PT",
    default=2,
    type=click.FloatRange(min=0.1, max=3.3),
    help="maximum span of Cq values of each standard dilution point.",
    metavar="",
)
@click.option(
    "--max-standard-dil-pts-removed-tech-var",
    "MAX_STANDARD_DIL_PTS_REMOVED_TECH_VAR",
    default=1,
    type=click.IntRange(min=0, max=4),
    help="maximum number of standard dilution points removed for technical replicate variation.",
    metavar="",
)
@click.option(
    "--min-cq-gap-conc-standards",
    "MIN_CQ_GAP_CONC_STANDARDS",
    default=3.11,
    type=click.FloatRange(min=2.81, max=3.41),
    help="minimum Cq gap between concentrated points of the standard curve (e.g. no plateau).",
    metavar="",
)
@click.option(
    "--cq-cutoff-conc-standards",
    "CQ_CUTOFF_CONC_STANDARDS",
    default=15,
    type=click.FloatRange(min=8, max=18),
    help="maximum Cq value to consider gap between concentrated points of the standard curve.",
    metavar="",
)
@click.option(
    "--cq-standards-sep-lob",
    "CQ_STANDARDS_SEP_LOB",
    default=2,
    type=click.FloatRange(min=1, max=6.6),
    help="minimum Cq separation between most dilute standard dilution point and limit of blank.",
    metavar="",
)
@click.option(
    "--max-fold-change-pvul-fpra",
    "MAX_FOLD_CHANGE_PVUL_FPRA",
    default=2,
    type=click.FloatRange(min=0.1, max=4),
    help="maximum fold change between P. vulgatus and F. prausnitzii standard curves.",
    metavar="",
)
@click.option(
    "--most-steep-slope-allowed",
    "MOST_STEEP_SLOPE_ALLOWED",
    default=-3.58,
    type=click.FloatRange(min=-3.98, max=-3.31),
    help="most steep slope allowed for the standard curve.",
    metavar="",
)
@click.option(
    "--least-steep-slope-allowed",
    "LEAST_STEEP_SLOPE_ALLOWED",
    default=-3.11,
    type=click.FloatRange(min=-3.29, max=-2.71),
    help="least steep slope allowed for the standard curve.",
    metavar="",
)
@click.option(
    "--min-r-squared",
    "MIN_R_SQUARED",
    default=0.98,
    type=click.FloatRange(min=0.93, max=0.999999),
    help="minimum R squared value for the standard curve.",
    metavar="",
)
@click.option(
    "--cq-non-standards-sep-lob",
    "CQ_NON_STANDARDS_SEP_LOB",
    default=2,
    type=click.FloatRange(min=1, max=6.6),
    help="minimum Cq separation between samples or controls and limit of blank.",
    metavar="",
)
@click.option(
    "--overhang-allowed",
    "OVERHANG_ALLOWED",
    default=False,
    type=bool,
    help="determines whether samples are allowed to overhang the dilute end of the standard curve, while remaining the minimum distance from the limit of blank.",
    metavar="",
)
@click.option(
    "--max-copies-rxn-lob",
    "MAX_COPIES_RXN_LOB",
    default=500,
    type=click.FloatRange(min=5, max=3000),
    help="maximum acceptable apparent 16S rRNA copies per reaction of the no template controls.",
    metavar="",
)
@click.option(
    "--max-cq-diff-sample-closest-two-reps",
    "MAX_CQ_DIFF_SAMPLE_CLOSEST_TWO_REPS",
    default=2,
    type=click.FloatRange(min=0.1, max=3.3),
    help="maximum difference between the Cq values of the two closest technical replicates for samples or controls.",
    metavar="",
)
@click.option(
    "--cq-low-conf-sep-lob",
    "CQ_LOW_CONF_SEP_LOB",
    default=3.3,
    type=click.FloatRange(min=1, max=6.6),
    help="Cq separation between samples or controls and limit of blank to decide a measurement is low confidence.",
    metavar="",
)
def main(
    qpcr_path, qpcr_sheet_name, layout_path, layout_sheet_name, format_conversion_path, output_folder_path, **param_dict
):

    print(f"\nThe files being loaded into the script are")
    print(f"qPCR data from {qpcr_path} on sheet {qpcr_sheet_name}")
    print(f"96-well layout from {layout_path} on sheet {layout_sheet_name}")
    print(f"96-well to 384-well format conversion TSV from {format_conversion_path}")
    print(f"output folder path is {output_folder_path}")

    print(f"\nThe parameters of the script are")
    for key, value in param_dict.items():
        click.echo(f"{key}: {value}")

    # Prepare matplotlib
    aa.prepare_matplotlib()

    # Load data
    format_conversion = pd.read_csv(format_conversion_path, sep="\t", names=["Well96", "Replica", "Well384"])
    layout96 = pd.read_excel(layout_path, sheet_name=layout_sheet_name)
    qpcr = pd.read_excel(qpcr_path, sheet_name=qpcr_sheet_name)

    # Check that Name is unique
    vc = layout96["Name"].value_counts()
    if vc.max() > 1:
        vc = vc[vc > 1]
        raise ValueError(f"{len(vc)} sample names occur more than once: {', '.join(vc.index)}")

    # Formatting
    qpcr2, layout96 = aa.qpcr_initial_formatting(layout96, format_conversion, qpcr, param_dict)

    # Step 72
    print(f"\nStep 72")
    qpcr3, fewer_than_two_reps = aa.manage_tech_rep_failures(qpcr2)
    print(
        f"{fewer_than_two_reps['Name'].nunique():.3g} measurements (e.g. wells on the 96-well plate) "
        f"were removed from analysis since there were fewer than two technical replicates."
        f" Please carefully review them."
    )

    # Step 73
    print(f"\nStep 73")
    Cq_NTC_span, param_dict = aa.assess_NTC(qpcr3, param_dict)
    if Cq_NTC_span < param_dict["MAX_CQ_SPAN_NTC"]:
        print(
            f"Cqs of the no template controls span {Cq_NTC_span:.3g}"
            f", which is less than the maximum Cq span of {param_dict['MAX_CQ_SPAN_NTC']:.3g}"
        )
        print(f"The Cq limit of blank is {param_dict['CQ_LIMIT_OF_BLANK']:.3g}")
    else:
        print(
            f"Cqs of the no template controls span {Cq_NTC_span:.3g}"
            f", which is greater than the maximum Cq span of {param_dict['MAX_CQ_SPAN_NTC']:.3g}."
        )
        raise RuntimeError("The script aborted. You may rerun it with different parameters, if appropriate.")

    # Step 75
    print(f"\nStep 75")
    pvul1 = qpcr3[qpcr3["Type"].isin(["Pvul"])].copy()
    fpra1 = qpcr3[qpcr3["Type"].isin(["Fpra"])].copy()
    pvul2 = aa.calculate_qubit_copies(pvul1, param_dict["PCOPY"])
    fpra2 = aa.calculate_qubit_copies(fpra1, param_dict["FCOPY"])
    print(f"The P. vulgatus stock plasmid copy number per reaction is {param_dict['PCOPY']:.3g}")
    print(f"The F. prausnitzii stock plasmid copy number per reaction is {param_dict['FCOPY']:.3g}")
    print(f"Check that these values are correct as any discrepancies will affect all downstream results.")

    # Step 76
    print(f"\nStep 76")
    aa.visualize_standards(df_p=pvul2, df_f=fpra2)
    plt.savefig(
        os.path.join(output_folder_path, "standard_curve_visualization_post_step_76.pdf"),
        format="pdf",
        dpi=300,
        bbox_inches="tight",
    )
    print(f"Plot outputted with all dilution points of the standard curve plasmids.")

    # Step 77
    print(f"\nStep 77")
    pvul3, p_var_dilutions_removed = aa.remove_dilutions_with_tech_rep_variation(pvul2, param_dict)
    fpra3, f_var_dilutions_removed = aa.remove_dilutions_with_tech_rep_variation(fpra2, param_dict)
    p_var_dilutions_removed_number = p_var_dilutions_removed["Dilution"].nunique()
    f_var_dilutions_removed_number = f_var_dilutions_removed["Dilution"].nunique()
    print(
        f"{p_var_dilutions_removed_number:.3g}"
        f" dilution points of the P. vulgatus standard were removed due to technical variation."
    )
    print(
        f"{f_var_dilutions_removed_number:.3g}"
        f" dilution points of the F. prausnitzii standard were removed due to technical variation."
    )
    if (p_var_dilutions_removed_number + f_var_dilutions_removed_number) > param_dict[
        "MAX_STANDARD_DIL_PTS_REMOVED_TECH_VAR"
    ]:
        print(f"This exceeds the maximum of {param_dict['MAX_STANDARD_DIL_PTS_REMOVED_TECH_VAR']:.3g}")
        raise RuntimeError("The script aborted. You may rerun it with different parameters, if appropriate.")
    else:
        print(
            f"Combined this is less than or equal to the maximum of "
            f"{param_dict['MAX_STANDARD_DIL_PTS_REMOVED_TECH_VAR']:.3g}"
        )
    aa.visualize_standards(df_p=pvul3, df_f=fpra3)
    plt.savefig(
        os.path.join(output_folder_path, "standard_curve_visualization_post_step_77.pdf"),
        format="pdf",
        dpi=300,
        bbox_inches="tight",
    )

    # Step 78
    print(f"\nStep 78")
    pvul4, p_conc_dilutions_removed = aa.remove_concentrated_standards(pvul3, param_dict)
    fpra4, f_conc_dilutions_removed = aa.remove_concentrated_standards(fpra3, param_dict)
    print(
        f"{p_conc_dilutions_removed['Dilution'].nunique():.3g}"
        f" dilution points of the P. vulgatus standard were removed due to being too concentrated."
    )
    print(
        f"{f_conc_dilutions_removed['Dilution'].nunique():.3g}"
        f" dilution points of the F. prausnitzii standard were removed due to being too concentrated."
    )
    aa.visualize_standards(df_p=pvul4, df_f=fpra4)
    plt.savefig(
        os.path.join(output_folder_path, "standard_curve_visualization_post_step_78.pdf"),
        format="pdf",
        dpi=300,
        bbox_inches="tight",
    )

    # Step 79
    print(f"\nStep 79")
    pvul5, p_dil_dilutions_removed = aa.remove_dilute_standards(pvul4, param_dict)
    fpra5, f_dil_dilutions_removed = aa.remove_dilute_standards(fpra4, param_dict)
    print(
        f"{p_dil_dilutions_removed['Dilution'].nunique():.3g}"
        f" dilution points of the P. vulgatus standard were removed due to being too dilute."
    )
    print(
        f"{f_dil_dilutions_removed['Dilution'].nunique():.3g}"
        f" dilution points of the F. prausnitzii standard were removed due to being too dilute."
    )
    aa.visualize_standards(df_p=pvul5, df_f=fpra5)
    plt.savefig(
        os.path.join(output_folder_path, "standard_curve_visualization_post_step_79.pdf"),
        format="pdf",
        dpi=300,
        bbox_inches="tight",
    )

    # Step 80
    print(f"\nStep 80")
    maximum_pvul_fpra_discrepancy = aa.compare_pvul_fpra(pvul5, fpra5)
    if maximum_pvul_fpra_discrepancy < param_dict["MAX_FOLD_CHANGE_PVUL_FPRA"]:
        print(
            f"The maximum discrepancy between Pvul and Fpra is: "
            f"{maximum_pvul_fpra_discrepancy:.3g}"
            f" , which is less than the maximum acceptable fold-discrepancy, which is "
            f"{param_dict['MAX_FOLD_CHANGE_PVUL_FPRA']:.3g}"
        )
    else:
        print(
            f"The maximum discrepancy between Pvul and Fpra is: {maximum_pvul_fpra_discrepancy:.3g}"
            f" , which is greater than the maximum acceptable fold-discrepancy, which is "
            f"{param_dict['MAX_FOLD_CHANGE_PVUL_FPRA']:.3g}"
        )
        raise RuntimeError("The script aborted. You may rerun it with different parameters, if appropriate.")

    # Step 81
    print(f"\nStep 81")
    smodel, param_dict = aa.final_linear_regression(pvul5, fpra5, param_dict)
    aa.plot_standard_model(pvul5, fpra5, smodel)
    plt.savefig(
        os.path.join(output_folder_path, "standard_curve_visualization_final_post_step_81.pdf"),
        format="pdf",
        dpi=300,
        bbox_inches="tight",
    )
    print(
        f"The median Cq of the most dilute standard dilution point {param_dict['CQ_OF_MOST_DILUTE_STANDARD_POINT']:.3g}"
    )
    print(
        f"The median Cq of the most concentrated standard dilution point {param_dict['CQ_OF_MOST_CONC_STANDARD_POINT']:.3g}"
    )
    print(f"The y-intercept is: {smodel.intercept:.3g}")
    print(f"The slope is: {smodel.slope:.3g}")
    if smodel.slope < param_dict["LEAST_STEEP_SLOPE_ALLOWED"] and smodel.slope > param_dict["MOST_STEEP_SLOPE_ALLOWED"]:
        print(
            f"The slope is between {param_dict['MOST_STEEP_SLOPE_ALLOWED']:.3g}"
            f" and {param_dict['LEAST_STEEP_SLOPE_ALLOWED']:.3g}."
        )
    else:
        print(
            f"The slope is not between {param_dict['MOST_STEEP_SLOPE_ALLOWED']:.3g}"
            f" and {param_dict['LEAST_STEEP_SLOPE_ALLOWED']:.3g}."
        )
        raise RuntimeError("The script aborted. You may rerun it with different parameters, if appropriate.")
    print(f"The PCR efficiency is: {aa.calculate_pcr_efficiency(smodel):.3g}")
    print(f"The R squared value is: {smodel.rvalue * smodel.rvalue:.3g}")
    if smodel.rvalue * smodel.rvalue > param_dict["MIN_R_SQUARED"]:
        print(f"The R squared value is above the minimum value of {param_dict['MIN_R_SQUARED']:.3g}.")
    else:
        print(f"The R squared value is below the minimum value of {param_dict['MIN_R_SQUARED']:.3g}.")
        raise RuntimeError("The script aborted. You may rerun it with different parameters, if appropriate.")

    # Step 82
    print(f"\nStep 82")
    param_dict = aa.calculate_cq_limit_of_quantification(param_dict)
    print(f"The Cq lower limit of quantification is {param_dict['CQ_LIMIT_OF_QUANTIFICATION']:.3g}")
    if param_dict["CQ_LIMIT_OF_QUANTIFICATION"] > param_dict["CQ_OF_MOST_DILUTE_STANDARD_POINT"]:
        print(
            f"Note that the Cq limit of quantification overhangs the Cq of the most dilute standard curve point, which is "
            f"{param_dict['CQ_OF_MOST_DILUTE_STANDARD_POINT']:.3g}."
        )
        print(f"In some situations, this may be permissible. Proceed with caution.")
    else:
        print(
            f"This is more conservative than the Cq of the most dilute standard curve point, which is "
            f"{param_dict['CQ_OF_MOST_DILUTE_STANDARD_POINT']:.3g}"
        )

    # Step 83
    print(f"\nStep 83")
    param_dict = aa.qpcr_calculate_copies_per_reaction_lob(smodel, param_dict)
    print(f"The copies per reaction limit of blank is: {param_dict['COPIES_RXN_LIMIT_OF_BLANK']:.3g}")
    if param_dict["COPIES_RXN_LIMIT_OF_BLANK"] <= param_dict["MAX_COPIES_RXN_LOB"]:
        print(
            f"This is less than the maximum acceptable value, which is {param_dict['MAX_COPIES_RXN_LOB']:.3g}"
        )
    else:
        print(
            f"This is greater than the maximum acceptable value, which is {param_dict['MAX_COPIES_RXN_LOB']:.3g}"
        )
        raise RuntimeError("The script aborted. You may rerun it with different parameters, if appropriate.")

    # Step 84
    print(f"\nStep 84")
    samples_controls = qpcr3[qpcr3["Type"].isin(["PCRPos", "DNAPos", "DNANeg", "Sample"])]
    samples_controls2, samples_controls_high_variation = aa.identify_samples_with_tech_rep_variation(
        samples_controls, param_dict
    )
    print(
        f"{samples_controls_high_variation['Name'].nunique():.3g} samples, NIST positive PCR controls, "
        f"negative DNA extraction controls, and/or Zymo mock positive DNA extraction controls were removed due "
        f"to technical replicate variation."
    )

    # Step 85
    print(f"\nStep 85")
    samples_controls3 = aa.calculate_medians(samples_controls2, layout96)
    print(f"Medians were calculated.")

    # Step 86
    print(f"\nStep 86")
    samples_controls4, samples_controls_too_dilute, samples_controls_too_conc = aa.cq_falls_in_quantifiable_range(
        samples_controls3, param_dict
    )
    print(
        f"{samples_controls_too_dilute.shape[0]:.3g} samples, NIST positive PCR controls, "
        f"negative DNA extraction controls, and/or Zymo mock positive DNA extraction controls were removed because they are too dilute "
        f"with a Cq greater than {param_dict['CQ_LIMIT_OF_QUANTIFICATION']:.3g}."
    )
    print(
        f"{samples_controls_too_conc.shape[0]:.3g} samples, NIST positive PCR controls, "
        f"negative DNA extraction controls, and/or Zymo mock positive DNA extraction controls were removed because they are too concentrated "
        f"with a Cq less than {param_dict['CQ_OF_MOST_CONC_STANDARD_POINT']:.3g}."
    )

    # Step 87
    print(f"\nStep 87")
    samples_controls5, low_confidence_not_undiluted = aa.is_low_confidence(samples_controls4, param_dict)
    print(
        f"{low_confidence_not_undiluted.shape[0]:.3g} were low confidence (e.g. Cq within "
        f"{param_dict['CQ_LOW_CONF_SEP_LOB']:.3g} Cq of the Cq limit of blank) and not assayed without dilution. "
        f"We suggest reassaying these samples with less or no dilution."
    )

    # Step 88
    print(f"\nStep 88")
    samples_controls6 = aa.qpcr_calculate_copies_from_df(samples_controls5, smodel)
    print(f"16S rRNA copies per reaction were calculated.")

    # Outputs
    standard_dilution_points_removed = pd.concat(
        [
            p_conc_dilutions_removed,
            f_conc_dilutions_removed,
            p_dil_dilutions_removed,
            f_dil_dilutions_removed,
            p_var_dilutions_removed,
            f_var_dilutions_removed,
        ]
    )
    low_confidence_not_undiluted2 = low_confidence_not_undiluted.drop("low_confidence", axis=1)
    samples_removed = pd.concat([samples_controls_too_dilute, samples_controls_too_conc, low_confidence_not_undiluted2])
    with pd.ExcelWriter(os.path.join(output_folder_path, "qpcr_specific_output.xlsx")) as writer:
        samples_controls6.to_excel(writer, sheet_name="for_universal_analysis", index=False)
        fewer_than_two_reps.to_excel(writer, sheet_name="removed_wells_tech_reps_fails", index=False)
        standard_dilution_points_removed.to_excel(writer, sheet_name="removed_standards", index=False)
        samples_controls_high_variation.to_excel(writer, sheet_name="removed_high_variation_samples", index=False)
        samples_removed.to_excel(writer, sheet_name="removed_samples", index=False)
    print(f"\nOutputs Complete\n")


if __name__ == "__main__":
    main()
