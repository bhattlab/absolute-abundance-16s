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
    "--ddpcr-path",
    "ddpcr_path",
    required=True,
    type=click.Path(exists=True),
    help="ddPCR data path to Excel file.",
    metavar="",
)
@click.option(
    "-ds",
    "--ddpcr-sheet",
    "ddpcr_sheet_name",
    required=True,
    type=str,
    help="ddPCR data Excel sheet name.",
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
    "-o",
    "--output-path",
    "output_folder_path",
    required=True,
    type=click.Path(exists=True),
    help="output folder path.",
    metavar="",
)
@click.option(
    "--droplet-volume",
    "DROPLET_VOLUME",
    default=0.795,
    type=click.FloatRange(0.75, 0.92),
    help="volume (nL) of each droplet.",
    metavar="",
)
@click.option(
    "--rxn-volume",
    "RXN_VOLUME",
    default=22,
    type=click.FloatRange(20, 30),
    help="volume (uL) of ddPCR reaction as setup on the pre-droplet plate.",
    metavar="",
)
@click.option(
    "--min-accepted-droplets",
    "MIN_ACCEPTED_DROPLETS",
    default=10000,
    type=click.IntRange(min=5000, max=20000),
    help="minimum number of accepted droplets per well.",
    metavar="",
)
@click.option(
    "--max-copies-rxn-lob",
    "MAX_COPIES_RXN_LOB",
    default=25,
    type=click.FloatRange(min=0, max=100),
    help="maximum 16S rRNA copies per reaction of the no template controls.",
    metavar="",
)
@click.option(
    "--max-copies-rxn-span-ntc",
    "MAX_COPIES_RXN_SPAN_NTC",
    default=4,
    type=click.FloatRange(min=1, max=10),
    help="maximum span (e.g. fold change) of 16S rRNA copies per reaction of all no template controls.",
    metavar="",
)
@click.option(
    "--min-negative-droplets",
    "MIN_NEGATIVE_DROPLETS",
    default=10,
    type=click.FloatRange(min=0, max=100),
    help="minimum number of negative droplets per well.",
    metavar="",
)
@click.option(
    "--copies-rxn-loq-mult",
    "COPIES_RXN_LOQ_MULT",
    default=4,
    type=click.FloatRange(min=2, max=100),
    help="multiplied by the limit of blank to define the limit of quantification, under which samples or controls are removed.",
    metavar="",
)
def main(ddpcr_path, ddpcr_sheet_name, layout_path, layout_sheet_name, output_folder_path, **param_dict):

    print(f"\nThe files being loaded into the script are")
    print(f"ddPCR data from {ddpcr_path} on sheet {ddpcr_sheet_name}")
    print(f"96-well layout from {layout_path} on sheet {layout_sheet_name}")
    print(f"output folder path is {output_folder_path}")

    print(f"\nThe parameters of the script are")
    for key, value in param_dict.items():
        click.echo(f"{key}: {value}")

    # Prepare matplotlib
    aa.prepare_matplotlib()

    # Load data
    layout96 = pd.read_excel(layout_path, sheet_name=layout_sheet_name)
    ddpcr = pd.read_excel(ddpcr_path, sheet_name=ddpcr_sheet_name)

    # Check that Name is unique
    vc = layout96["Name"].value_counts()
    if vc.max() > 1:
        vc = vc[vc > 1]
        raise ValueError(f"{len(vc)} sample names occur more than once: {', '.join(vc.index)}")

    # Formatting
    ddpcr2 = aa.ddpcr_formatting(layout96, ddpcr)

    # Step 103
    print(f"\nStep 103")
    ddpcr3 = aa.ddpcr_calculate_copies_per_reaction_as_setup(ddpcr2, param_dict)
    print(f"16S rRNA copies per reaction were calculated.")

    # Step 104
    print(f"\nStep 104")
    ddpcr4, droplets_low = aa.check_num_of_droplets(ddpcr3, param_dict)
    print(f"{droplets_low.shape[0]:.3g} is the number of wells with a low number of droplets.")

    # Step 105
    print(f"\nStep 105")
    copies_rxn_NTC_span, param_dict = aa.ddpcr_determine_limit_of_blank(ddpcr4, param_dict)
    if copies_rxn_NTC_span < param_dict["MAX_COPIES_RXN_SPAN_NTC"]:
        print(
            f"Copies per reaction of the no template controls span {copies_rxn_NTC_span:.3g}-fold"
            f", which is less than the maximum copies per reaction span of {param_dict['MAX_COPIES_RXN_SPAN_NTC']:.3g}-fold."
        )
        if param_dict["COPIES_RXN_LIMIT_OF_BLANK"] < param_dict["MAX_COPIES_RXN_LOB"]:
            print(
                f"Copies per reaction limit of blank is {param_dict['COPIES_RXN_LIMIT_OF_BLANK']:.3g},"
                f" which is less than the maximum allowed value of {param_dict['MAX_COPIES_RXN_LOB']:.3g}."
            )
        else:
            print(
                f"Copies per reaction limit of blank is {param_dict['COPIES_RXN_LIMIT_OF_BLANK']:.3g},"
                f" which is more than the maximum allowed value of {param_dict['MAX_COPIES_RXN_LOB']:.3g}."
            )
            raise RuntimeError("The script aborted. You may rerun it with different parameters, if appropriate.")
    else:
        print(
            f"Copies per reaction of the no template controls span {copies_rxn_NTC_span:.3g}-fold"
            f", which is more than the maximum copies per reaction span of {param_dict['MAX_COPIES_RXN_SPAN_NTC']:.3g}-fold."
        )
        raise RuntimeError("The script aborted. You may rerun it with different parameters, if appropriate.")

    # Step 106
    print(f"\nStep 106")
    ddpcr5 = ddpcr4[ddpcr4["Type"].isin(["PCRPos", "DNAPos", "DNANeg", "Sample"])].copy()
    ddpcr6, too_conc = aa.find_low_negative_droplets(ddpcr5, param_dict)
    print(
        f"{too_conc.shape[0]:.3g} is the number of samples, positive PCR NIST controls, negative DNA extraction controls,"
        f" and positive DNA extraction controls that are too concentrated, with fewer than {param_dict['MIN_NEGATIVE_DROPLETS']:.3g}"
        f" negative droplets."
    )

    # Step 107
    print(f"\nStep 107")
    ddpcr7, too_dilute = aa.ddpcr_find_too_dilute(ddpcr6, param_dict)
    print(
        f"{too_dilute.shape[0]:.3g} is the number of samples, positive PCR NIST controls, negative DNA extraction controls,"
        f" and positive DNA extraction controls that are too dilute, with fewer than "
        f"{param_dict['COPIES_RXN_LIMIT_OF_BLANK'] * param_dict['COPIES_RXN_LOQ_MULT']:.3g}"
        f" 16S rRNA copies per reaction."
    )

    # Outputs
    with pd.ExcelWriter(os.path.join(output_folder_path, "ddpcr_specific_output.xlsx")) as writer:
        ddpcr7.to_excel(writer, sheet_name="for_universal_analysis", index=False)
        droplets_low.to_excel(writer, sheet_name="too_few_droplets", index=False)
        too_conc.to_excel(writer, sheet_name="too_concentrated", index=False)
        too_dilute.to_excel(writer, sheet_name="too_dilute", index=False)
    print(f"\nOutputs Complete\n")


if __name__ == "__main__":
    main()
