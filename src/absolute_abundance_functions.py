import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import seaborn as sns


# Function to configure matplotlib settings for consistent font sizes and properties
def prepare_matplotlib(sm=12, md=14, bg=16):
    SMALL_SIZE = sm
    MEDIUM_SIZE = md
    BIGGER_SIZE = bg

    plt.rc("font", size=SMALL_SIZE)  # controls default text sizes
    plt.rc("axes", titlesize=MEDIUM_SIZE)  # fontsize of the axes title
    plt.rc("axes", labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
    plt.rc("xtick", labelsize=SMALL_SIZE - 1)  # fontsize of the tick labels
    plt.rc("ytick", labelsize=SMALL_SIZE - 1)  # fontsize of the tick labels
    plt.rc("legend", fontsize=SMALL_SIZE)  # legend fontsize
    plt.rc("figure", titlesize=BIGGER_SIZE)  # fontsize of the figure title

    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = "Liberation Sans"


# Sorts by column number then row letter
def sort_by_well_96(df):
    cols = df.columns
    df["Letter"] = df["Well96"].str[0]
    df["Number"] = df["Well96"].str[1:]
    df.sort_values(by=["Number", "Letter"], ascending=[True, True], inplace=True)
    return df[cols].copy()


# Turns A1 to A01 while A10 stays as A10
def fix_well_96(st):
    if len(st) < 3:
        return st[0] + "0" + st[1]
    else:
        return st


def qpcr_initial_formatting(layout, format_con, qpcr_initial, param_dict):
    layout = layout.copy()
    format_con = format_con.copy()
    qpcr_initial = qpcr_initial.copy()

    layout["Dilution"] = layout["LiquidHandlerDilution"] * layout["SinglePipettorDilution"]
    layout384 = pd.merge(format_con, layout, how="left", on="Well96")

    if param_dict["NUM_OF_TECH_REPS"] == 3:
        layout384 = layout384[layout384.Replica.isin(["Rep1", "Rep2", "Rep3"])]
    elif param_dict["NUM_OF_TECH_REPS"] == 2:
        layout384 = layout384[layout384.Replica.isin(["Rep1", "Rep2"])]
    else:
        raise ValueError(f"{param_dict['NUM_OF_TECH_REPS']:.3g} number of tech reps is not an allowed value of 2 or 3.")

    qpcr_initial.rename(columns={"Well": "Well384"}, inplace=True)
    qpcr2 = pd.merge(layout384, qpcr_initial[["Well384", "Cq"]].copy(), how="left", on="Well384")

    qpcr2["Well96"] = qpcr2["Well96"].apply(fix_well_96)
    layout["Well96"] = layout["Well96"].apply(fix_well_96)
    return qpcr2, layout


def manage_tech_rep_failures(df1):
    df1 = df1.copy()
    df = df1[df1["Name"].notna()].copy()
    df["Cq"] = pd.to_numeric(df["Cq"], errors="coerce")
    counts = df.groupby("Name")["Cq"].count()
    filtered_names = counts[counts >= 2].index
    at_least_two_reps = df[df["Name"].isin(filtered_names)]
    fewer_than_two_reps = df[~df["Name"].isin(filtered_names)]
    return at_least_two_reps, fewer_than_two_reps


def assess_NTC(df, param_dict):
    df = df.copy()
    neg = df[df["Type"] == "PCRNeg"].copy()
    temp = neg.groupby(["Name"])["Cq"].median().reset_index()
    param_dict["CQ_LIMIT_OF_BLANK"] = temp["Cq"].min()
    Cq_NTC_span = temp["Cq"].max() - temp["Cq"].min()
    return Cq_NTC_span, param_dict


def calculate_qubit_copies(df, copynum):
    df = df.copy()
    df["Q_copies_rxn"] = (1.0 / (df["Dilution"])) * copynum
    return df


def special_y_axis_Cq(ax):
    y_min, y_max = ax.get_ylim()
    tick_positions = []
    current_tick = y_min
    while current_tick <= y_max:
        tick_positions.append(current_tick)
        current_tick += 3.3
    ax.set_yticks(tick_positions)


def visualize_standards(df_p, df_f):
    df_p = df_p.copy()
    df_f = df_f.copy()
    fig, ax1 = plt.subplots(1, 1, figsize=(5, 4))
    sns.scatterplot(
        data=pd.concat([df_p, df_f]),
        x="Q_copies_rxn",
        y="Cq",
        hue="Type",
        palette=sns.color_palette("colorblind", 2),
        ax=ax1,
    )
    plt.xscale("log")
    sns.despine()
    plt.xlabel("16S rRNA copies/rxn from Qubit HS")
    plt.ylabel("Cq")
    special_y_axis_Cq(ax1)
    plt.grid()


def remove_dilutions_with_tech_rep_variation(df_s, param_dict):
    df_s = df_s.copy()
    max_and_min = df_s.groupby("Dilution")["Cq"].agg(["max", "min"]).reset_index()
    max_and_min["Cq_span"] = max_and_min["max"] - max_and_min["min"]
    valid_dilutions = max_and_min[max_and_min["Cq_span"] < param_dict["MAX_CQ_SPAN_STANDARD_DIL_PT"]].Dilution.values

    thrown_out = df_s[~df_s["Dilution"].isin(valid_dilutions)].copy()
    thrown_out["reason_for_removal"] = "standard_plasmid_tech_rep_variation"
    return df_s[df_s["Dilution"].isin(valid_dilutions)].copy(), thrown_out


def remove_concentrated_standards(df_s, param_dict):
    df_s = df_s.copy()
    median_Cq_by_dilution = df_s.groupby("Dilution")["Cq"].median().to_frame().reset_index()
    median_Cq_by_dilution = median_Cq_by_dilution.sort_values(by="Dilution")
    median_Cq_by_dilution["Cq_diffs"] = median_Cq_by_dilution["Cq"].diff()
    diffs = pd.concat([median_Cq_by_dilution["Dilution"].shift(1), median_Cq_by_dilution["Dilution"]], axis=1)
    diffs.columns = ["Dilution_before", "Dilution"]
    diffs = diffs.merge(median_Cq_by_dilution, on="Dilution", how="left")
    diffs_select = diffs[diffs["Cq"] < param_dict["CQ_CUTOFF_CONC_STANDARDS"]].copy()
    diffs_select = diffs_select.drop(diffs_select.index[0])

    dilutions_to_remove = []
    for index, row in diffs_select.iterrows():
        if row["Cq_diffs"] < param_dict["MIN_CQ_GAP_CONC_STANDARDS"]:
            dilutions_to_remove.append(row["Dilution_before"])
        else:
            break

    thrown_out = df_s[df_s["Dilution"].isin(dilutions_to_remove)].copy()
    thrown_out["reason_for_removal"] = "standard_plasmid_too_concentrated"
    return df_s[~df_s["Dilution"].isin(dilutions_to_remove)].copy(), thrown_out


def remove_dilute_standards(df_s, param_dict):
    df_s = df_s.copy()
    max_Cq_by_dilution = df_s.groupby("Dilution")["Cq"].median().to_frame().reset_index()
    valid_dilutions = max_Cq_by_dilution[
        max_Cq_by_dilution["Cq"] < (param_dict["CQ_LIMIT_OF_BLANK"] - param_dict["CQ_STANDARDS_SEP_LOB"])
    ].Dilution.values

    thrown_out = df_s[~df_s["Dilution"].isin(valid_dilutions)].copy()
    thrown_out["reason_for_removal"] = "standard_plasmid_too_dilute"
    return df_s[df_s["Dilution"].isin(valid_dilutions)].copy(), thrown_out


def abs_reciprocal(x):
    if x < 1.0:
        return 1.0 / x
    else:
        return x


def compare_pvul_fpra(df_p, df_f):
    df_p = df_p.copy()
    df_f = df_f.copy()
    smodel_p = stats.linregress(np.log10(df_p["Q_copies_rxn"]), df_p["Cq"])
    smodel_f = stats.linregress(np.log10(df_f["Q_copies_rxn"]), df_f["Cq"])
    combined = pd.concat([df_p, df_f])
    min_Cq = combined["Cq"].min()
    max_Cq = combined["Cq"].min()
    min_Cq_p_to_f = (10 ** ((min_Cq - smodel_p.intercept) / smodel_p.slope)) / (
        10 ** ((min_Cq - smodel_f.intercept) / smodel_f.slope)
    )
    max_Cq_p_to_f = (10 ** ((max_Cq - smodel_p.intercept) / smodel_p.slope)) / (
        10 ** ((max_Cq - smodel_f.intercept) / smodel_f.slope)
    )
    return max(abs_reciprocal(min_Cq_p_to_f), abs_reciprocal(max_Cq_p_to_f))


def final_linear_regression(df_p, df_f, param_dict):
    df_p = df_p.copy()
    df_f = df_f.copy()
    combined = pd.concat([df_p, df_f])
    smodel = stats.linregress(np.log10(combined["Q_copies_rxn"]), combined["Cq"])
    param_dict["CQ_OF_MOST_DILUTE_STANDARD_POINT"] = (
        combined.groupby(["Type", "Dilution"])["Cq"].median().reset_index().Cq.max()
    )
    param_dict["CQ_OF_MOST_CONC_STANDARD_POINT"] = (
        combined.groupby(["Type", "Dilution"])["Cq"].median().reset_index().Cq.min()
    )
    return smodel, param_dict


def plot_standard_model(df_p, df_f, model):
    df = pd.concat([df_p, df_f])
    fig, ax1 = plt.subplots(1, 1, figsize=(5, 4))
    xsl = np.log10(df["Q_copies_rxn"])
    xs = df["Q_copies_rxn"]
    sns.scatterplot(data=df, x="Q_copies_rxn", y="Cq", hue="Type", palette=sns.color_palette("colorblind", 2), ax=ax1)
    sns.lineplot(x=xs, y=(model.intercept + model.slope * xsl), color="black", ax=ax1)
    plt.xscale("log")
    sns.despine()
    plt.xlabel("16S rRNA copies/rxn from Qubit HS")
    plt.ylabel("Cq")
    special_y_axis_Cq(ax1)
    plt.grid()


def calculate_pcr_efficiency(model):
    return 10**(-1/model.slope) - 1

def calculate_cq_limit_of_quantification(param_dict):
    if param_dict["OVERHANG_ALLOWED"] == True:
        param_dict["CQ_LIMIT_OF_QUANTIFICATION"] = (
            param_dict["CQ_LIMIT_OF_BLANK"] - param_dict["CQ_NON_STANDARDS_SEP_LOB"]
        )
    else:
        param_dict["CQ_LIMIT_OF_QUANTIFICATION"] = min(
            param_dict["CQ_LIMIT_OF_BLANK"] - param_dict["CQ_NON_STANDARDS_SEP_LOB"],
            param_dict["CQ_OF_MOST_DILUTE_STANDARD_POINT"],
        )
    return param_dict


def qpcr_calculate_copies_per_reaction_lob(model, param_dict):
    param_dict["COPIES_RXN_LIMIT_OF_BLANK"] = 10 ** ((param_dict["CQ_LIMIT_OF_BLANK"] - model.intercept) / model.slope)
    return param_dict


def identify_samples_with_tech_rep_variation(df, param_dict):
    df = df.copy()
    max_cq_sample_closest_two_reps = param_dict["MAX_CQ_DIFF_SAMPLE_CLOSEST_TWO_REPS"]
    if param_dict["NUM_OF_TECH_REPS"] == 2:
        max_cq_sample_closest_two_reps = max_cq_sample_closest_two_reps / 2.0

    max_median_min = df.groupby("Name")["Cq"].agg(["max", "median", "min"]).reset_index()
    max_median_min["max_to_median"] = max_median_min["max"] - max_median_min["median"]
    max_median_min["median_to_min"] = max_median_min["median"] - max_median_min["min"]
    max_median_min["closest_reps"] = max_median_min[["max_to_median", "median_to_min"]].min(axis=1)
    valid_names = max_median_min[max_median_min["closest_reps"] < max_cq_sample_closest_two_reps].Name.values

    thrown_out = df[~df["Name"].isin(valid_names)].copy()
    thrown_out["reason_for_removal"] = "tech_rep_variation"
    return df[df["Name"].isin(valid_names)].copy(), thrown_out


def calculate_medians(df, layout):
    df = df.copy()
    max_median_min = df.groupby("Name")["Cq"].agg(["median", "max", "min"]).reset_index()
    max_median_min.columns = ["Name", "Cq_median", "Cq_max", "Cq_min"]
    return max_median_min.merge(layout, on="Name", how="left")


def cq_falls_in_quantifiable_range(df, param_dict):
    df = df.copy()
    df2 = df[df["Cq_median"] <= param_dict["CQ_LIMIT_OF_QUANTIFICATION"]].copy()
    df2 = df2[df2["Cq_median"] >= param_dict["CQ_OF_MOST_CONC_STANDARD_POINT"]].copy()
    too_dilute = df[df["Cq_median"] > param_dict["CQ_LIMIT_OF_QUANTIFICATION"]].copy()
    too_conc = df[df["Cq_median"] < param_dict["CQ_OF_MOST_CONC_STANDARD_POINT"]].copy()

    too_dilute["reason_for_removal"] = "too_dilute"
    too_conc["reason_for_removal"] = "too_conc"
    return df2, too_dilute, too_conc


def is_low_confidence(df, param_dict):
    df = df.copy()
    cutoff = param_dict["CQ_LIMIT_OF_BLANK"] - param_dict["CQ_LOW_CONF_SEP_LOB"]
    low = df[df["Cq_median"] > cutoff].copy()
    low["low_confidence"] = True
    low_keep = low[low["Dilution"] <= 1.0].copy()
    low_reassay = low[low["Dilution"] > 1.0].copy()
    high = df[df["Cq_median"] <= cutoff].copy()
    high["low_confidence"] = False
    low_reassay["reason_for_removal"] = "low_conf_assay_with_less_dilution_possible"
    return sort_by_well_96(pd.concat([low_keep, high])), low_reassay


def qpcr_calculate_copies_from_df(df, model):
    df = df.copy()
    df["copies_reaction"] = 10 ** ((df["Cq_median"] - model.intercept) / model.slope)
    return df


def ddpcr_formatting(layout, ddpcr):
    layout = layout.copy()
    ddpcr = ddpcr.copy()

    layout["Dilution"] = layout["LiquidHandlerDilution"] * layout["SinglePipettorDilution"]
    layout["Well96"] = layout["Well96"].apply(fix_well_96)

    ddpcr = ddpcr.rename(columns={"Well": "Well96"}, inplace=False)
    ddpcr = ddpcr[["Well96", "Positives", "Negatives", "Accepted Droplets"]].copy()
    ddpcr2 = pd.merge(layout, ddpcr, how="left", on="Well96")
    return ddpcr2


def ddpcr_calculate_copies_per_reaction_as_setup(df, param_dict):
    df = df.copy()
    neg_is_zero = df[df["Negatives"] == 0].copy()
    neg_is_zero["copies_reaction"] = 20000000
    neg_not_zero = df[df["Negatives"] > 0].copy()
    neg_not_zero["copies_reaction"] = (
        np.log(neg_not_zero["Accepted Droplets"] / neg_not_zero["Negatives"])
        * (1.0 / param_dict["DROPLET_VOLUME"])
        * (1000.0)
        * param_dict["RXN_VOLUME"]
    )
    return sort_by_well_96(pd.concat([neg_is_zero, neg_not_zero]))


def check_num_of_droplets(df, param_dict):
    df = df.copy()
    enough_droplets = df[df["Accepted Droplets"] >= param_dict["MIN_ACCEPTED_DROPLETS"]].copy()
    not_enough_droplets = df[df["Accepted Droplets"] < param_dict["MIN_ACCEPTED_DROPLETS"]].copy()
    return enough_droplets, not_enough_droplets


def ddpcr_determine_limit_of_blank(df, param_dict):
    df = df.copy()
    param_dict["COPIES_RXN_LIMIT_OF_BLANK"] = df[df["Type"] == "PCRNeg"].copies_reaction.max()
    copies_rxn_NTC_span = (
        df[df["Type"] == "PCRNeg"].copies_reaction.max() / df[df["Type"] == "PCRNeg"].copies_reaction.min()
    )
    return copies_rxn_NTC_span, param_dict


def find_low_negative_droplets(df, param_dict):
    df = df.copy()
    too_conc = df[df["Negatives"] < param_dict["MIN_NEGATIVE_DROPLETS"]].copy()
    keep = df[df["Negatives"] >= param_dict["MIN_NEGATIVE_DROPLETS"]].copy()
    return keep, too_conc


def ddpcr_find_too_dilute(df, param_dict):
    df = df.copy()
    cutoff = param_dict["COPIES_RXN_LIMIT_OF_BLANK"] * param_dict["COPIES_RXN_LOQ_MULT"]
    too_dilute = df[df["copies_reaction"] < cutoff].copy()
    keep = df[df["copies_reaction"] >= cutoff].copy()
    return keep, too_dilute


def format_nist_expected(df):
    df = df.copy()
    df["copies_uL_expected"] = pd.to_numeric(df["copies_uL_expected"])
    return df


def calculate_nist_copies_uL(df):
    df = df.copy()
    df["copies_uL"] = df["copies_reaction"] * df["Dilution"] * (1.0 / df["uLAdded"])
    return df


def compare_nist_to_expected(df, nist_expect, param_dict):
    df = df.copy()
    nist_expect = nist_expect.copy()

    df = df.merge(nist_expect, on="Name", how="left")
    df["measured_to_expected"] = df["copies_uL"] / df["copies_uL_expected"]
    df_too_low = df[df["measured_to_expected"] < (1.0 / param_dict["NIST_MAX_FOLD_DIFF"])].copy()
    df_too_high = df[df["measured_to_expected"] > param_dict["NIST_MAX_FOLD_DIFF"]].copy()
    df_just_right = df[df["measured_to_expected"] >= (1.0 / param_dict["NIST_MAX_FOLD_DIFF"])].copy()
    df_just_right = df_just_right[df_just_right["measured_to_expected"] <= param_dict["NIST_MAX_FOLD_DIFF"]].copy()

    df_just_right["within_desired_range"] = "yes"
    df_too_low["within_desired_range"] = "no; measured is too low"
    df_too_high["within_desired_range"] = "no; measured is too high"

    return pd.concat([df_just_right, df_too_low, df_too_high])[
        ["Name", "ATCC ID", "copies_uL", "copies_uL_expected", "measured_to_expected", "within_desired_range"]
    ].copy()


def calculate_copies_per_dna_extraction(df):
    df = df.copy()
    df["copies_dna_extraction"] = df["copies_reaction"] * df["Dilution"] * (1.0 / df["uLAdded"]) * df["ElutionVolume"]
    return df


def assess_positive_DNA_extraction_controls(df):
    df = df.copy()
    posDNA = df[df["Type"] == "DNAPos"].copy()
    fold_diff = posDNA["copies_dna_extraction"].max() / posDNA["copies_dna_extraction"].min()
    return posDNA, fold_diff


def calculate_wet_mass_extracted_from(df, param_dict):
    df = df.copy()
    df["amt_dna_wet"] = df["filled_PB"] - df["empty_PB"]
    outside_amount_extracted_range = pd.concat(
        [
            df[df["amt_dna_wet"] > param_dict["EXTRACT_MAX_INPUT"]].copy(),
            df[df["amt_dna_wet"] < param_dict["EXTRACT_MIN_INPUT"]].copy(),
        ]
    )
    return df, outside_amount_extracted_range


def calculate_wet_and_dry_drying_mass(df, param_dict):
    df = df.copy()
    df["amt_dried_wet"] = df["filled_wt"] - df["empty_wt"]
    df["amt_dried_dry"] = df["dry_wt"] - df["empty_wt"]
    outside_amount_dried_range = pd.concat(
        [
            df[df["amt_dried_wet"] > param_dict["DRYING_MAX_INPUT"]].copy(),
            df[df["amt_dried_wet"] < param_dict["DRYING_MIN_INPUT"]].copy(),
        ]
    )
    dry_amount_of_stool_low = df[df["amt_dried_dry"] < param_dict["MIN_DRIED_DRY_MASS"]].copy()
    return df, outside_amount_dried_range, dry_amount_of_stool_low


def calculate_water_fraction(df, param_dict):
    df = df.copy()
    df["water_fraction"] = (df["amt_dried_wet"] - df["amt_dried_dry"]) / df["amt_dried_wet"]
    return df, df[df["water_fraction"] > param_dict["WATER_FRACTION_CUTOFF"]].copy()


def cutoff_water_fraction_helper(x, cutoff):
    if x > cutoff:
        return cutoff
    else:
        return x


def cutoff_water_fraction(df, param_dict):
    df = df.copy()
    df["cutoff_water_fraction"] = df["water_fraction"].apply(
        lambda x: cutoff_water_fraction_helper(x, cutoff=param_dict["WATER_FRACTION_CUTOFF"])
    )
    return df


def calculate_effective_dry_stool_extracted_from(df):
    df = df.copy()
    df["effective_amt_dna_dry"] = df["amt_dna_wet"] * (1.0 - df["water_fraction"])
    df["effective_cutoff_amt_dna_dry"] = df["amt_dna_wet"] * (1.0 - df["cutoff_water_fraction"])
    return df


def calculate_copies_per_wet_stool_g(df):
    df = df.copy()
    df["copies_wet_g"] = df["copies_dna_extraction"] / df["amt_dna_wet"]
    return df


def calculate_copies_per_dry_stool_g(df):
    df = df.copy()
    df["copies_dry_g"] = df["copies_wet_g"] * (1.0 / (1.0 - df["cutoff_water_fraction"]))
    return df
