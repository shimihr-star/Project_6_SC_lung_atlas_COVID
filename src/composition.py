
# ============================================================
# CELL-TYPE COMPOSITION ANALYSIS
# ============================================================
#
# PURPOSE:
# Calculate cell-type counts and proportions per donor.
#
# IMPORTANT:
# Each donor is treated as one biological replicate.
# The original AnnData object is not modified.
# ============================================================

import pandas as pd


# ======================================================================================================================
def calculate_cell_type_composition(adata):
    """
    Calculate the number and percentage of nuclei
    belonging to each cell type for every donor.

    Returns
    -------
    counts : pandas.DataFrame
        Rows: donors
        Columns: cell types
        Values: number of nuclei

    percentages : pandas.DataFrame
        Rows: donors
        Columns: cell types
        Values: percentage of nuclei per donor
    """

    # Check required metadata columns.
    required_columns = ["donor_id", "cell_type"]

    for column in required_columns:
        if column not in adata.obs.columns:
            raise ValueError(f"Missing metadata column: {column}")

    # Count nuclei of each cell type within each donor.
    counts = pd.crosstab(adata.obs["donor_id"],adata.obs["cell_type"],)

    # Convert categorical column labels to a regular Index.
    # This prevents compatibility issues in later pandas joins.
    counts.columns = pd.Index(counts.columns.astype(str))

    # Calculate the total number of nuclei per donor.
    total_nuclei = counts.sum(axis=1)

    # Convert counts to percentages within each donor.
    percentages = counts.div(total_nuclei, axis=0,) * 100
    return counts, percentages
# ======================================================================================================================


# ======================================================================================================================
def get_donor_metadata(adata):
    """
    Extract one clinical metadata record per donor.
    Verify consistent group and disease assignments.
    """

    required_columns = ["donor_id", "group", "disease"]

    for column in required_columns:
        if column not in adata.obs.columns:
            raise ValueError(f"Missing metadata column: {column}")

    donor_metadata = (adata.obs[required_columns].drop_duplicates().copy())

    if donor_metadata["donor_id"].duplicated().any():
        raise ValueError("A donor has inconsistent clinical metadata.")

    return donor_metadata.set_index("donor_id")
# ======================================================================================================================



# ======================================================================================================================
def summarize_composition_by_group(percentages, donor_metadata):
    """
    Calculate median cell-type percentages for COVID-19
    and control donors.

    Each donor contributes one observation.
    """

    # Match clinical metadata to the composition table.
    data = percentages.join(donor_metadata[["group"]],how="left",validate="one_to_one",)
    if data["group"].isna().any():
        raise ValueError("Some donors have no disease-group assignment.")

    # Calculate the median percentage per cell type and group.
    group_medians = data.groupby("group",observed=True,).median(numeric_only=True)

    # Transpose: rows become cell types.
    summary = group_medians.T

    # Calculate the COVID-19 minus Control difference.
    summary["difference"] = (summary["COVID-19"] - summary["Control"])

    # Sort by the largest positive difference.
    return summary.sort_values("difference", ascending=False,)
# ======================================================================================================================



# ======================================================================================================================
def plot_group_composition_barplot(summary, output_path):
    """
    Plot median cell-type percentages for COVID-19 and controls.
    Each bar represents the median across donors in that group.
    """

    import matplotlib.pyplot as plt
    import numpy as np

    # Select cell types to display.
    selected_cell_types = ["Macrophage", "AT2", "AT1 - putative", "T cell", "Fibroblast", "Alveolar macrophage",]

    # Extract group medians in the desired order.
    plot_data = summary.loc[selected_cell_types]
    y = np.arange(len(selected_cell_types))
    bar_height = 0.35
    fig, ax = plt.subplots(figsize=(11, 6))

    # COVID-19 bars.
    ax.barh(y - bar_height / 2, plot_data["COVID-19"], height=bar_height, label="COVID-19", color="#5799F7",)

    # Control bars.
    ax.barh(y + bar_height / 2, plot_data["Control"], height=bar_height, label="Control", color="#65C887",)
    ax.set_yticks(y)
    ax.set_yticklabels(selected_cell_types)

    # Display first cell type at the top.
    ax.invert_yaxis()
    ax.set_xlabel("Median percentage of nuclei per donor (%)")
    ax.set_title("Cell-Type Composition: COVID-19 vs Control")
    ax.legend(frameon=False)
    ax.grid(axis="x", alpha=0.2)
    ax.set_axisbelow(True)

    # Remove unnecessary borders.
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Group composition bar chart saved:", output_path)
# ======================================================================================================================


# ======================================================================================================================
def plot_composition_by_group(percentages, donor_metadata, output_path):
    """
    Plot donor-level cell-type percentages by disease group.
    Each dot represents one donor.
    """

    import matplotlib.pyplot as plt
    import seaborn as sns

    selected_cell_types = ["Macrophage","Alveolar macrophage","AT2","AT1 - putative",]

    # Combine cell-type proportions with clinical metadata.
    data = percentages[selected_cell_types].join(donor_metadata[["group"]],validate="one_to_one",)

    # Convert wide table into long format for plotting.
    data = data.reset_index().melt(id_vars=["donor_id", "group"], var_name="cell_type",value_name="percentage",)
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    for ax, cell_type in zip(axes.flat, selected_cell_types):
        subset = data[data["cell_type"] == cell_type]
        sns.boxplot(data=subset, x="group", y="percentage", order=["Control", "COVID-19"], ax=ax, showfliers=False,)
        sns.stripplot(data=subset, x="group", y="percentage", order=["Control", "COVID-19"], ax=ax, color="black", alpha=0.7, jitter=0.15,)
        ax.set_title(cell_type)
        ax.set_xlabel("")
        ax.set_ylabel("Percentage of nuclei")

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
# ======================================================================================================================



# ======================================================================================================================
def test_cell_type_composition(percentages, donor_metadata):
    """
    Compare cell-type percentages between COVID-19 and control donors.

    Uses:
    - Two-sided Mann-Whitney U test
    - Benjamini-Hochberg FDR correction

    Each donor is one independent biological observation.
    """

    import pandas as pd
    from scipy.stats import mannwhitneyu
    from statsmodels.stats.multitest import multipletests

    # Match clinical group labels to donor-level percentages.
    groups = donor_metadata["group"].reindex(percentages.index)
    if groups.isna().any():
        raise ValueError("Some donors are missing group information.")

    if not groups.isin(["COVID-19", "Control"]).all():
        raise ValueError("Unexpected disease-group label.")

    results = []

    for cell_type in percentages.columns:

        # Extract donor-level percentages for each group.
        covid_values = percentages.loc[groups == "COVID-19", cell_type].astype(float)
        control_values = percentages.loc[groups == "Control", cell_type].astype(float)

        # Two-sided Mann-Whitney U test.
        statistic, p_value = mannwhitneyu(covid_values, control_values, alternative="two-sided",method="auto",)

        # Calculate descriptive effect size.
        covid_median = covid_values.median()
        control_median = control_values.median()

        results.append({"cell_type": cell_type, "n_covid": len(covid_values), "n_control": len(control_values),
            "covid_median": covid_median, "control_median": control_median, "median_difference_pp": covid_median - control_median,
            "u_statistic": statistic,"p_value": p_value,})

    results_df = pd.DataFrame(results)

    # Correct for testing multiple cell types.
    results_df["fdr"] = multipletests(results_df["p_value"],alpha=0.05, method="fdr_bh",)[1]

    # Sort by statistical evidence.
    results_df = results_df.sort_values("fdr",ascending=True,)

    return results_df
# ======================================================================================================================



# ======================================================================================================================
def summarize_donor_characteristics(adata, counts):
    """
    Extract donor-level clinical characteristics and nucleus counts.

    Verify that each donor has consistent clinical metadata.
    """

    clinical_columns = ["donor_id", "group", "age", "sex", "pmi_h",]

    for column in clinical_columns:
        if column not in adata.obs.columns:
            raise ValueError(f"Missing metadata column: {column}")

    clinical = (adata.obs[clinical_columns].drop_duplicates().copy())

    # Verify one clinical record per donor.
    if clinical["donor_id"].duplicated().any():
        raise ValueError("A donor has inconsistent clinical metadata.")

    clinical = clinical.set_index("donor_id")

    # Add total recovered nuclei for each donor.
    clinical["n_nuclei"] = counts.sum(axis=1).reindex(clinical.index)

    if clinical["n_nuclei"].isna().any():
        raise ValueError("Missing nucleus counts for some donors.")

    return clinical
# ======================================================================================================================


def adjusted_composition_analysis(percentages, clinical):
    """
    Test associations between disease group and cell-type
    proportions while adjusting for age and sex.

    Uses OLS regression on logit-transformed proportions.
    """

    import numpy as np
    import pandas as pd
    import statsmodels.api as sm
    results = []

    # Focus on the two main biological findings.
    cell_types = ["Macrophage", "AT2"]
    metadata = clinical.copy()
    metadata["age"] = pd.to_numeric(metadata["age"], errors="coerce",)
    metadata["covid"] = (metadata["group"] == "COVID-19").astype(int)
    metadata["male"] = (metadata["sex"] == "male").astype(int)

    for cell_type in cell_types:
        data = metadata[["age", "covid", "male"]].join(percentages[[cell_type]], how="inner",)
        data = data.dropna()

        # Convert percentage to proportion.
        p = data[cell_type] / 100

        # Avoid undefined logit values at 0 and 1.
        epsilon = 1e-4
        p = p.clip(epsilon, 1 - epsilon)
        y = np.log(p / (1 - p))
        X = sm.add_constant(data[["covid", "age", "male"]])
        model = sm.OLS(y, X).fit(cov_type="HC3")
        results.append({"cell_type": cell_type, "n_donors": len(data), "covid_coefficient": model.params["covid"],
                        "covid_p_value": model.pvalues["covid"], "covid_ci_lower": model.conf_int().loc["covid", 0],
                        "covid_ci_upper": model.conf_int().loc["covid", 1], "r_squared": model.rsquared,})

    return pd.DataFrame(results)
