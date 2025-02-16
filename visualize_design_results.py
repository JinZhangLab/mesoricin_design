import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import itertools
import re
from antifungal.design import segment, single_point_mutate, predict_MIC
from matplotlib.patches import Rectangle

PROTEINOGENIC_AA = "ACDEFGHIKLMNPQRSTVWY"


def vis_segment(segment_prediction, original_seq, original_AFI, save_fig_path):
    """
    Visualizes the impact of peptide sequence segmentation on antifungal activity through an exhaustive search.

    Parameters:
    ----------
    segment_prediction : dict
        A dictionary containing the prediction results for each segment of the peptide sequence.
    original_seq : str
        The original peptide sequence being analyzed.
    original_AFI : float
        The antifungal activity index (AFI) of the original peptide sequence.
    save_fig_path : str or Path
        The file path where the heatmap plot will be saved.
    """
    
    seg_pred_df = pd.DataFrame(segment_prediction)

    # 1. Calculate the antifungal activity improvement (relative AFI decrease) for each mutation
    activity_improvement = pd.DataFrame(columns=["start", "end", "activity_improvement"])

    pattern = re.compile(r"segment_(\d+)_(\d+)")  # Regex to extract segment start and end position
    for idx, row in seg_pred_df.iterrows():
        start, end = pattern.search(row["seq_name"]).groups()
        AFI_decrease = original_AFI - float(row["AFI"])
        AFI_decrease = ((AFI_decrease) / (original_AFI + 1e-12)) * 100  # Prevent division by zero
        activity_improvement = pd.concat(
            (
                activity_improvement,
                pd.DataFrame(
                    [[start, end, AFI_decrease]], 
                    columns=["start", "end", "activity_improvement"]
                )
            ),
            ignore_index=True
        )

    # Convert start and end to numeric
    activity_improvement["start"] = pd.to_numeric(activity_improvement["start"])
    activity_improvement["end"] = pd.to_numeric(activity_improvement["end"])

    # Get unique sorted start and end positions
    unique_starts = sorted(activity_improvement["start"].unique())
    unique_ends = sorted(activity_improvement["end"].unique())

    # Create mappings for positions to indices
    start_to_index = {start: index for index, start in enumerate(unique_starts)}
    end_to_index = {end: index for index, end in enumerate(unique_ends)}

    # Create a matrix to store the AFI decrease values, initialized with NaN
    heatmap_data = np.full((len(unique_ends), len(unique_starts)), np.nan)

    # Populate the matrix with AFI decrease values using the mappings
    for _, row in activity_improvement.iterrows():
        start_index = start_to_index[row["start"]]
        end_index = end_to_index[row["end"]]
        heatmap_data[end_index, start_index] = row["activity_improvement"]  # end is y, start is x


    # 2. Visualization
    if save_fig_path:
        plt.figure(figsize=(10, 8))  # Adjust figure size for better readability

        # Determine the symmetric range for color scale, centered at 0
        abs_max_value = np.nanmax(np.abs(heatmap_data))
        vmin = -abs_max_value
        vmax = abs_max_value

        # Identify the optimal segmentation: the cell with the maximum AFI decrease
        max_value = np.nanmax(heatmap_data)
        if np.isfinite(max_value):
            max_indices = np.where(heatmap_data == max_value)
            if len(max_indices[0]) > 0:
                max_end_index = max_indices[0][0]  # row index
                max_start_index = max_indices[1][0]  # column index

        # Generate the heatmap with proper labels, using unique start and end positions.
        # Reduced annotation font size using annot_kws.
        ax = sns.heatmap(
            heatmap_data,
            cmap='coolwarm',  # or 'RdBu_r'
            annot=True,
            fmt=".1f",
            annot_kws={'fontsize': 8},
            xticklabels=unique_starts,  # Unique Start positions
            yticklabels=unique_ends,    # Unique End positions
            cbar_kws={"label": "Relative AFI decrease (%)"},
            vmin=vmin, vmax=vmax,  # Center the color scale at 0
            linewidths=.5, linecolor="lightgrey",  # Optional: add lines to separate cells
            mask=np.isnan(heatmap_data)  # Mask NaN values to make them transparent
        )

        # Highlight the optimal segmentation cell with a rectangle (requires matplotlib.patches)
        if 'max_start_index' in locals() and 'max_end_index' in locals():
            ax.add_patch(Rectangle((max_start_index, max_end_index), 1, 1,
                                   fill=False, edgecolor='black', lw=3))

            # Add an annotation outside the plot indicating the optimal segmentation cell
            plt.figtext(0.95, 0.01, "Optimal segmentation", horizontalalignment="right",
                        fontsize=12, bbox=dict(facecolor='yellow', alpha=0.5, edgecolor='black'))

        plt.xlabel("Start Position", fontsize=14)
        plt.ylabel("End Position", fontsize=14)
        plt.title("Segmentation Analysis", fontsize=16)

        plt.savefig(save_fig_path, bbox_inches='tight', dpi=300)
        plt.close()  # close the figure to prevent display in notebook



def vis_single_point_mutation(mutation_prediction, original_seq, original_AFI, save_fig_path):
    """
    Analyze amino acid substitutions across a peptide sequence to identify positions and 
    mutations that enhance antifungal activity.

    This function performs simulations of all possible single amino acid substitutions 
    within a peptide sequence and evaluates their impact on antifungal activity (AFI). 
    It is specifically designed to analyze single-point mutations only, providing insights 
    into residue-specific contributions to activity enhancement. A higher percentage 
    decrease in AFI indicates a mutation with greater potential to optimize activity.

    Parameters:
    ----------
    mutation_prediction : dict
        A dictionary containing the prediction results of single point mutation for the 
        original_seq.
    original_seq : str
        The original peptide sequence to be analyzed.
    original_AFI : float
        The antifungal activity index (AFI) of the original peptide sequence.
    save_fig_path : str or Path
        Path to save the plot.
    
    Returns:
    -------
    None
        Saves a violin plot with overlayed scatter points showing the impact of amino acid 
        substitutions on antifungal activity.
    """

    mut_pred_df = pd.DataFrame(mutation_prediction)

    # 1. Calculate the antifungal activity improvement for each mutation
    activity_improvement = pd.DataFrame(columns=["original residue", "position", "mutant residue", "activity_improvement"])

    pattern = re.compile(r"mutate_(\w+)_(\d+)_(\w+)")  # Regex to extract mutation position and amino acid
    for _, row in mut_pred_df.iterrows():
        # Original amino acid, mutation position, and mutated amino acid
        original_aa, pos, mut_aa = pattern.search(row["seq_name"]).groups()
        
        # Since lower values indicate stronger antifungal activity, a decreased AFI indicates improved antifungal activity
        # Formula: ((original_index - mutated_index ) / original_index) * 100%
        AFI_decrease = original_AFI - float(row["AFI"])
        AFI_decrease = ((AFI_decrease) / (original_AFI + 1e-12)) * 100  # Prevent division by zero
        activity_improvement = pd.concat(
            (
                activity_improvement,
                pd.DataFrame(
                    [[original_aa, int(pos), mut_aa, AFI_decrease]], 
                    columns=[
                        "original residue",
                        "position",
                        "mutant residue",
                        "activity_improvement"
                    ]
                )
            ),
            ignore_index=True
        )

    # 2. Visualization
    if save_fig_path:
        colors = sns.color_palette("Set2", 4)
        markers = ['o', 's', 'D', '^', 'P']
        color_marker_combinations = list(itertools.product(colors, markers))

        amino_acid_styles = {
            aa: style for aa, style in zip(PROTEINOGENIC_AA, color_marker_combinations)
        }

        plt.figure(figsize=(14, 8))
        sns.violinplot(x="position", y="activity_improvement", data=activity_improvement,
                        inner=None, color="skyblue", linewidth=0, alpha=0.5)

        # Plot scatter plot with colors and markers, and add some jitter
        for aa, group in activity_improvement.groupby("mutant residue"):
            color, marker = amino_acid_styles[aa]
            jittered_positions = group["position"]-1 + np.random.uniform(-0.1, 0.1, size=len(group))
            plt.scatter(
                jittered_positions, group["activity_improvement"],
                label=aa, color=color, marker=marker, s=20, alpha=1.0
            )

        plt.legend(
            title="Mutant Residue",
            fontsize=10, ncol=4
        )

        plt.axhline(y=0, color='r', linestyle='--', linewidth=1)

        plt.xlabel("Sequence Position", fontsize=14)
        plt.ylabel(f"Relative AFI decrease (%)", fontsize=14)
        plt.title(f"Mutation Impact Analysis for Mesoricin2", fontsize=14)

        plt.xticks(range(len(original_seq)), list(original_seq))
        plt.grid(axis='y', linestyle='--', alpha=0.7)
                
        plt.savefig(save_fig_path, bbox_inches='tight', dpi=300)

if __name__ == "__main__":
    mesoricin1 = "RRYCRTYWRYGRLRRRCYRRRVWIWFRL"
    original_AFI_mesoricin1 = predict_MIC([mesoricin1])["AFI"][0]
    segment_instance = segment(mesoricin1)
    segment_predictions = segment_instance.get_candidate_sequences().predict()
    vis_segment(segment_predictions, mesoricin1, original_AFI_mesoricin1, "segmentation_analysis.png")


    mesoricin2 = "YCRTYWRYGRLRRRCYRRR"
    original_AFI_mesoricin2 = predict_MIC([mesoricin2])["AFI"][0]
    single_point_mutate_instance = single_point_mutate(mesoricin2)
    single_point_mutate_predictions = single_point_mutate_instance.get_candidate_sequences().predict()

    vis_single_point_mutation(single_point_mutate_predictions, mesoricin2, original_AFI_mesoricin2, "single_point_mutation.png")