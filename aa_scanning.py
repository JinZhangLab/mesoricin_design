import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import itertools
import re
from antifungal.design import single_point_mutate, predict_MIC

PROTEINOGENIC_AA = "ACDEFGHIKLMNPQRSTVWY"

def vis_mutation(mutation_prediction, original_seq, original_AFI, save_fig_path="scanning_profile.png"):
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
    save_fig_path : str or Path, optional
        Path to save the plot, default is "scanning_profile.png".
    
    Returns:
    -------
    None
        Saves a violin plot with overlayed scatter points showing the impact of amino acid 
        substitutions on antifungal activity.
    """

    mp = pd.DataFrame(mutation_prediction)

    # 1. Calculate the antifungal activity improvement for each mutation
    activity_improvement = pd.DataFrame(columns=["original residue", "position", "mutant residue", "activity_improvement"])

    pattern = re.compile(r"mutate_(\w+)_(\d+)_(\w+)")  # Regex to extract mutation position and amino acid
    for _, row in mp.iterrows():
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
    mesoricin2 = "YCRTYWRYGRLRRRCYRRR"
    original_AFI = predict_MIC([mesoricin2])["AFI"][0]

    single_point_mutate_instance = single_point_mutate(mesoricin2)
    single_point_mutate_predictions = single_point_mutate_instance.get_candidate_sequences().predict()

    vis_mutation(single_point_mutate_predictions, mesoricin2, original_AFI)