"""
This script analyzes the distribution of mesoricin structures generated through AF-Cluster,
t-SNE for dimensionality reduction, and KMeans for clustering. The results are visualized in a scatter plot
and saved to a CSV file.
"""
from pathlib import Path
import subprocess
import numpy as np
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import matplotlib
import re
import csv
from sklearn.cluster import KMeans
import concurrent.futures
from matplotlib.patches import Ellipse

matplotlib.use('Agg')

CURRENT_DIR = Path(__file__).resolve().parent
# Replace with the path to the TMalign executable
TMalign_PATH = Path("/Data/app/TMalign/TMalign").resolve()

# Verify TMalign executable exists
if not TMalign_PATH.is_file():
    raise FileNotFoundError(f"TMalign executable not found at {TMalign_PATH}")

def run_tmalign(pair):
    i, j, p1, p2 = pair
    if i == j:
        return i, j, 1.0
    try:
        res = subprocess.run([str(TMalign_PATH), str(p1), str(p2)], capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running TMalign on {p1} and {p2}: {e}")
        return i, j, 0.0
    matches = re.findall(r"TM-score=\s*([\d.]+)", res.stdout)
    if len(matches) >= 2:
        score = (float(matches[0]) + float(matches[1])) / 2
    elif len(matches) == 1:
        score = float(matches[0])
    else:
        score = 0.0
    return i, j, score

if __name__ == "__main__":
    mesoricin1_dir = CURRENT_DIR / "mesoricin1_dd7d1" / "pdb"
    mesoricin2_dir = CURRENT_DIR / "mesoricin2_b4283" / "pdb"
    mesoricin3_dir = CURRENT_DIR / "mesoricin3_b519a" / "pdb"
    mesoricin4_dir = CURRENT_DIR / "mesoricin4_b6c6e" / "pdb"

    # Collect all PDB files from the directories
    mesoricin1_files = list(mesoricin1_dir.glob("*.pdb"))
    mesoricin2_files = list(mesoricin2_dir.glob("*.pdb"))
    mesoricin3_files = list(mesoricin3_dir.glob("*.pdb"))
    mesoricin4_files = list(mesoricin4_dir.glob("*.pdb"))

    all_pdbs =  mesoricin1_files + mesoricin2_files + mesoricin3_files + mesoricin4_files
    labels = ["Mesoricin 1"]*len(mesoricin1_files) + ["Mesoricin 2"]*len(mesoricin2_files) + ["Mesoricin 3"]*len(mesoricin3_files) + ["Mesoricin 4"]*len(mesoricin4_files)

    # Check key data
    print(f"Number of Mesoricin 1 files: {len(mesoricin1_files)}")
    print(f"Number of Mesoricin 2 files: {len(mesoricin2_files)}")
    print(f"Number of Mesoricin 3 files: {len(mesoricin3_files)}")
    print(f"Number of Mesoricin 4 files: {len(mesoricin4_files)}")
    print(f"Total number of labels: {len(labels)}")
    print(f"Total number of PDB files: {len(all_pdbs)}")

    n = len(all_pdbs)
    tm_data = np.zeros((n, n))
    pairs = [(i, j, all_pdbs[i], all_pdbs[j]) for i in range(n) for j in range(i, n)]

    # Run TMalign in parallel to compute TM-scores
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(run_tmalign, pairs)
        for i, j, score in results:
            tm_data[i][j] = score
            if i != j:
                tm_data[j][i] = score

    # Perform t-SNE for dimensionality reduction
    tsne_input = 1 - tm_data  # Convert similarity to distance
    embeddings = TSNE(n_components=2, metric="precomputed", init="random", random_state=42).fit_transform(tsne_input)
    
    # Visualize the t-SNE embeddings
    fig, ax = plt.subplots()
    color_map = {
        "Mesoricin 1": "red",
        "Mesoricin 2": "blue",
        "Mesoricin 3": "green",
        "Mesoricin 4": "orange"
    }
    markers = {
        "Mesoricin 1": "o",
        "Mesoricin 2": "s",
        "Mesoricin 3": "^",
        "Mesoricin 4": "P"
    }

    for j, point in enumerate(embeddings):
        ax.scatter(
            point[0], 
            point[1], 
            color=color_map[labels[j]], 
            marker=markers[labels[j]], 
            s=20
        )

    # Explicitly define legend order
    label_order = ["Mesoricin 1", "Mesoricin 2", "Mesoricin 3", "Mesoricin 4"]
    handles = [
        plt.Line2D(
            [0], [0],
            marker=markers[label],
            color='w',
            label=label,
            markerfacecolor=color_map[label],
            markersize=10
        ) for label in label_order
    ]
    ax.legend(handles=handles, loc='best')
    plt.xlabel("t-SNE 1")
    plt.ylabel("t-SNE 2")

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Save the plot to a file
    plt.savefig(CURRENT_DIR / "figures" / "mesoricin_distribution.png", dpi=300)
    plt.close()