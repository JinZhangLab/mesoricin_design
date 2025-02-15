from time import time
import matplotlib
import numpy as np
from antifungal.design import segment, single_point_mutate, globally_optimize

matplotlib.use("Agg")

def main():
    # Define the initial prototype peptide sequence, mesoricin1
    mesoricin1 = "RRYCRTYWRYGRLRRRCYRRRVWIWFRL"
    print(f"The prototype mesoricin1: {mesoricin1}")

    # Step 1: Segment mesoricin1 into smaller peptides and identify the one with the lowest AFI
    start_time = time()
    segment_instance = segment(mesoricin1)
    segment_predictions = segment_instance.get_candidate_sequences().predict()
    opt_mesoricin2_idx = np.argmin(segment_predictions["AFI"])
    mesoricin2 = segment_predictions["peptide_seq"][opt_mesoricin2_idx]
    print(f"The segmentation of mesoricin1 results in mesoricin2: {mesoricin2}")
    print(f"Design step 1 searched {len(segment_predictions['peptide_seq'])} sequences, taking {time() - start_time:.2f} seconds")

    # Step 2: Apply single-point mutation to mesoricin2 and find the variant with the lowest AFI
    start_time = time()
    single_point_mutate_instance = single_point_mutate(mesoricin2)
    single_point_mutate_predictions = single_point_mutate_instance.get_candidate_sequences().predict()
    opt_mesoricin3_idx = np.argmin(single_point_mutate_predictions["AFI"])
    mesoricin3 = single_point_mutate_predictions["peptide_seq"][opt_mesoricin3_idx]
    print(f"The single-point mutation of mesoricin2 results in mesoricin3: {mesoricin3}")
    print(f"Design step 2 searched {len(single_point_mutate_predictions['peptide_seq'])} sequences, taking {time() - start_time:.2f} seconds")

    # Step 3: Perform global optimization on mesoricin3 to further minimize AFI
    start_time = time()
    globally_optimize_instance = globally_optimize(mesoricin3)
    mesoricin4, _ = globally_optimize_instance.optimize()
    print(f"The global optimization of mesoricin3 results in mesoricin4: {mesoricin4}")
    print(f"Design step 3 took {time() - start_time:.2f} seconds")

if __name__ == "__main__":
    main()
