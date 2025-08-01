import ase
from ase.io import read
import time
import numpy as np
import glob


def CCC(y_true, y_pred):
    '''Concordance Correlation Coefficient (CCC)'''


    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)

    mean_true = np.mean(y_true)
    mean_pred = np.mean(y_pred)

    std_true = np.std(y_true, ddof=1)
    std_pred = np.std(y_pred, ddof=1)

    if std_true == 0 and std_pred == 0:  # Handles identical constant lists
        return -5  # Defined behavior when both lists are constant

    pearson_corr = np.corrcoef(y_true, y_pred)[0, 1]

    if np.isnan(pearson_corr):  # If correlation fails (e.g., due to NaN values)
        return -5  

    ccc = (2 * pearson_corr * std_true * std_pred) / (
        std_true**2 + std_pred**2 + (mean_true - mean_pred) ** 2
    )

    return ccc if not np.isnan(ccc) else -5  # Ensure no NaN result


def parallel_function_call(func, arg):
    start_time = time.time()
    result = func(arg)
    end_time = time.time()
    duration = end_time - start_time
    
    # Log execution time
    with open("execution_times.log", "a") as log_file:
        log_file.write(f"{func.__name__} took {duration:.4f} seconds\n")
    
    return result

def read_lowest_energy_frame(filename):
    frames = read(filename, index=":")
    energy_values = [(i, atoms.get_potential_energy()/len(atoms)) for i, atoms in enumerate(frames)]
    min_energy_frame_idx, min_energy = min(energy_values, key=lambda x: x[1])
    print(f"Lowest energy frame: {min_energy_frame_idx} with energy {min_energy:.6f}")
    return frames[min_energy_frame_idx]



def get_lowest_score_id():
    lowest_uuid = None
    lowest_score = float("inf")
    best_iteration = None

    # Iterate over all optimization log files in iteration directories
    for file_path in glob.glob("iteration_*/optimization_log.txt"):
        iteration = int(file_path.split("_")[1].split("/")[0])  # Extract iteration number
        
        with open(file_path, 'r') as file:
            data = file.read()

        for entry in filter(None, map(str.strip, data.split("-----------------------------------------"))):
            lines = [line.strip() for line in entry.split("\n") if line.strip()]
            try:
                uid = lines[0].split(": ")[1]  # Extract UUID
                score = float(next(line.split(": ")[1] for line in lines if "Total Score" in line))
                
                if score < lowest_score:
                    lowest_uuid, lowest_score, best_iteration = uid, score, iteration

            except (IndexError, ValueError, StopIteration):
                pass  # Ignore malformed entries

    return lowest_uuid, best_iteration
