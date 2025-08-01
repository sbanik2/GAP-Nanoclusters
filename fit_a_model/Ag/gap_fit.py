import numpy as np
import subprocess
from calculator import PropertyCalculator
from utilities import read_lowest_energy_frame



class GAPModelFitter:
    def __init__(self, modalities, training_data, test_data_file, evaluation_files):
        self.modalities = modalities
        self.training_data = training_data
        self.test_data_file = test_data_file
        self.evaluation_files = evaluation_files


    def fit_gap_model(self, hyperparams):
        cutoff = hyperparams[0]
        sigma_energy = {self.modalities[i]: hyperparams[i + 1] for i in range(len(self.modalities))}
        force_factor = hyperparams[1 + len(self.modalities)]
        sigma_force = {modality: sigma_energy[modality] * force_factor for modality in self.modalities}
        
        gp_filename = "gap_model.xml"
        gap_params = (
            f"gap={{ soap l_max=8 n_max=8 atom_sigma=0.5 cutoff={cutoff} "
            "radial_scaling=-0.5 cutoff_transition_width=1.0 central_weight=1.0 "
            "n_sparse=200 delta=0.2 covariance_type=dot_product zeta=4 sparse_method=cur_points }} "
            f"default_sigma={{ {sigma_energy[self.modalities[0]]} {sigma_force[self.modalities[0]]} 0.2 0.0 }} "
            "config_type_sigma={ " +
            " ".join([f"{mod}:{sigma_energy[mod]}:{sigma_force[mod]}:0.2:0.0" for mod in self.modalities]) + " } "
            "e0=-0.198095 energy_parameter_name=energy force_parameter_name=forces "
            "virial_parameter_name=stress sparse_jitter=1.0e-8 do_copy_at_file=F sparse_separate_file=F "
            f"gp_file={gp_filename}"
        )
        
        cmd = f"gap_fit atoms_filename={self.training_data} {gap_params}"
        subprocess.run(cmd, shell=True, check=True)


    def objective_function(self, objectives, priorities):
        gp_filename = "gap_model.xml"
        calculator = PropertyCalculator(
            self.evaluation_files['ground_state']
        )
        
        delE_cohe, ordering_coheE, del_Latt = calculator.compute_lattice_errors(gp_filename, save_file='lattice_GS.json')
        ccc_errors = calculator.ef_test_set(gp_filename, 
                                            self.test_data_file,
                                            self.modalities,
                                            save_e="E_pred.json",
                                            save_f="F_pred.json")
        delE = sum([1 - x for x in list(ccc_errors["energy"].values())])
        delF = sum([1 - x for x in list(ccc_errors["force"].values())])

        result = {
            'ordering_coheE': ordering_coheE,
            'delF': delF,
            'del_Latt': del_Latt,
            'delE_cohe': delE_cohe,
            'delE': delE,
        }

        evaluation = np.array([result[key] for key in objectives])
        weights = np.array([priorities[key] for key in objectives])
        
        score = np.sum(evaluation * weights)
        
        self.log_optimization_step(result, score)

    def log_optimization_step(self, result, score):
        log_file = "optimization_log.txt"
        with open(log_file, "a") as f:
            f.write("Optimization Results:\n")
            for key, value in result.items():
                f.write(f"{key}: {value}\n")
            f.write(f"Total Score: {score}\n")
            f.write("-----------------------------------------\n")
        
        print(f"Logged optimization results in {log_file}")


# Main Script

modalities = ['cluster', 'bulk']
training_data = "train.xyz"
test_data_file = "test.xyz"

evaluation_files = {
    'ground_state': 'ground_state.xyz', 
}


ref_atoms = read_lowest_energy_frame('ground_state.xyz')

# Objective Function Setup
priorities = {'ordering_coheE': 2, 'del_SE': 2, 'delF': 2, 'del_Latt': 1, 'delE_cohe': 1, 'delE': 1}
objectives = ['ordering_coheE', 'delF', 'del_Latt', 'delE_cohe', 'delE']

# Initialize GAP Model Fitter
gap_fitter = GAPModelFitter(modalities, training_data, test_data_file, evaluation_files)
hyperparams = [float(x) for x in open("param.txt").read().strip().split(",")]
gap_fitter.fit_gap_model(hyperparams)
gap_fitter.objective_function(objectives, priorities)
