# define a pure component NVT MC Lennard-Jones simulation

default_sim_params = {
    "seed": "time",
    "length": 8,
    "num_particles": 50,
    "fstprt": "/feasst/forcefield/lj.fstprt",
    "beta": 1.2,
    "steps_per": 1e5,
    "equilibration": 1e6,
    "production": 1e6,
    "sim": 0,
}

def mc_lj(params=default_sim_params, file_name="tutorial.txt"):
    with open(file_name, "w") as myfile: myfile.write("""
# high temperature gcmc to generate initial configuration
RandomMT19937 seed {seed}
Configuration cubic_box_length {length} particle_type0 {fstprt}
Potential Model LennardJones
Potential VisitModel LongRangeCorrections
ThermoParams beta 0.1 chemical_potential 10
Metropolis
TrialTranslate tunable_param 2 tunable_target_acceptance 0.2
TrialAdd particle_type 0
Run until_num_particles {num_particles}

# nvt equilibration
RemoveTrial name TrialAdd
ThermoParams beta {beta}
Tune steps_per {steps_per}
CheckEnergy steps_per {steps_per} tolerance 1e-8
Run num_attempts {equilibration}

# nvt production
RemoveModify name Tune
Log steps_per {steps_per} file_name lj{sim}.txt
Movie steps_per {steps_per} file_name lj{sim}.xyz
Energy steps_per_write {steps_per} file_name en{sim}.txt
Run num_attempts {production}
""".format(**params))

if __name__ == "__main__":
    import subprocess
    mc_lj()
    subprocess.call("~/feasst/build/bin/fst < tutorial.txt", shell=True, executable='/bin/bash')
