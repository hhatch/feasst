import argparse
import feasst as fst

print(fst.version())
parser = argparse.ArgumentParser()
parser.add_argument("--seed", type=str, help="random number generator seed",
    default="time")
parser.add_argument("--length", type=float, help="cubic box length", default=8)
parser.add_argument("--num", type=int, help="number of particles", default=50)
parser.add_argument("--data", type=str, help="LMP forcefield data file",
    default=fst.install_dir() + "/forcefield/data.lj")
parser.add_argument("--beta", type=float, help="inverse temperature",
    default=1.2)
parser.add_argument("--equilibration", type=int, help="number of equilibration trials",
    default=int(1e6))
parser.add_argument("--production", type=int, help="number of production trials",
    default=int(1e6))
args = parser.parse_args()
print("args:", args)

mc = fst.MakeMonteCarlo(fst.arglist([
    ["RandomMT19937", {"seed": args.seed}],
    ["Configuration", {"cubic_box_length": str(args.length),
                      "particle_type": args.data}],
    ["Potential", {"Model": "LennardJones"}],
    ["Potential", {"VisitModel": "LongRangeCorrections"}],
    ["ThermoParams", {"beta": "0.1", "chemical_potential": "10"}],
    ["Metropolis", {}],
    ["TrialTranslate", {"tunable_param": "2.", "tunable_target_acceptance": "0.2"}],
    ["TrialAdd", {"particle_type": "0"}],
    ["Run", {"until_num_particles": str(args.num)}],
    ["ThermoParams", {"beta": str(args.beta)}],
    ["RemoveTrial", {"index": "1"}],
    ["Tune", {"steps_per": str(1e5)}],
    ["CheckEnergy", {"steps_per": str(1e5), "tolerance": str(1e-8)}],
    ["Run", {"num_attempts": str(args.equilibration)}],
    ["Log", {"steps_per": str(1e5), "file_name": "lj.txt"}],
    ["Movie", {"steps_per": str(1e5), "file_name": "lj.xyz"}],
    ["Run", {"num_attempts": str(args.production)}],
]))
