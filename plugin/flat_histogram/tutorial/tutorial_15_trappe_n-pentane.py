import argparse
import feasst as fst

print(fst.version())
parser = argparse.ArgumentParser()
parser.add_argument("--task", type=int, help="SLURM job array index", default=0)
parser.add_argument("--num_procs", type=int, help="number of processors", default=12)
parser.add_argument("--num_hours", type=float, help="number of hours before restart", default=5*24)
parser.add_argument("--dccb_cutoff", type=int, help="cutoff for dual cut configurational bias", default=4)
parser.add_argument("--dccb_begin", type=int, help="number of molecules before using DCCB", default=300)
parser.add_argument("--lx", type=float, help="box length in x", default=33.0)
parser.add_argument("--ly", type=float, help="box length in y", default=33.0)
parser.add_argument("--lz", type=float, help="box length in z", default=33.0)
parser.add_argument("--max_particles", type=int, help="maximum number of particles", default=400)
parser.add_argument("--temperature", type=float, help="temperature in Kelvin", default=263.15)
parser.add_argument("--particle", "-p", type=str, help="data file of alkane", required=True)
parser.add_argument("--collect_flatness", type=int, help="number of WL flatness to begin collection", default=18)
parser.add_argument("--min_flatness", type=int, help="number of WL flatness to switch to TM", default=22)
parser.add_argument("--beta_mu", type=int, help="baseline chemical potential of each species", default=-7)
args = parser.parse_args()
print("args:", args)

def mc(thread, mn, mx):
    steps_per=int(1e4)
    mc = fst.MakeMonteCarlo()
    mc.add(fst.MakeConfiguration(fst.args({"side_length0": str(args.lx), "side_length1": str(args.ly), "side_length2": str(args.lz),
        "particle_type0": args.particle})))
    mc.add(fst.MakePotential(fst.MakeLennardJones()))
    mc.add(fst.MakePotential(fst.MakeLennardJones(),
                             fst.MakeVisitModelIntra(fst.args({"cutoff": "4"}))))
    mc.add(fst.MakePotential(fst.MakeBondVisitor()))
    mc.add(fst.MakePotential(fst.MakeLongRangeCorrections()))
    if mx > args.dccb_begin:
        reference = fst.Potential(fst.MakeLennardJones(), fst.MakeVisitModelCell(fst.args({"min_length": str(args.dccb_cutoff)})))
        reference.set_model_params(mc.configuration())
        for site_type in range(mc.configuration().num_site_types()):
            reference.set_model_param("cutoff", site_type, args.dccb_cutoff)
        mc.add_to_reference(reference)
        mc.add_to_reference(fst.MakePotential(fst.MakeLennardJones(),
                            fst.MakeVisitModelIntra(fst.args({"cutoff": "4"}))))
        ref = "0"
        num = "4"
    else:
        ref = "-1"
        num = "1"
    mc.set(fst.MakeThermoParams(fst.args({"beta": str(1./args.temperature), "chemical_potential0": str(args.beta_mu*args.temperature)})))
    mc.set(fst.MakeFlatHistogram(
        fst.MakeMacrostateNumParticles(
            fst.Histogram(fst.args({"width": "1", "max": str(mx), "min": str(mn)}))),
        # fst.MakeTransitionMatrix(fst.args({"min_sweeps": str(args.sweeps)})),
        fst.MakeWLTM(fst.args({
            "collect_flatness": str(args.collect_flatness),
            "min_flatness": str(args.min_flatness),
            "min_sweeps": "1000"}))))
    #mc.add(fst.MakeTrialTranslate(fst.args({"particle_type": "0", "weight": "1.",
    #    "tunable_param": "1."})))
    #mc.add(fst.MakeTrialRotate(fst.args({"particle_type": "0", "weight": "1.",
    #    "tunable_param": "1."})))
    #mc.add(fst.MakeTrialTransfer(fst.args({"particle_type": "0", "weight": "4"})))
    mc.add(fst.MakeTrialGrow(fst.ArgsVector([
        {"transfer": "true", "particle_type": "0", "site": "0", "weight": "4"},
        {"bond": "true", "mobile_site": "1", "anchor_site": "0"},
        {"angle": "true", "mobile_site": "2", "anchor_site": "1", "anchor_site2": "0"},
        {"dihedral": "true", "mobile_site": "3", "anchor_site": "2", "anchor_site2": "1", "anchor_site3": "0"},
        {"dihedral": "true", "mobile_site": "4", "anchor_site": "3", "anchor_site2": "2", "anchor_site3": "1"},
    ]), fst.args({"reference_index": ref, "num_steps": num})))
    mc.add(fst.MakeTrialGrow(fst.ArgsVector([
        {"particle_type": "0", "weight": "2",
         "bond": "true", "mobile_site": "1", "anchor_site": "0"},
        {"angle": "true", "mobile_site": "2", "anchor_site": "1", "anchor_site2": "0"},
        {"dihedral": "true", "mobile_site": "3", "anchor_site": "2", "anchor_site2": "1", "anchor_site3": "0"},
        {"dihedral": "true", "mobile_site": "4", "anchor_site": "3", "anchor_site2": "2", "anchor_site3": "1"},
    ]), fst.args({"reference_index": ref, "num_steps": num})))
    mc.add(fst.MakeTrialGrow(fst.ArgsVector([
        {"particle_type": "0", "weight": "1",
         "angle": "true", "mobile_site": "2", "anchor_site": "1", "anchor_site2": "0"},
        {"dihedral": "true", "mobile_site": "3", "anchor_site": "2", "anchor_site2": "1", "anchor_site3": "0"},
        {"dihedral": "true", "mobile_site": "4", "anchor_site": "3", "anchor_site2": "2", "anchor_site3": "1"},
    ]), fst.args({"reference_index": ref, "num_steps": num})))
    mc.add(fst.MakeTrialGrow(fst.ArgsVector([
        {"particle_type": "0", "weight": "0.5",
         "dihedral": "true", "mobile_site": "3", "anchor_site": "2", "anchor_site2": "1", "anchor_site3": "0"},
        {"dihedral": "true", "mobile_site": "4", "anchor_site": "3", "anchor_site2": "2", "anchor_site3": "1"},
    ]), fst.args({"reference_index": ref, "num_steps": num})))
    mc.add(fst.MakeTrialGrow(fst.ArgsVector([
        {"particle_type": "0", "weight": "0.25",
         "dihedral": "true", "mobile_site": "4", "anchor_site": "3", "anchor_site2": "2", "anchor_site3": "1"},
    ]), fst.args({"reference_index": ref, "num_steps": num})))

    mc.add(fst.MakeCheckEnergy(fst.args({"steps_per": str(steps_per), "tolerance": "0.0001"})))
    mc.add(fst.MakeTune(fst.args({"steps_per": str(steps_per), "stop_after_phase": "0"})))
    mc.add(fst.MakeLogAndMovie(fst.args({"steps_per": str(steps_per),
                                         "file_name": "clones" + str(thread),
                                         "file_name_append_phase": "True"})))
    mc.add(fst.MakeEnergy(fst.args({
        "file_name": "en" + str(thread) + '.txt',
        "file_name_append_phase": "True",
        "start_after_phase": "0",
        "steps_per_write": str(steps_per),
        "steps_per_update": "1",
        "multistate": "True"})))
    mc.add(fst.MakeCriteriaUpdater(fst.args({"steps_per": str(steps_per)})))
    mc.add(fst.MakeCriteriaWriter(fst.args({"steps_per": str(steps_per),
                                            "file_name": "clones" + str(thread) + "_crit.txt",
                                            "file_name_append_phase": "True"})))
    mc.set(fst.MakeCheckpoint(fst.args({"file_name": "checkpoint" + str(thread) + ".fst",
                                        "num_hours_terminate": str(0.9*args.num_procs*args.num_hours)})))
    return mc

windows=fst.WindowExponential(fst.args({
    "alpha": "1.75",
    "num": str(args.num_procs),
    "maximum": str(args.max_particles),
    "extra_overlap": "2"})).boundaries()
print(windows)

if args.task == 0:
    clones = fst.MakeClones()
    for proc, win in enumerate(windows):
        clones.add(mc(proc, win[0], win[1]))
    clones.set(fst.MakeCheckpoint(fst.args({"file_name": "checkpoint.fst"})))
else:
    clones = fst.MakeClones("checkpoint", args.num_procs)
#clones.initialize_and_run_until_complete()
clones.initialize_and_run_until_complete(fst.args({"ln_prob_file": "ln_prob.txt"}))
print(clones.ln_prob().values())
open('clones.fst', 'w').write(clones.serialize())
