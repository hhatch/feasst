#include <iostream>
#include "feasst.h"

static feasst::ArgumentParse args("A canonical ensemble Metropolis Monte Carlo simulation of a bulk Lennard Jones fluid.", {
  {"--seed", "random number generator seed", "time"},
  {"--length", "cubic periodic boundary length", "8"},
  {"--num", "number of particles", "50"},
  {"--data", "LMP forcefield data file",
    feasst::install_dir() + "/forcefield/data.lj"},
  {"--beta", "inverse temperature", "1.2"},
  {"--steps_per", "number of trials per analysis or check", feasst::str(1e5)},
  {"--equilibrate", "number of trials for equilibration", feasst::str(1e6)},
  {"--production", "number of trials for production", feasst::str(1e6)}
});

int main(int argc, char ** argv) {
  std::cout << feasst::version() << std::endl
            << args.parse(argc, argv) << std::endl;
  feasst::MakeMonteCarlo({{
    {"RandomMT19937", {{"seed", args.get("--seed")}}},
    {"Configuration", {{"cubic_box_length", args.get("--length")},
                       {"particle_type", args.get("--data")}}},
    {"Potential", {{"Model", "LennardJones"}}},
    {"Potential", {{"VisitModel", "LongRangeCorrections"}}},
    {"ThermoParams", {{"beta", args.get("--beta")},
                      {"chemical_potential", "10"}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"tunable_param", "2."},
                        {"tunable_target_acceptance", "0.2"}}},
    {"TrialAdd", {{"particle_type", "0"}}},
    {"Run", {{"until_num_particles", "50"}}},
    {"RemoveTrial", {{"name", "TrialAdd"}}},
    {"Tune", {{"steps_per", args.get("--steps_per")}}},
    {"CheckEnergy", {{"steps_per", args.get("--steps_per")}, {"tolerance", "1e-8"}}},
    {"Run", {{"num_attempts", args.get("--equilibrate")}}},
    {"Log", {{"steps_per", args.get("--steps_per")}, {"file_name", "lj.txt"}}},
    {"Movie", {{"steps_per", args.get("--steps_per")}, {"file_name", "lj.xyz"}}},
    {"Run", {{"num_attempts", args.get("--production")}}},
  }});
}
