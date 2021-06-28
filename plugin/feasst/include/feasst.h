#include "configuration/include/properties.h"
#include "configuration/include/typed_entity.h"
#include "configuration/include/bond.h"
#include "monte_carlo/include/tunable.h"
#include "system/include/model.h"
#include "system/include/synchronize_data.h"
#include "threads/include/thread.h"
#include "threads/include/thread_omp.h"
#include "utils/include/timer.h"
#include "utils/include/file.h"
#include "utils/include/argument_parse.h"
#include "utils/include/io.h"
#include "utils/include/utils.h"
#include "utils/include/custom_exception.h"
#include "utils/include/debug.h"
#include "utils/include/arguments.h"
#include "utils/include/checkpoint.h"
#include "utils/include/serialize.h"
#include "system/include/energy_map.h"
#include "system/include/visit_model_inner.h"
#include "system/include/thermo_params.h"
#include "monte_carlo/include/action.h"
#include "monte_carlo/include/run.h"
#include "configuration/include/physical_constants.h"
#include "utils/include/progress_report.h"
#include "utils/include/cache.h"
#include "math/include/table.h"
#include "math/include/accumulator.h"
#include "monte_carlo/include/stepper.h"
#include "math/include/histogram.h"
#include "math/include/formula.h"
#include "shape/include/formula_sine_wave.h"
#include "math/include/formula_polynomial.h"
#include "math/include/utils_math.h"
#include "math/include/position.h"
#include "shape/include/shape.h"
#include "shape/include/half_space_tilted.h"
#include "shape/include/sphere.h"
#include "shape/include/half_space.h"
#include "shape/include/half_space_sine.h"
#include "shape/include/cuboid.h"
#include "shape/include/cylinder.h"
#include "shape/include/shape_union.h"
#include "shape/include/shape_intersect.h"
#include "shape/include/slab_sine.h"
#include "shape/include/finite_cylinder.h"
#include "shape/include/slab.h"
#include "system/include/bond_three_body.h"
#include "system/include/angle_square_well.h"
#include "system/include/bond_two_body.h"
#include "system/include/bond_visitor.h"
#include "system/include/bond_square_well.h"
#include "system/include/neighbor_criteria.h"
#include "cluster/include/energy_map_all.h"
#include "cluster/include/energy_map_neighbor.h"
#include "cluster/include/energy_map_neighbor_criteria.h"
#include "cluster/include/energy_map_all_criteria.h"
#include "system/include/visit_model.h"
#include "system/include/model_two_body.h"
#include "system/include/lennard_jones.h"
#include "models/include/yukawa.h"
#include "models/include/square_well.h"
#include "example/include/model_example.h"
#include "system/include/model_two_body_factory.h"
#include "system/include/hard_sphere.h"
#include "system/include/ideal_gas.h"
#include "system/include/dont_visit_model.h"
#include "system/include/model_one_body.h"
#include "system/include/model_empty.h"
#include "system/include/model_three_body.h"
#include "system/include/visit_model_intra_map.h"
#include "system/include/visit_model_intra.h"
#include "system/include/visit_model_bond.h"
#include "system/include/long_range_corrections.h"
#include "configuration/include/site.h"
#include "configuration/include/particle.h"
#include "configuration/include/model_params.h"
#include "models/include/lennard_jones_alpha.h"
#include "models/include/lennard_jones_cut_shift.h"
#include "models/include/lennard_jones_force_shift.h"
#include "models/include/mie.h"
#include "system/include/model_two_body_table.h"
#include "configuration/include/group.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/select.h"
#include "system/include/cells.h"
#include "system/include/visit_model_cell.h"
#include "monte_carlo/include/rosenbluth.h"
#include "monte_carlo/include/acceptance.h"
#include "configuration/include/visit_particles.h"
#include "configuration/include/configuration.h"
#include "system/include/potential.h"
#include "system/include/potential_factory.h"
#include "system/include/system.h"
#include "system/include/utils.h"
#include "monte_carlo/include/trial_select.h"
#include "chain/include/select_site_of_type.h"
#include "chain/include/select_perturbed.h"
#include "beta_expanded/include/select_nothing.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "chain/include/select_segment.h"
#include "chain/include/select_end_segment.h"
#include "chain/include/select_reptate.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/select_particle_avb_divalent.h"
#include "cluster/include/select_cluster.h"
#include "monte_carlo/include/trial_select_bond.h"
#include "monte_carlo/include/trial_select_angle.h"
#include "chain/include/select_branch.h"
#include "monte_carlo/include/perturb.h"
#include "monte_carlo/include/perturb_volume.h"
#include "beta_expanded/include/perturb_beta.h"
#include "monte_carlo/include/perturb_move.h"
#include "monte_carlo/include/perturb_distance.h"
#include "chain/include/perturb_reptate.h"
#include "monte_carlo/include/perturb_translate.h"
#include "cluster/include/perturb_point_reflect.h"
#include "morph/include/perturb_particle_type.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_factory.h"
#include "monte_carlo/include/analyze.h"
#include "chain/include/check_rigid_bonds.h"
#include "steppers/include/density_profile.h"
#include "chain/include/analyze_bonds.h"
#include "monte_carlo/include/analyze_factory.h"
#include "steppers/include/log_and_movie.h"
#include "steppers/include/seek_analyze.h"
#include "steppers/include/wall_clock_limit.h"
#include "steppers/include/profile_trials.h"
#include "steppers/include/volume.h"
#include "steppers/include/chirality_2d.h"
#include "steppers/include/check_physicality.h"
#include "steppers/include/check.h"
#include "steppers/include/extensive_moments.h"
#include "steppers/include/cpu_time.h"
#include "steppers/include/mean_squared_displacement.h"
#include "steppers/include/num_particles.h"
#include "steppers/include/log.h"
#include "steppers/include/energy.h"
#include "steppers/include/criteria_writer.h"
#include "monte_carlo/include/modify.h"
#include "steppers/include/pair_distribution.h"
#include "steppers/include/check_energy.h"
#include "monte_carlo/include/modify_factory.h"
#include "monte_carlo/include/monte_carlo.h"
#include "prefetch/include/prefetch.h"
#include "steppers/include/check_energy_and_tune.h"
#include "steppers/include/increment_phase.h"
#include "steppers/include/seek_modify.h"
#include "steppers/include/check_properties.h"
#include "steppers/include/criteria_updater.h"
#include "steppers/include/tune.h"
#include "monte_carlo/include/seek_num_particles.h"
#include "cluster/include/trial_avb2.h"
#include "cluster/include/trial_rigid_cluster.h"
#include "cluster/include/trial_transfer_avb_divalent.h"
#include "cluster/include/trial_transfer_avb.h"
#include "chain/include/trials.h"
#include "chain/include/trial_grow.h"
#include "morph/include/trial_morph.h"
#include "morph/include/trial_morph_expanded.h"
#include "cluster/include/trial_avb4.h"
#include "beta_expanded/include/trial_beta.h"
#include "monte_carlo/include/trial_compute_volume.h"
#include "monte_carlo/include/trial_compute_remove.h"
#include "monte_carlo/include/trial_compute_add.h"
#include "monte_carlo/include/trial_compute_move.h"
#include "monte_carlo/include/trial_move.h"
#include "monte_carlo/include/trials.h"
#include "beta_expanded/include/compute_beta.h"
#include "cluster/include/compute_move_cluster.h"
#include "cluster/include/compute_add_avb.h"
#include "cluster/include/compute_remove_avb.h"
#include "cluster/include/compute_gca.h"
#include "cluster/include/compute_add_avb_divalent.h"
#include "morph/include/compute_morph.h"
#include "cluster/include/compute_remove_avb_divalent.h"
#include "cluster/include/compute_avb2.h"
#include "cluster/include/compute_avb4.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/constraint.h"
#include "monte_carlo/include/constrain_num_particles.h"
#include "egce/include/a_equal_b.h"
#include "egce/include/a_half_b.h"
#include "configuration/include/visit_configuration.h"
#include "configuration/include/file_xyz.h"
#include "steppers/include/movie.h"
#include "configuration/include/utils.h"
#include "configuration/include/file_lmp.h"
#include "configuration/include/domain.h"
#include "math/include/random.h"
#include "math/include/random_modulo.h"
#include "math/include/constants.h"
#include "math/include/minimize.h"
#include "math/include/golden_search.h"
#include "math/include/quadratic_equation.h"
#include "math/include/formula_exponential.h"
#include "math/include/matrix.h"
#include "monte_carlo/include/perturb_rotate.h"
#include "chain/include/perturb_pivot.h"
#include "chain/include/perturb_crankshaft.h"
#include "cluster/include/perturb_rotate_com.h"
#include "cluster/include/perturb_move_avb.h"
#include "cluster/include/perturb_add_avb.h"
#include "monte_carlo/include/perturb_anywhere.h"
#include "chain/include/perturb_site_type.h"
#include "chain/include/trial_grow_linear.h"
#include "monte_carlo/include/perturb_add.h"
#include "monte_carlo/include/perturb_remove.h"
#include "monte_carlo/include/perturb_distance_angle.h"
#include "chain/include/perturb_branch.h"
#include "math/include/random_mt19937.h"
#include "math/include/solver.h"
#include "math/include/solver_newton_raphson.h"
#include "math/include/solver_bisection.h"
#include "math/include/solver_brent_dekker.h"
#include "opt_lj/include/visit_model_opt_lj.h"
#include "opt_lj/include/visit_model_opt_rpm.h"
#include "ewald/include/slab_correction.h"
#include "ewald/include/trial_remove_multiple.h"
#include "ewald/include/trial_add_multiple.h"
#include "ewald/include/charge_self.h"
#include "ewald/include/charge_screened.h"
#include "ewald/include/utils.h"
#include "ewald/include/trial_transfer_multiple.h"
#include "ewald/include/coulomb.h"
#include "ewald/include/compute_remove_multiple.h"
#include "ewald/include/electric_field.h"
#include "ewald/include/compute_add_multiple.h"
#include "ewald/include/charge_screened_intra.h"
#include "ewald/include/ewald.h"
#include "ewald/include/check_net_charge.h"
#include "mayer/include/criteria_mayer.h"
#include "confinement/include/model_lj_shape.h"
#include "confinement/include/trial_anywhere.h"
#include "confinement/include/model_square_well_shape.h"
#include "confinement/include/always_reject.h"
#include "confinement/include/henry_coefficient.h"
#include "confinement/include/model_hard_shape.h"
#include "confinement/include/model_table_cartesian.h"
#include "patch/include/patch_angle.h"
#include "patch/include/visit_model_inner_patch.h"
#include "flat_histogram/include/window.h"
#include "flat_histogram/include/window_custom.h"
#include "flat_histogram/include/window_exponential.h"
#include "flat_histogram/include/ln_probability.h"
#include "flat_histogram/include/collection_matrix.h"
#include "flat_histogram/include/bias.h"
#include "flat_histogram/include/transition_matrix.h"
#include "flat_histogram/include/wang_landau.h"
#include "flat_histogram/include/wltm.h"
#include "flat_histogram/include/macrostate.h"
#include "beta_expanded/include/macrostate_beta.h"
#include "morph/include/macrostate_morph.h"
#include "flat_histogram/include/macrostate_num_particles.h"
#include "flat_histogram/include/ensemble.h"
#include "flat_histogram/include/flat_histogram.h"
#include "flat_histogram/include/clones.h"
#include "flat_histogram/include/macrostate_energy.h"
std::shared_ptr<feasst::ComputeBeta> __feasst__ComputeBeta = std::make_shared<feasst::ComputeBeta>();
std::shared_ptr<feasst::Tune> __feasst__Tune = std::make_shared<feasst::Tune>();
