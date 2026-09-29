#include "neo/problems/fixed_angle.hpp"
#include "neo/problems/fixed_angle_mitm.hpp"
#include "neo/problems/obtuse.hpp"
#include "neo/problems/parabola.hpp"
#include "neo/problems/circumcircle.hpp"

#include <cstdlib>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

namespace {

void usage() {
    std::cout << "Usage: neo-euclid fixed-angle [--target N] [--max-cost N] "
                 "[--beam N] [--max-points N] [--no-macros] [--quiet]\n"
                 "       neo-euclid fixed-angle-mitm [--target N] [--max-cost N] "
                 "[--shared-cost N] [--arm-cost N] [--beam N] [--max-points N]\n"
                 "       neo-euclid obtuse [--max-cost N] [--seconds N] "
                 "[--expansions N] [--max-points N] [--seed N] [--quiet]\n"
                 "       neo-euclid obtuse-mitm [--max-cost N] [--seconds N] "
                 "[--shared-min N] [--shared-max N] [--records N] "
                 "[--branch-attempts N] [--raw] [--quiet]\n"
                 "       neo-euclid parabola-search [--max-cost N] [--seconds N] "
                 "[--threads N] [--prefix-depth N] [--samples 2|3] [--sample-degree N] [--mask LCLC] "
                 "[--beam-search] [--beam N] [--tangent] [--center-target] [--quiet]\n";
}

int integer_value(char** argv, int& index, int argc) {
    if (++index >= argc) throw std::runtime_error("missing option value");
    return std::stoi(argv[index]);
}

double double_value(char** argv, int& index, int argc) {
    if (++index >= argc) throw std::runtime_error("missing option value");
    return std::stod(argv[index]);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 2) {
            usage();
            return EXIT_SUCCESS;
        }
        if (std::string(argv[1]) == "parabola-search") {
            neo::parabola::SearchConfig config; bool tangent=false,mitm=false,beam_search=false; int samples=3; std::vector<double> explicit_samples;
            for(int i=2;i<argc;++i){const std::string option=argv[i];
                if(option=="--max-cost")config.max_cost=integer_value(argv,i,argc);
                else if(option=="--seconds")config.time_limit_seconds=double_value(argv,i,argc);
                else if(option=="--threads")config.threads=integer_value(argv,i,argc);
                else if(option=="--beam")config.beam_width=static_cast<std::size_t>(integer_value(argv,i,argc));
                else if(option=="--prefix-depth")config.prefix_depth=integer_value(argv,i,argc);
                else if(option=="--samples")samples=integer_value(argv,i,argc);
                else if(option=="--sample-degree")explicit_samples.push_back(double_value(argv,i,argc));
                else if(option=="--max-points")config.max_points=static_cast<std::size_t>(integer_value(argv,i,argc));
                else if(option=="--max-states")config.max_states=static_cast<std::size_t>(integer_value(argv,i,argc));
                else if(option=="--state-cache")config.state_cache_entries=static_cast<std::size_t>(integer_value(argv,i,argc));
                else if(option=="--geometry-cache")config.geometry_cache_entries=static_cast<std::size_t>(integer_value(argv,i,argc));
                else if(option=="--mask")config.masks.push_back(argv[++i]);
                else if(option=="--tangent")tangent=true;
                else if(option=="--center-target"){config.target_circle_center=true;config.angle_at_parabola_vertex=true;}
                else if(option=="--e-target"){config.target_e_point=true;config.angle_at_parabola_vertex=true;}
                else if(option=="--euclidea-cost")config.cost_policy=neo::engine::CostPolicy::Euclidea;
                else if(option=="--no-transfer")config.allow_transfer_circle=false;
                else if(option=="--exhaustive-311"){config.exhaustive_311=true;config.target_e_point=true;config.angle_at_parabola_vertex=true;config.cost_policy=neo::engine::CostPolicy::Euclidea;}
                else if(option=="--five-step-e-prefix"){config.retain_five_step_e_prefix=true;config.angle_at_parabola_vertex=true;config.cost_policy=neo::engine::CostPolicy::Euclidea;}
                else if(option=="--exhaustive-three-after-e"){beam_search=true;config.exhaustive_three_after_e=true;config.retain_five_step_e_prefix=true;config.allow_transfer_circle=false;config.angle_at_parabola_vertex=true;config.cost_policy=neo::engine::CostPolicy::Euclidea;config.max_cost=3;}
                else if(option=="--exhaustive-cost4-e-mitm"){mitm=true;config.exhaustive_cost4_e_mitm=true;config.target_e_point=true;config.allow_transfer_circle=false;config.angle_at_parabola_vertex=true;config.cost_policy=neo::engine::CostPolicy::Euclidea;config.max_cost=4;}
                else if(option=="--exhaustive-cost4-e-sequential"){beam_search=true;config.exhaustive_cost4_e_sequential=true;config.target_e_point=true;config.allow_transfer_circle=false;config.angle_at_parabola_vertex=true;config.cost_policy=neo::engine::CostPolicy::Euclidea;config.max_cost=4;}
                else if(option=="--require-bisector")config.require_bisector_use=true;
                else if(option=="--k2-circle"){config.target_k2_circle=true;config.angle_at_parabola_vertex=true;config.validation_max_degrees=53.;}
                else if(option=="--mitm")mitm=true;
                else if(option=="--beam-search")beam_search=true;
                else if(option=="--vertex-angle")config.angle_at_parabola_vertex=true;
                else if(option=="--structured-pruning"){config.require_first_step_uses_p=true;config.require_circle_parabola_use=true;}
                else if(option=="--quiet")config.verbose=false;
                else if(option=="--help"){usage();return EXIT_SUCCESS;}
                else throw std::runtime_error("unknown option: "+option);
            }
            if(!explicit_samples.empty())config.search_degrees=std::move(explicit_samples);
            else if(config.target_k2_circle)config.search_degrees=samples==2?std::vector<double>{21,47}:std::vector<double>{17,33,49};
            else if(config.angle_at_parabola_vertex)config.search_degrees=samples==2?std::vector<double>{23,71}:std::vector<double>{17,43,79};
            else if(samples==2)config.search_degrees={38,137};else if(samples!=3)throw std::runtime_error("samples must be 2 or 3");
            const auto report=tangent?neo::parabola::search_tangent(config):(mitm?neo::parabola::search_mitm(config):(beam_search?neo::parabola::search_beam(config):neo::parabola::search(config)));
            std::cout<<(report.found?"FOUND":"No construction found")<<" cost="<<report.cost
                     <<" verified="<<report.densely_verified<<" schema="<<report.terminal_schema
                     <<" mask="<<report.mask<<" prefixes="<<report.prefixes
                     <<" expanded="<<report.expanded<<" generated="<<report.generated
                     <<" duplicates="<<report.duplicates<<" dense-rejections="<<report.dense_rejections
                     <<" state-cache="<<report.state_cache_entries<<" state-hits="<<report.state_cache_hits
                     <<" geometry-cache="<<report.geometry_cache_entries
                     <<" l1-parabola-hits="<<report.local_parabola_cache_hits
                     <<" l1-pair-hits="<<report.local_pair_cache_hits
                     <<" l2-parabola-hits="<<report.parabola_cache_hits<<" l2-pair-hits="<<report.pair_cache_hits
                     <<" mitm-records="<<report.mitm_records<<" mitm-matches="<<report.mitm_matches
                     <<" elapsed="<<report.elapsed_seconds<<"s\n";
            for(std::size_t i=0;i<report.steps.size();++i)std::cout<<"  "<<i+1<<". "<<report.steps[i]<<'\n';
            if(!report.found&&!report.sample_only_steps.empty()){
                std::cout<<"SAMPLE-ONLY candidate schema="<<report.sample_only_schema
                         <<" mask="<<report.sample_only_mask<<'\n';
                for(std::size_t i=0;i<report.sample_only_steps.size();++i)
                    std::cout<<"  "<<i+1<<". "<<report.sample_only_steps[i]<<'\n';
            }
            return report.found?EXIT_SUCCESS:2;
        }
        if (std::string(argv[1]) == "circumcircle-mitm") {
            double seconds=1800.;std::size_t states=200000000,points=40;
            for(int i=2;i<argc;++i){const std::string option=argv[i];
                if(option=="--seconds")seconds=double_value(argv,i,argc);
                else if(option=="--max-states")states=static_cast<std::size_t>(integer_value(argv,i,argc));
                else if(option=="--max-points")points=static_cast<std::size_t>(integer_value(argv,i,argc));
                else throw std::runtime_error("unknown option: "+option);
            }
            const auto r=neo::circumcircle::search_mitm(seconds,states,points);
            std::cout<<(r.found?"FOUND":"No construction found")<<" cost="<<r.cost<<" generated="<<r.generated<<" states="<<r.states<<" target-curves="<<r.target_curves<<" elapsed="<<r.elapsed_seconds<<"s\n";
            if(!r.left.empty())std::cout<<"  left: "<<r.left<<'\n';if(!r.right.empty())std::cout<<"  right: "<<r.right<<'\n';if(r.found&& !r.right.empty())std::cout<<"  final: circle(intersection(left,right); A)\n";
            return r.found?EXIT_SUCCESS:2;
        }
        if (std::string(argv[1]) == "fixed-angle-mitm") {
            neo::problems::FixedAngleMitmConfig config;
            double target = 144.0;
            for (int i = 2; i < argc; ++i) {
                const std::string option = argv[i];
                if (option == "--target") target = double_value(argv, i, argc);
                else if (option == "--max-cost") config.max_cost = integer_value(argv, i, argc);
                else if (option == "--shared-cost") config.shared_cost = integer_value(argv, i, argc);
                else if (option == "--arm-cost") config.arm_cost = integer_value(argv, i, argc);
                else if (option == "--beam") config.prefix_beam = static_cast<std::size_t>(integer_value(argv, i, argc));
                else if (option == "--max-points") config.max_points = static_cast<std::size_t>(integer_value(argv, i, argc));
                else if (option == "--max-states") config.max_states = static_cast<std::size_t>(integer_value(argv, i, argc));
                else if (option == "--seconds") config.time_limit_seconds = double_value(argv, i, argc);
                else if (option == "--macros") config.use_macros = true;
                else if (option == "--euclidea-cost") config.cost_policy = neo::engine::CostPolicy::Euclidea;
                else if (option == "--help") { usage(); return EXIT_SUCCESS; }
                else throw std::runtime_error("unknown option: " + option);
            }
            const auto result = neo::problems::fixed_angle_mitm(target, config);
            if (!result) {
                std::cout << "No MITM construction found.\n";
                return 2;
            }
            std::cout << "FOUND by MITM at cost " << result->goal.total_cost << '\n';
            int number = 1;
            for (const auto recipe : result->state.steps)
                std::cout << "  " << number++ << ". "
                          << neo::describe_recipe(result->state, recipe) << '\n';
            std::cout << "  " << number << ". " << result->goal.final_step << '\n';
            std::cout << "expanded-prefixes=" << result->expanded
                      << " generated=" << result->generated << '\n';
            return EXIT_SUCCESS;
        }
        if (std::string(argv[1]) == "obtuse" || std::string(argv[1]) == "obtuse-mitm") {
            const bool use_mitm = std::string(argv[1]) == "obtuse-mitm";
            neo::obtuse::MeetConfig config;
            double range_min = 0.0, range_max = 0.0;
            int sample_count = 8;
            for (int i = 2; i < argc; ++i) {
                const std::string option = argv[i];
                if (option == "--max-cost") config.max_cost = integer_value(argv, i, argc);
                else if (option == "--seconds") config.time_limit_seconds = double_value(argv, i, argc);
                else if (option == "--expansions") config.expansions_per_restart = static_cast<std::size_t>(integer_value(argv, i, argc));
                else if (option == "--max-points") config.max_points = static_cast<std::size_t>(integer_value(argv, i, argc));
                else if (option == "--seed") config.seed = static_cast<std::uint64_t>(integer_value(argv, i, argc));
                else if (option == "--range-min") range_min = double_value(argv, i, argc);
                else if (option == "--range-max") range_max = double_value(argv, i, argc);
                else if (option == "--samples") sample_count = integer_value(argv, i, argc);
                else if (option == "--shared-min") config.shared_cost_min = integer_value(argv, i, argc);
                else if (option == "--shared-max") config.shared_cost_max = integer_value(argv, i, argc);
                else if (option == "--records") config.branch_records_per_base = static_cast<std::size_t>(integer_value(argv, i, argc));
                else if (option == "--branch-attempts") config.branch_attempts_per_base = static_cast<std::size_t>(integer_value(argv, i, argc));
                else if (option == "--raw") config.reflected_prefix = false;
                else if (option == "--quiet") config.verbose = false;
                else if (option == "--help") { usage(); return EXIT_SUCCESS; }
                else throw std::runtime_error("unknown option: " + option);
            }
            if (range_min != 0.0 || range_max != 0.0) {
                if (!(range_min > 90.0 && range_max < 180.0 && range_min < range_max && sample_count >= 2))
                    throw std::runtime_error("range must satisfy 90 < min < max < 180 with samples >= 2");
                for (int i=0;i<sample_count;++i)
                    config.sample_degrees.push_back(range_min+(range_max-range_min)*i/(sample_count-1));
            }
            const auto report = use_mitm ? neo::obtuse::meet_in_the_middle(config)
                                         : neo::obtuse::search(config);
            if (!report.found) {
                std::cout << "No construction found. restarts=" << report.restarts
                          << " expanded=" << report.expanded << " generated=" << report.generated
                          << " duplicates=" << report.duplicates << " elapsed=" << report.elapsed_seconds
                          << "s\n" << report.best_description << '\n';
                if (!report.best_steps.empty()) {
                    std::cout << "Best candidate steps (not yet verified as a construction):\n";
                    for (std::size_t i=0;i<report.best_steps.size();++i)
                        std::cout << "  " << i+1 << ". " << report.best_steps[i] << '\n';
                    for (const auto& definition : report.best_point_definitions)
                        std::cout << "    " << definition << '\n';
                }
                return 2;
            }
            std::cout << "FOUND obtuse construction at cost " << report.cost << '\n';
            for (std::size_t i=0;i<report.steps.size();++i)
                std::cout << "  " << i+1 << ". " << report.steps[i] << '\n';
            return EXIT_SUCCESS;
        }
        if (std::string(argv[1]) != "fixed-angle") {
            usage(); return EXIT_FAILURE;
        }
        neo::SearchConfig config;
        double target = 72.0;
        for (int i = 2; i < argc; ++i) {
            const std::string option = argv[i];
            if (option == "--target") target = double_value(argv, i, argc);
            else if (option == "--max-cost") config.max_cost = integer_value(argv, i, argc);
            else if (option == "--beam") config.beam_width = static_cast<std::size_t>(integer_value(argv, i, argc));
            else if (option == "--max-points") config.max_points = static_cast<std::size_t>(integer_value(argv, i, argc));
            else if (option == "--no-macros") config.use_macros = false;
            else if (option == "--euclidea-cost") config.cost_policy = neo::engine::CostPolicy::Euclidea;
            else if (option == "--quiet") config.verbose = false;
            else if (option == "--help") { usage(); return EXIT_SUCCESS; }
            else throw std::runtime_error("unknown option: " + option);
        }
        neo::problems::FixedAngle problem(target);
        const auto result = neo::beam_search(problem, config);
        if (!result) {
            std::cout << "No construction found inside this heuristic beam.\n";
            return 2;
        }
        std::cout << "FOUND at cost " << result->goal.total_cost << "\n";
        int number = 1;
        for (const auto recipe : result->state.steps)
            std::cout << "  " << number++ << ". " << neo::describe_recipe(result->state, recipe) << '\n';
        if (result->goal.total_cost > result->state.cost)
            std::cout << "  " << number << ". " << result->goal.final_step << '\n';
        std::cout << "expanded=" << result->expanded << " generated=" << result->generated << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
