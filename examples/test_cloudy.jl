### Buildkite job `single_column_precipitation_test`

using Revise; import ClimaAtmos as CA;

config_dict = Dict();
config_dict["initial_condition"] = "PrecipitatingColumnNM";
config_dict["implicit_diffusion"] = true;
config_dict["approximate_linear_solve_iters"] = 2;
config_dict["z_elem"] = 6;
config_dict["dt"] = "10secs";
config_dict["dt_cloud_fraction"] = "10secs";
config_dict["surface_setup"] = "DefaultMoninObukhov";
config_dict["t_end"] = "10secs";
config_dict["z_stretch"] = false;
config_dict["vert_diff"] = "FriersonDiffusion";
config_dict["dt_save_state_to_disk"] = "10secs";
config_dict["check_precipitation"] = false;
config_dict["config"] = "column";
config_dict["z_max"] = 250.0;
config_dict["precip_model"] = "NM";
config_dict["regression_test"] = false;
config_dict["toml"] = ["toml/single_column_precipitation_test.toml"];
config_dict["diagnostics"] = Dict{Any, Any}[Dict("short_name" => ["hus", "clw", "cli", "husra", "ta", "wa"], "period" => "500secs")];
config_dict["job_id"] = "single_column_cloudy_test";
config_dict["moist"] = "nonequil";

config = CA.AtmosConfig(config_dict);

include("hybrid/driver.jl")