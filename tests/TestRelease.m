classdef TestRelease < matlab.unittest.TestCase
    % End-to-end checks for the supported v0.1.0 public workflow.

    properties
        RepoDir string
        WorkDir string
        InputDir string
        OutputDir string
    end

    methods (TestClassSetup)
        function buildFreshOutput(testCase)
            testCase.RepoDir = string(fileparts(fileparts(mfilename("fullpath"))));
            testCase.WorkDir = string(tempname);
            testCase.InputDir = fullfile(testCase.WorkDir, "inputs");
            testCase.OutputDir = fullfile(testCase.WorkDir, "output");
            mkdir(testCase.WorkDir);
            mkdir(testCase.InputDir);
            example_inputs = fullfile(testCase.RepoDir, "examples", "synthetic", "inputs");
            copyfile(fullfile(example_inputs, "CatchmentA"), ...
                fullfile(testCase.InputDir, "CatchmentA"));
            copyfile(fullfile(example_inputs, "CatchmentC"), ...
                fullfile(testCase.InputDir, "CatchmentC"));
            addpath(testCase.RepoDir);

            K_map = containers.Map({'CatchmentC'}, {3});
            run_detrital_pipeline(testCase.InputDir, testCase.OutputDir, ...
                Kmax=5, ...
                K_override_map=K_map, ...
                NSigma=1, ...
                P_thresh=0.65, ...
                Delta=8, ...
                Nmc=2, ...
                run_sensitivity=true);
        end
    end

    methods (TestClassTeardown)
        function removeFreshOutput(testCase)
            rmpath(testCase.RepoDir);
            if isfolder(testCase.WorkDir)
                rmdir(testCase.WorkDir, "s");
            end
        end
    end

    methods (Test)
        function versionIsReleaseCandidate(testCase)
            testCase.verifyEqual(detrital_chron_filter_version(), "0.1.0-rc1");
        end

        function returnedFieldNamesMatchPublicTerminology(testCase)
            catchment_dir = fullfile(testCase.InputDir, "CatchmentA");
            result = filter_detrital_thermo( ...
                fullfile(catchment_dir, "HblAr.csv"), ...
                fullfile(catchment_dir, "ApHeApPb.csv"), ...
                fullfile(catchment_dir, "ZrnHeZrnPb.csv"), ...
                fullfile(testCase.WorkDir, "direct_filter_check"), ...
                [85.4 94.7], WriteOutputs=false);
            testCase.verifyTrue(isfield(result, "filter_results"));
            testCase.verifyTrue(isfield(result, "output_summary"));
            testCase.verifyFalse(isfield(result, "screening_results"));
            testCase.verifyFalse(isfield(result, "screening_summary"));
        end

        function runLevelFilesExist(testCase)
            names = ["README.txt", "pipeline_summary.csv", "output_summary.csv", ...
                "filter_code_lookup.csv", ...
                "reference_boundary_sensitivity_summary.csv"];
            for name = names
                testCase.verifyTrue(isfile(fullfile(testCase.OutputDir, name)), ...
                    "Missing run-level output: " + name);
            end

            pipeline = readtable(fullfile(testCase.OutputDir, "pipeline_summary.csv"), ...
                Delimiter=",", VariableNamingRule="preserve", TextType="string");
            c = pipeline(string(pipeline.Catchment) == "CatchmentC", :);
            testCase.verifyEqual(c.K_used, 3);
            testCase.verifyEqual(string(c.K_selection_method), "manual_override");
            testCase.verifyTrue(logical(c.K_override_applied));
        end

        function catchmentOutputsAreConsistent(testCase)
            catchments = ["CatchmentA", "CatchmentC"];
            required = ["filter_results_full.csv", "filter_results_coded.csv", ...
                "model_input_ages.csv", "excluded_ages.csv", ...
                "review_flags.csv", "output_summary.csv"];
            lookup = readtable(fullfile(testCase.OutputDir, "filter_code_lookup.csv"), ...
                Delimiter=",", VariableNamingRule="preserve", TextType="string");

            for catchment = catchments
                filter_dir = fullfile(testCase.OutputDir, catchment, "filter_output");
                for name = required
                    testCase.verifyTrue(isfile(fullfile(filter_dir, name)), ...
                        "Missing catchment output: " + name);
                end
                testCase.verifyEqual(numel(dir(fullfile(filter_dir, "*.csv"))), 6);

                full = readtable(fullfile(filter_dir, "filter_results_full.csv"), ...
                    Delimiter=",", VariableNamingRule="preserve", TextType="string");
                coded = readtable(fullfile(filter_dir, "filter_results_coded.csv"), ...
                    Delimiter=",", VariableNamingRule="preserve", TextType="string");
                model = readtable(fullfile(filter_dir, "model_input_ages.csv"), ...
                    Delimiter=",", VariableNamingRule="preserve", TextType="string");
                excluded = readtable(fullfile(filter_dir, "excluded_ages.csv"), ...
                    Delimiter=",", VariableNamingRule="preserve", TextType="string");
                review = readtable(fullfile(filter_dir, "review_flags.csv"), ...
                    Delimiter=",", VariableNamingRule="preserve", TextType="string");
                summary = readtable(fullfile(filter_dir, "output_summary.csv"), ...
                    Delimiter=",", VariableNamingRule="preserve", TextType="string");

                testCase.verifyEqual(height(coded), height(full));
                testCase.verifyEqual(string(coded.AnalysisID), string(full.AnalysisID));
                testCase.verifyEqual(height(model), nnz(logical(full.ModelInclude)));
                testCase.verifyTrue(all(logical(model.ModelInclude)));
                testCase.verifyEqual(height(excluded), nnz(string(full.Action) == "exclude"));
                testCase.verifyTrue(all(ismember(string(excluded.ReferenceCode), ["OR1", "OR2"])));
                testCase.verifyTrue(all(~logical(excluded.ModelInclude)));
                testCase.verifyEqual(height(review), nnz(logical(full.ReviewRecommended)));
                testCase.verifyTrue(all(logical(review.ReviewRecommended)));
                review_excluded = string(review.Action) == "exclude";
                testCase.verifyTrue(all(ismember( ...
                    string(review.ReferenceCode(review_excluded)), ["OR1", "OR2"])));
                testCase.verifyEqual(sum(summary.N), height(full));

                lookup_codes = string(lookup.Code);
                reference_codes = string(full.ReferenceCode);
                review_codes = string(full.ReviewCode);
                expected_reference_ids = nan(height(full), 1);
                expected_review_ids = nan(height(full), 1);
                for code_row = 1:height(lookup)
                    expected_reference_ids(reference_codes == lookup_codes(code_row)) = ...
                        lookup.CodeID(code_row);
                    expected_review_ids(review_codes == lookup_codes(code_row)) = ...
                        lookup.CodeID(code_row);
                end
                testCase.verifyFalse(any(isnan(expected_reference_ids)));
                testCase.verifyFalse(any(isnan(expected_review_ids)));
                testCase.verifyEqual(coded.ReferenceResultID, expected_reference_ids);
                testCase.verifyEqual(coded.ReviewFlagID, expected_review_ids);

                component_dir = fullfile(testCase.OutputDir, catchment, ...
                    "youngest_zircon_component");
                testCase.verifyTrue(isfile(fullfile(component_dir, ...
                    "youngest_zircon_component_plot.png")));
                testCase.verifyTrue(isfile(fullfile(component_dir, ...
                    "youngest_zircon_component_summary.csv")));
                testCase.verifyTrue(isfile(fullfile(testCase.OutputDir, catchment, ...
                    "sensitivity", "reference_boundary_comparison.csv")));
            end
        end

        function templatesUseOneSigma(testCase)
            template_dir = fullfile(testCase.RepoDir, "input_templates");
            files = dir(fullfile(template_dir, "*.csv"));
            for i = 1:numel(files)
                headers = string(readcell(fullfile(files(i).folder, files(i).name), ...
                    Range="1:1"));
                testCase.verifyTrue(any(contains(headers, "1sigerr")));
                testCase.verifyFalse(any(contains(headers, "2sigerr")));
            end
        end
    end
end
