function results = run_release_tests()
%RUN_RELEASE_TESTS Run the public release test suite.
repo_dir = string(fileparts(mfilename("fullpath")));
tests_dir = fullfile(repo_dir, "tests");
addpath(repo_dir, tests_dir);
cleanup = onCleanup(@() rmpath(tests_dir)); %#ok<NASGU>

suite = testsuite(tests_dir, "IncludeSubfolders", true);
results = run(suite);
assertSuccess(results);
end
