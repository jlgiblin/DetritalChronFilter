function manifest = convert_fixed_inputs_to_two_file(source_root, destination_root)
% CONVERT_FIXED_INPUTS_TO_TWO_FILE
% Converts the original four-file DetritalChronFilter input layout into the
% generic two-file layout used by v0.1.0-rc2.
%
% The source folder may be one sample containing the four files directly,
% or a parent folder with one subfolder per sample. Source files are read
% only. Existing converted CSVs are never overwritten.
%
% Original files:
%   ZrnPb.csv, HblAr.csv, ApHeApPb.csv, ZrnHeZrnPb.csv
%
% Converted files:
%   ReferenceDistribution.csv, ChronometerData.csv
%
% Conversion policy:
%   - zircon U-Pb distribution -> reference distribution
%   - hornblende Ar -> unpaired model candidate
%   - apatite He and U-Pb -> ordered pair; both are model candidates
%   - zircon He and U-Pb -> ordered pair; U-Pb is context only

arguments
    source_root      (1,1) string
    destination_root (1,1) string
end

if ~isfolder(source_root)
    error("DetritalChronFilter:MissingSourceFolder", ...
        "Source folder not found: %s", source_root);
end
if paths_match(source_root, destination_root)
    error("DetritalChronFilter:UnsafeConversionTarget", ...
        "Choose a destination folder different from the source folder.");
end

[sample_names, sample_dirs] = discover_fixed_layout(source_root);
converted = cell(numel(sample_names), 1);
manifest_rows = cell(numel(sample_names), 1);

for i = 1:numel(sample_names)
    sample_name = sample_names(i);
    sample_dir = sample_dirs(i);

    T_reference_old = read_required_table( ...
        fullfile(sample_dir, "ZrnPb.csv"), ...
        ["ZrnPbDate","ZrnPb1sigerr"]);
    T_hbl = read_required_table( ...
        fullfile(sample_dir, "HblAr.csv"), ...
        ["HblGrain","HblArDate","HblAr1sigerr"]);
    T_ap = read_required_table( ...
        fullfile(sample_dir, "ApHeApPb.csv"), ...
        ["ApGrain","ApHeDate","ApHe1sigerr","ApPbDate","ApPb1sigerr"]);
    T_zrn = read_required_table( ...
        fullfile(sample_dir, "ZrnHeZrnPb.csv"), ...
        ["ZrnGrain","ZrnHeDate","ZrnHe1sigerr","ZrnPbDate","ZrnPb1sigerr"]);

    reference_age = positive_numeric(T_reference_old.ZrnPbDate, ...
        "ZrnPbDate", sample_name);
    reference_uncertainty = positive_numeric(T_reference_old.ZrnPb1sigerr, ...
        "ZrnPb1sigerr", sample_name);
    reference_id = compose("Ref%03d", (1:height(T_reference_old))');
    T_reference = table(repmat("ZrnUPb", height(T_reference_old), 1), ...
        reference_id, reference_age, reference_uncertainty, ...
        'VariableNames', {'ReferenceSystem','GrainID','Age_Ma','Age_1sigma_Ma'});

    hbl_id = valid_ids(T_hbl.HblGrain, "HblGrain", sample_name);
    ap_id = valid_ids(T_ap.ApGrain, "ApGrain", sample_name);
    zrn_id = valid_ids(T_zrn.ZrnGrain, "ZrnGrain", sample_name);
    ap_pair_id = "Ap:" + ap_id;
    zrn_pair_id = "Zrn:" + zrn_id;

    T_chronometers = [ ...
        make_rows("HblAr", hbl_id, ...
            positive_numeric(T_hbl.HblArDate, "HblArDate", sample_name), ...
            positive_numeric(T_hbl.HblAr1sigerr, "HblAr1sigerr", sample_name), ...
            repmat("", height(T_hbl), 1), repmat("", height(T_hbl), 1), ...
            true(height(T_hbl), 1));
        make_rows("ApHe", ap_id, ...
            positive_numeric(T_ap.ApHeDate, "ApHeDate", sample_name), ...
            positive_numeric(T_ap.ApHe1sigerr, "ApHe1sigerr", sample_name), ...
            ap_pair_id, repmat("expected_younger", height(T_ap), 1), ...
            true(height(T_ap), 1));
        make_rows("ApUPb", ap_id, ...
            positive_numeric(T_ap.ApPbDate, "ApPbDate", sample_name), ...
            positive_numeric(T_ap.ApPb1sigerr, "ApPb1sigerr", sample_name), ...
            ap_pair_id, repmat("expected_older", height(T_ap), 1), ...
            true(height(T_ap), 1));
        make_rows("ZrnHe", zrn_id, ...
            positive_numeric(T_zrn.ZrnHeDate, "ZrnHeDate", sample_name), ...
            positive_numeric(T_zrn.ZrnHe1sigerr, "ZrnHe1sigerr", sample_name), ...
            zrn_pair_id, repmat("expected_younger", height(T_zrn), 1), ...
            true(height(T_zrn), 1));
        make_rows("ZrnUPb", zrn_id, ...
            positive_numeric(T_zrn.ZrnPbDate, "ZrnPbDate", sample_name), ...
            positive_numeric(T_zrn.ZrnPb1sigerr, "ZrnPb1sigerr", sample_name), ...
            zrn_pair_id, repmat("expected_older", height(T_zrn), 1), ...
            false(height(T_zrn), 1))];

    output_dir = fullfile(destination_root, sample_name);
    output_reference = fullfile(output_dir, "ReferenceDistribution.csv");
    output_chronometers = fullfile(output_dir, "ChronometerData.csv");
    if isfile(output_reference) || isfile(output_chronometers)
        error("DetritalChronFilter:ConvertedFileExists", ...
            "Converted CSV already exists for sample %s in %s.", ...
            sample_name, output_dir);
    end

    converted{i} = struct( ...
        "output_dir", output_dir, ...
        "reference", T_reference, ...
        "chronometers", T_chronometers);
    manifest_rows{i} = table(sample_name, height(T_reference), ...
        height(T_hbl), height(T_ap), height(T_zrn), height(T_chronometers), ...
        'VariableNames', {'Sample','N_ReferenceAges','N_HblAr', ...
        'N_ApatitePairs','N_ZirconPairs','N_ChronometerRows'});
end

for i = 1:numel(converted)
    item = converted{i};
    if ~isfolder(item.output_dir), mkdir(item.output_dir); end
    writetable(item.reference, ...
        fullfile(item.output_dir, "ReferenceDistribution.csv"));
    writetable(item.chronometers, ...
        fullfile(item.output_dir, "ChronometerData.csv"));
end

manifest = vertcat(manifest_rows{:});
fprintf("Converted %d sample(s) to: %s\n", height(manifest), destination_root);
disp(manifest);
end

function [sample_names, sample_dirs] = discover_fixed_layout(source_root)
required = ["ZrnPb.csv","HblAr.csv","ApHeApPb.csv","ZrnHeZrnPb.csv"];
root_has = arrayfun(@(name) isfile(fullfile(source_root, name)), required);
if any(root_has)
    if ~all(root_has)
        error("DetritalChronFilter:IncompleteFixedInput", ...
            "The source folder contains only part of the original four-file layout.");
    end
    normalized_root = char(source_root);
    while numel(normalized_root) > 1 && normalized_root(end) == filesep
        normalized_root(end) = [];
    end
    [~, sample_name] = fileparts(normalized_root);
    sample_names = string(sample_name);
    sample_dirs = source_root;
    return
end

entries = dir(source_root);
is_dot = strcmp({entries.name}, '.') | strcmp({entries.name}, '..');
subfolders = entries([entries.isdir] & ~is_dot);
names = string({subfolders.name});
dirs = fullfile(source_root, names);
valid = false(size(names));
for i = 1:numel(names)
    present = arrayfun(@(name) isfile(fullfile(dirs(i), name)), required);
    if any(present) && ~all(present)
        error("DetritalChronFilter:IncompleteFixedInput", ...
            "Sample %s contains only part of the original four-file layout.", ...
            names(i));
    end
    valid(i) = all(present);
end
sample_names = names(valid);
sample_dirs = dirs(valid);
if isempty(sample_names)
    error("DetritalChronFilter:NoFixedInputs", ...
        "No complete original four-file sample inputs found in %s.", source_root);
end
end

function T = read_required_table(filepath, required)
if ~isfile(filepath)
    error("DetritalChronFilter:MissingFixedInput", ...
        "Required source file not found: %s", filepath);
end
T = readtable(filepath, "VariableNamingRule","preserve", "TextType","string");
found = string(T.Properties.VariableNames);
missing = required(~ismember(required, found));
if ~isempty(missing)
    error("DetritalChronFilter:MissingFixedColumn", ...
        "File %s is missing required column(s): %s.", ...
        filepath, strjoin(missing, ", "));
end
end

function values = positive_numeric(values, column_name, sample_name)
if isnumeric(values)
    values = double(values(:));
else
    values = str2double(string(values(:)));
end
if any(~isfinite(values) | values <= 0)
    error("DetritalChronFilter:InvalidFixedValue", ...
        "Sample %s column %s contains a missing, nonnumeric, or nonpositive value.", ...
        sample_name, column_name);
end
end

function ids = valid_ids(values, column_name, sample_name)
ids = strtrim(string(values(:)));
ids(ismissing(ids)) = "";
if any(strlength(ids) == 0)
    error("DetritalChronFilter:InvalidFixedID", ...
        "Sample %s column %s contains a blank identifier.", ...
        sample_name, column_name);
end
if numel(unique(ids)) ~= numel(ids)
    error("DetritalChronFilter:DuplicateFixedID", ...
        "Sample %s column %s contains duplicate identifiers.", ...
        sample_name, column_name);
end
end

function T = make_rows(chronometer, grain_id, age, uncertainty, ...
        pair_id, pair_role, use_for_model)
n = numel(age);
T = table(repmat(string(chronometer), n, 1), string(grain_id), ...
    age(:), uncertainty(:), string(pair_id), string(pair_role), ...
    logical(use_for_model), ...
    'VariableNames', {'Chronometer','GrainID','Age_Ma','Age_1sigma_Ma', ...
    'PairID','PairRole','UseForModel'});
end

function tf = paths_match(first, second)
first = string(java.io.File(char(first)).getCanonicalPath());
second = string(java.io.File(char(second)).getCanonicalPath());
tf = first == second;
end
