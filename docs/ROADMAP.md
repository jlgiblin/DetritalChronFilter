# Roadmap

DetritalChronFilter now uses a generic two-file input design:

- `ReferenceDistribution.csv` supplies the complete age distribution used to estimate the target component.
- `ChronometerData.csv` supplies one or more user-labelled chronometer systems in long format.
- Optional `PairID`, `PairRole`, and `UseForModel` fields support paired-age review and context-only observations without interpreting chronometer names.

The items below are possible future additions, not features of the current release candidate.

## 1. Additional component-estimation methods

The current implementation fits Gaussian mixture models and chooses the youngest eligible component, using either BIC or a reviewed manual K override. Future versions could separate three choices more explicitly:

1. `ComponentMethod`: how candidate components are estimated.
2. `ComponentSelection`: which candidate becomes the target.
3. `WindowMethod`: how the selected component becomes an age window.

Possible methods include KDE peak finding, a directly supplied user window, or externally supplied component memberships. Each would require its own validation and transparent provenance in the output; these methods should not be presented as scientifically interchangeable.

## 2. Additional target-selection rules

Possible future selection rules include the component nearest a user-supplied age, the largest-weight eligible component, or a component chosen after inspecting the QA plot. The current `TargetComponentAgeRange` and manual K override already allow users to constrain and review the GMM result without silently replacing it.

## 3. More complex analytical relationships

The current pair design intentionally supports either:

- two observations with blank roles, which records a relationship without testing age order; or
- exactly one `expected_younger` and one `expected_older` observation.

A future version could support groups containing more than two analyses or other user-defined relationship tests. Such an extension should retain explicit roles and avoid inferring geologic meaning from chronometer labels.

## 4. Scope boundary with MultichronFitTSF

DetritalChronFilter does not use closure-temperature metadata or impose an ordering among unpaired chronometer distributions. Its reference screen is the same for every row marked `UseForModel=true`; paired-age order is evaluated only when the user explicitly supplies pair roles.

Chronometer temperature ordering belongs in MultichronFitTSF, where it can inform construction and review of synthetic age-elevation transects. Keeping that logic downstream avoids implying that this filter applies different reference-age rules to different chronometer systems.

## Supporting literature

- Vermeesch, P. (2012), *On the visualisation of detrital age distributions*, Chemical Geology, 312–313, 190–194. https://doi.org/10.1016/j.chemgeo.2012.04.021
- Coutts, D.S., Matthews, W.A., and Hubbard, S.M. (2019), *Assessment of widely used methods to derive depositional ages from detrital zircon populations*, Geoscience Frontiers, 10, 1421–1435. https://doi.org/10.1016/j.gsf.2018.11.002
