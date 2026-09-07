# Roadmap: modular component selection and chronometer input design

Status: future design only. The features below are **not implemented in v0.1.0**. The supported release uses the documented four-file input layout, GMM component estimation, BIC or a manual K override, and youngest-eligible component selection.

## 1. Separate estimation from selection

The public interface should distinguish:

1. `ComponentMethod`: how candidate age components are estimated.
2. `ComponentSelection`: which candidate is used as the target.
3. `WindowMethod`: how the selected component becomes an age window.

This prevents “youngest component” from concealing multiple scientific choices.

### Proposed component-estimation methods

| Option | Intended use | Important limitation |
|---|---|---|
| `gmm` | Current parametric Gaussian-mixture workflow | Component number and Gaussian form require checking |
| `kde_peak` | Exploratory mode finding without assuming Gaussian components | Peaks depend on bandwidth and are not automatically geological populations |
| `user_window` | Reproduce a published or independently justified target window | Not a statistical estimate from the supplied sample |
| `external_membership` | Import component labels or probabilities from another reviewed workflow | Interpretation depends on the external method |

For GMM, `ModelCriterion` should allow `bic` (default) or `aic`, while preserving an explicit component-count override.

### Proposed component-selection rules

| Option | Selection rule |
|---|---|
| `youngest_eligible` | Youngest component whose mean lies in `TargetComponentAgeRange` (current behavior) |
| `nearest_age` | Component mean nearest a user-supplied `TargetAge_Ma` |
| `largest_weight` | Highest-weight component within `TargetComponentAgeRange` |
| `component_index` | User selects a component after inspecting QA plots |
| `user_window` | No component selection; supplied bounds are used directly |

Every output should record the method, selection rule, parameter values, selected component, and whether the choice was automatic or user-specified.

Youngest-grain and youngest-cluster estimators used for maximum depositional age should not be presented as interchangeable target-component estimators. They answer a different question unless the study objective is specifically maximum depositional age.

## 2. Generic chronometer schema

Use one long-format observation table rather than a separate hard-coded file for every mineral–method combination. Keep the current four-file layout as an import adapter for existing projects.

Required observation fields:

- `GrainID`
- `SystemID`
- `Age_Ma`
- `Age_1SE_Ma`

Optional paired-age fields:

- `PairedAge_Ma`
- `PairedAge_1SE_Ma`
- `PairType`

An unpaired apatite or zircon (U–Th)/He result leaves the paired fields blank. It can still be evaluated against the reference age, but it cannot be assigned a crystallization-to-cooling interval class.

A `SystemID` is only a display/grouping label; it does not imply a closure temperature or change the reference-age calculation. The presence of a valid paired crystallization age determines whether the paired-age checks apply. This supports hornblende Ar/Ar, paired or unpaired apatite He, paired or unpaired zircon He, and user-defined systems without changing the core classifier.

## 3. Scope boundary with MultichronFitTSF

DetritalChronFilter does not need closure-temperature metadata. Its common operation is age-based: all cooling-age observations are assessed relative to the selected zircon U–Pb reference, and observations with paired crystallization ages receive the additional age-order and interval checks.

Chronometer temperature ordering belongs in MultichronFitTSF, where it is used to generate a synthetic age–elevation transect. Keeping that information downstream avoids implying that DetritalChronFilter applies different reference-age rules to different chronometers.

If closure-temperature ranges or kinetic models are later made user-adjustable, they should be implemented and documented in MultichronFitTSF rather than added to this filter.

## 4. Implementation order

1. Add and validate the generic long-format importer while preserving current numerical results.
2. Add paired/unpaired handling; determine interval-screen applicability directly from the presence of paired ages.
3. Refactor GMM estimation and component selection into separate functions.
4. Add `user_window`, GMM criterion choice, and component-selection rules.
5. Add KDE peak selection only with bandwidth sensitivity and bootstrap stability outputs.
6. Hand the screened ages and system labels to MultichronFitTSF; keep temperature-ordering logic there.

This order keeps the public code flexible without implying that all user-selectable settings are scientifically equivalent.

## Supporting literature

- Vermeesch, P. (2012), *On the visualisation of detrital age distributions*, Chemical Geology, 312–313, 190–194. https://doi.org/10.1016/j.chemgeo.2012.04.021
- Coutts, D.S., Matthews, W.A., and Hubbard, S.M. (2019), *Assessment of widely used methods to derive depositional ages from detrital zircon populations*, Geoscience Frontiers, 10, 1421–1435. https://doi.org/10.1016/j.gsf.2018.11.002
