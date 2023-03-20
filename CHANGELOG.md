# Change Log

## Unreleased
### Changed
* fixed readme
* fixed missing text with unimplemented `_offer_skip_option`

## v0.4.1
### Changed
* Readme reflects new behavior
* Readme improved
* improved setup.py
* fixed a helper function

## v0.4.1
### Changed
* `./example` dir is actually included, now

## v0.4.0
### Added
* Made config.py/template.sbt accessible
  * Can now create (semi-)generic config & template files via `ncbi_submit example`
* `./example` dir should now be part of package
### Changed
* Now check that outdir is not a file and create before doing much else
### Removed
* `./config` directory - merged contents into `./example`

## v0.3.2 (First fully-functional version)
### Changed
* Fixed some issues with logging
* Edited to README to fit with recent changes
* Fixed some field name access issues (how config variable `sars_cov_2_diag_pcr_ct_values` is used)

## v0.3.1
### Changed
* Reconfigured a few arguments in `arguments.py` for better command line usage

## v0.3.0
### Added
* Can get version from within scripts
* Logging
### Changed
* How version number was accessed

## v0.2.2 - v0.2.8
### Changed
* Dependency edits because PyPI and test.PyPI had me confused

## v0.2.0 / v0.2.1
### Changed
* Restructured code to fit poetry's typical packaging format

## v0.1.0
### Added
* First attempt at reconfiguring and packaging this project
  * All features exist within the main code
  * Still plenty of bugs
* Versioning
