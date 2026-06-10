# Changelog

All notable changes to the published `fountain_engine` crate are documented here.

## [1.3.1] - 2026-06-08

### Added

- **Native padding API** for schemes where the application source block size K is smaller than the internal block size K′ (`CodeParams.k`). Implicit zero padding symbols are not transmitted; the engine installs them during encode/decode setup.

  Consumers can pass the application source count directly:

  ```rust
  Encoder::new_with_num_source(&scheme, k_app);
  Decoder::new_with_num_source(&scheme, k_app);
  ```

  Existing constructors (`Encoder::new`, `Decoder::new`, and operator variants) default to `num_source == CodeParams.k` (no padding).

- **`DataManager`**: `set_num_source`, `num_source`, `num_padding`, `has_padding`, and `solver_type` for session-level padding configuration.
- **`Solver::install_padding`**: registers implicit padding for systematic and ordinary decode paths after solver construction.
- **Ordinary encode/decode**: updated `prepare_for_ordinary` and `restore_for_ordinary` to account for padding in the active message range.
