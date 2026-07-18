# Changelog

All notable changes to the published `fountain_engine` crate are documented here.

## [2.0.1] - 2026-07-18

### Fixed

- **Systematic `decode_status` with padding:** compare received source count to
  `DataManager::num_source()` (application \(K\)) instead of `CodeParams::num_message()`
  (block \(K'\)). With `num_source < K'`, `decode_status()` previously never reported
  `Decoded` from the systematic early path even when every source symbol had arrived.
- **Post-solve message recovery:** iterate `received_msg_vectors` by its natural length
  (\(K\)) instead of `.take(num_message())` (\(K'\)).

## [2.0.0] - 2026-07-17

### Changed

- **Solver cutover (breaking):** `Encoder` / `Decoder` are built on `core::system_solver::SystemSolver` (sparse master system + inactivation + Gaussian elimination) instead of the legacy belief-propagation / `Solver` stack from 1.x.
- **Public API:** the crate root re-exports `BinaryMatrix` alongside the existing encoder, decoder, traits, and types.

### Removed

- **Legacy v1 stack:** `encoder_v1`, `decoder_v1`, and the old BP/`Solver` core modules are not part of the package.
- **Feature flags:** `legacy_solver`, `next_solver`, and `profiling` are not published. Profiling and staging solvers remain monorepo development facilities only.

### Migration

- Applications already using the default (non-`legacy_solver`) 1.x `Encoder` / `Decoder` API typically need only a version bump to `2.0`.
- Applications that depended on `legacy_solver` or on published `profiling` APIs must migrate to the `SystemSolver`-based API; those features are not available in the published 2.x crate.

## [1.3.2] - 2026-06-19

### Changed

- **Ordinary precoding:** refactored `core/precode.rs` — HDPC factorization of \(I' + D_s S_h\) moved to [`HDPC::lu_idssh`](src/traits/hdpc.rs); entry point renamed to `ordinary_precode_encode` (called from `Encoder::precode_encode`).
- **LDPC encode:** default implementation uses `broadcast_add_owned` to avoid extra slice clones.
- **Encoder:** systematic padding slots use `ensure_zero_one`; LT encoding uses `add_to_vector_owned`.

### Added

- **`HDPC::lu_idssh`:** default LU of \(I' + D_s S_h\) for ordinary precoding (override for custom factorizations).
- **`DataManager::broadcast_add_owned`:** take ownership of target ID lists when recording broadcast ops.
- **`Operation::EnsureZeroOne`** / **`DataManager::ensure_zero_one`:** zero a single vector without allocating a `Vec`.

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
