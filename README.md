# fountain_engine

Core algorithms for fountain code encoding and decoding.

## What It Provides

- **Encoder / Decoder**: Generic fountain code encoder and decoder that work with any code scheme implementing the `CodeScheme` trait.
- **Traits**: `CodeScheme` and `DataOperator` — implement these to define custom fountain codes and data backends.
- **BinaryMatrix**: Packed GF(2) matrix utilities used by the solver stack.

## Architecture

The engine is data-free by design. Instead of operating on data directly, the encoder and decoder emit a sequence of `Operation` values that a `DataOperator` executes on actual data vectors. This separation enables:

- In-memory operation (`VecDataOperater` in `fountain_utility`)
- I/O logging (`IoDataOperator` in `fountain_utility`)
- Network streaming or any custom backend

## Solver

Encoding and decoding use `SystemSolver`: belief propagation, inactivation, and Gaussian elimination over a sparse master-system representation. This is the only solver in the published 2.x crate.

Profiling / timing introspection and staging solver experiments are monorepo development facilities and are **not** part of this published package.

## Dependency

```toml
[dependencies]
fountain_engine = "2.0.0"
```

## License

This crate is dual-licensed:

1. **AGPL-3.0** — Free for open-source use. See [LICENSE-AGPL](LICENSE-AGPL).
2. **Commercial** — For proprietary use without AGPL obligations. See [LICENSE-COMMERCIAL](LICENSE-COMMERCIAL).

Copyright (c) 2025 Shenghao Yang. All rights reserved.
