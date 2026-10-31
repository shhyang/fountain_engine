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

## What's new in 2.1

- `Operation::{EnsureZeroRange, CopyToN, MoveToN}` plus matching `DataManager` / `DataOperator` helpers
- LDPC `for_each_*` streaming enumerators (Vec adjacency APIs still work)
- HDPC `mul_binary_streamed` / `mul_binary_words` / `mul_binary_supports`; `mul_data` takes contiguous ID ranges
- MIT license (from 2.0.2)

## Solver

Encoding and decoding use `SystemSolver`: belief propagation, inactivation, and Gaussian elimination over a sparse master-system representation. This is the only solver in the published 2.x crate.

Profiling / timing introspection and staging solver experiments are monorepo development facilities and are **not** part of this published package.

## Dependency

```toml
[dependencies]
fountain_engine = "2.1.0"
```

## License

MIT License. See [LICENSE](LICENSE).

AI / machine-learning training on this code requires a separate license
from the copyright holder; see [NOTICE](NOTICE). That notice does not
change the MIT grant for ordinary use, modification, and redistribution.

Copyright (c) 2025, 2026 Shenghao Yang.
