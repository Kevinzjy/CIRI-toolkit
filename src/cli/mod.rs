//! Command-line front ends for the two public binaries.
//!
//! `ciri` is the user-facing analysis entry point, while `ciri-simulator` is a
//! development fixture generator.  Keeping the CLI glue here prevents the thin
//! `src/bin/*` wrappers from accumulating pipeline or simulator logic.

pub mod ciri;
pub mod simulator;
