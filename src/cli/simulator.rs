//! CLI entry point for the `ciri-simulator` development binary.
//!
//! The simulator is intentionally separated from the user-facing `ciri`
//! pipeline. This module only parses command-line options and prints the stable
//! run summary; generation logic lives in `crate::simulator`.

use anyhow::Result;
use clap::Parser;

use crate::simulator::{run, SimulateArgs};

/// Parses simulator arguments, runs the generator, and prints its summary.
///
/// Keeping this wrapper thin ensures tests can call `simulator::run` directly
/// while the binary exposes the same argument contract to development scripts.
pub fn main() -> Result<()> {
    let summary = run(SimulateArgs::parse())?;
    eprintln!("{summary}");
    Ok(())
}
