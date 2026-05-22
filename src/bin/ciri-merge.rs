//! Thin binary wrapper for the cohort-level `ciri-merge` command.

use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

fn main() -> anyhow::Result<()> {
    ciri_toolkit::cli::merge::main()
}
