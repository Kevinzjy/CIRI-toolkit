//! Thin binary wrapper for the development-only `ciri-simulator` command.

use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

fn main() -> anyhow::Result<()> {
    ciri_toolkit::cli::simulator::main()
}
