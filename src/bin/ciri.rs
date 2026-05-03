//! Thin binary wrapper for the user-facing `ciri` command.

use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

fn main() -> anyhow::Result<()> {
    ciri_toolkit::cli::ciri::main()
}
