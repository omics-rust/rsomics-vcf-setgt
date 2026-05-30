//! Genotype rewriter matching `bcftools +setGT` semantics.
//!
//! Supported `-t` (target) selectors:
//!   - `.`   any missing (partial or fully missing)
//!   - `./x` partially missing (≥1 allele missing, not all)
//!   - `./.` fully missing (all alleles missing)
//!   - `a`   all genotypes
//!   - `q`   filter-expression query — requires `-i` or `-e` via `SetGtConfig`
//!
//! Supported `-n` (new-GT) forms:
//!   - `.`   set all alleles to missing
//!   - `0`   set all alleles to REF (unphased)
//!   - `p`   phase existing genotype
//!   - `u`   unphase and sort alleles ascending
//!   - `c:GT` custom literal genotype (e.g. `0/0`, `0|1`)
//!
//! The `b`, `nb`, `np`, `miss` target forms and `m`/`M`/`X`/`i` new-GT forms
//! require FORMAT tag lookups beyond what this crate implements — they exit
//! non-zero with a clear message.
//!
//! ## Origin
//!
//! This crate is an independent Rust reimplementation of `bcftools +setGT`
//! (plugin in the bcftools suite) based on:
//! - Reading the MIT-licensed upstream source (setGT.c, bcftools develop branch)
//! - The public VCF 4.3 format specification
//! - Black-box behaviour testing against `bcftools +setGT`
//!
//! Source consulted: <https://github.com/samtools/bcftools/blob/develop/plugins/setGT.c>
//! License: MIT OR Apache-2.0
//! Upstream credit: bcftools <https://github.com/samtools/bcftools> (MIT)

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::too_many_lines
)]

pub mod driver;
pub mod newgt;
pub mod record;
pub mod target;

pub use driver::{SetGtConfig, SetGtStats, set_gt};
pub use newgt::NewGt;
pub use record::rewrite_record_into;
pub use target::Target;
