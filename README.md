# rsomics-vcf-setgt

Rewrite genotypes (the `GT` FORMAT field) in a VCF — a Rust port of
`bcftools +setGT`. Set genotypes to missing/ref/custom, (un)phase them, or
target specific samples by an expression. Record-exact with `bcftools +setGT`
and ~2.5× faster.

## Install

```sh
cargo install rsomics-vcf-setgt
```

## Usage

```sh
rsomics-vcf-setgt --target . -n 0 input.vcf                    # missing GTs → ref
rsomics-vcf-setgt -t q -i 'FMT/DP<5' -n . input.vcf -o out.vcf # DP<5 samples → missing
rsomics-vcf-setgt --target a -n c:0/0 input.vcf -o out.vcf     # all GTs → 0/0
```

| flag | meaning |
|---|---|
| `--target <T>` | `.` any-missing, `./x` partial-missing, `./.` full-missing, `a` all, `q` expression |
| `-n, --new-gt <V>` | `.` missing, `0` ref, `p` phase, `u` unphase, `c:GT` custom |
| `-i, --include <EXPR>` / `-e, --exclude <EXPR>` | for `-t q`: sample-level expression (via `rsomics-vcf-expr`) |
| `-o, --output <FILE>` | output VCF (default stdout) |

## Origin

Independent Rust reimplementation of `bcftools +setGT` based on the public plugin
documentation, the VCF spec, and black-box testing against `bcftools +setGT`
(1.23.1). Expression targeting uses the [`rsomics-vcf-expr`](https://crates.io/crates/rsomics-vcf-expr)
Layer-A engine. No GPL/MIT upstream source was used as reference.

License: MIT OR Apache-2.0.
Upstream credit: [bcftools](https://www.htslib.org/) (MIT/Expat).
