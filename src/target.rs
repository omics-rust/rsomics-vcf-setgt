use rsomics_common::{Result, RsomicsError};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Target {
    AnyMissing,
    PartialMissing,
    FullyMissing,
    All,
    Query,
}

impl Target {
    pub fn parse(s: &str) -> Result<Self> {
        match s {
            "." => Ok(Self::AnyMissing),
            "./x" => Ok(Self::PartialMissing),
            "./." => Ok(Self::FullyMissing),
            "a" => Ok(Self::All),
            "q" => Ok(Self::Query),
            "b" | "nb" | "np" | "miss" => Err(RsomicsError::InvalidInput(format!(
                "target '{s}' requires bcftools binomial/missing engine — not supported by rsomics-vcf-setgt"
            ))),
            _ => Err(RsomicsError::InvalidInput(format!(
                "unknown target selector '{s}'. Supported: ., ./x, ./., a, q"
            ))),
        }
    }
}

/// Returns `(n_missing, n_alleles)` for a GT field string.
#[inline]
pub(crate) fn count_missing(gt: &str) -> (usize, usize) {
    let bytes = gt.as_bytes();
    let mut n_allele = 1usize;
    let mut n_miss = 0usize;

    if bytes.first() == Some(&b'.') {
        n_miss += 1;
    }

    for &b in bytes {
        if b == b'/' || b == b'|' {
            n_allele += 1;
        }
    }

    let mut after_sep = false;
    for &b in bytes.iter().skip(1) {
        if after_sep {
            if b == b'.' {
                n_miss += 1;
            }
            after_sep = false;
        }
        if b == b'/' || b == b'|' {
            after_sep = true;
        }
    }

    (n_miss, n_allele)
}

/// `Target::Query` always returns false; the caller uses `EvalContext` per sample.
#[inline]
pub(crate) fn matches_target(gt: &str, target: Target) -> bool {
    match target {
        Target::All => true,
        Target::Query => false,
        Target::AnyMissing | Target::PartialMissing | Target::FullyMissing => {
            let (n_miss, n_allele) = count_missing(gt);
            match target {
                Target::AnyMissing => n_miss > 0,
                Target::PartialMissing => n_miss > 0 && n_miss < n_allele,
                Target::FullyMissing => n_miss == n_allele,
                Target::All | Target::Query => unreachable!(),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn target_parse_ok() {
        assert_eq!(Target::parse(".").unwrap(), Target::AnyMissing);
        assert_eq!(Target::parse("./x").unwrap(), Target::PartialMissing);
        assert_eq!(Target::parse("./.").unwrap(), Target::FullyMissing);
        assert_eq!(Target::parse("a").unwrap(), Target::All);
        assert_eq!(Target::parse("q").unwrap(), Target::Query);
    }

    #[test]
    fn target_parse_unsupported() {
        assert!(Target::parse("b").is_err());
        assert!(Target::parse("nb").is_err());
        assert!(Target::parse("miss").is_err());
    }

    #[test]
    fn target_parse_unknown() {
        assert!(Target::parse("z").is_err());
    }

    #[test]
    fn count_missing_basic() {
        assert_eq!(count_missing("0/1"), (0, 2));
        assert_eq!(count_missing("./."), (2, 2));
        assert_eq!(count_missing("./1"), (1, 2));
        assert_eq!(count_missing("."), (1, 1));
        assert_eq!(count_missing("0/./1"), (1, 3));
    }

    #[test]
    fn matches_target_any() {
        assert!(matches_target("./.", Target::AnyMissing));
        assert!(matches_target("./1", Target::AnyMissing));
        assert!(!matches_target("0/1", Target::AnyMissing));
    }

    #[test]
    fn matches_target_partial() {
        assert!(!matches_target("./.", Target::PartialMissing));
        assert!(matches_target("./1", Target::PartialMissing));
    }

    #[test]
    fn matches_target_full() {
        assert!(matches_target("./.", Target::FullyMissing));
        assert!(!matches_target("./1", Target::FullyMissing));
    }

    #[test]
    fn matches_target_all() {
        assert!(matches_target("0/1", Target::All));
        assert!(matches_target("./.", Target::All));
    }

    #[test]
    fn matches_target_query_always_false() {
        assert!(!matches_target("0/1", Target::Query));
    }
}
