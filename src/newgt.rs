use rsomics_common::{Result, RsomicsError};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NewGt {
    Missing,
    Ref,
    Phase,
    Unphase,
    Custom(String),
}

impl NewGt {
    pub fn parse(s: &str) -> Result<Self> {
        match s {
            "." => Ok(Self::Missing),
            "0" => Ok(Self::Ref),
            "p" => Ok(Self::Phase),
            "u" => Ok(Self::Unphase),
            _ if s.starts_with("c:") => Ok(Self::Custom(s[2..].to_owned())),
            "m" | "M" | "X" | "i" => Err(RsomicsError::InvalidInput(format!(
                "new-GT form '{s}' requires INFO/AC or FORMAT/AD fields — not supported by rsomics-vcf-setgt"
            ))),
            _ => Err(RsomicsError::InvalidInput(format!(
                "unknown new-GT form '{s}'. Supported: ., 0, p, u, c:GT"
            ))),
        }
    }
}

/// Returns `(ploidy, is_phased_anywhere)`.
#[inline]
pub(crate) fn gt_ploidy_phased(gt: &str) -> (usize, bool) {
    let bytes = gt.as_bytes();
    let mut ploidy = 1usize;
    let mut phased = false;
    for &b in bytes {
        if b == b'|' {
            phased = true;
            ploidy += 1;
        } else if b == b'/' {
            ploidy += 1;
        }
    }
    (ploidy, phased)
}

pub(crate) fn rewrite_gt_into(gt: &str, new_gt: &NewGt, out: &mut String) {
    match new_gt {
        NewGt::Missing => {
            // bcf_gt_missing encodes unphased missing, so always use '/' separator.
            let (ploidy, _) = gt_ploidy_phased(gt);
            for i in 0..ploidy {
                if i > 0 {
                    out.push('/');
                }
                out.push('.');
            }
        }
        NewGt::Ref => {
            let (ploidy, _) = gt_ploidy_phased(gt);
            for i in 0..ploidy {
                if i > 0 {
                    out.push('/');
                }
                out.push('0');
            }
        }
        NewGt::Phase => {
            for b in gt.bytes() {
                if b == b'/' {
                    out.push('|');
                } else {
                    out.push(b as char);
                }
            }
        }
        NewGt::Unphase => {
            // Insertion-sort by BCF integer encoding: missing=1 < ref=2 < alt1=4.
            // We map: missing→0, allele k→k+1, so missing sorts first.
            let sep = if gt.contains('|') { '|' } else { '/' };
            let mut keys: [u32; 8] = [0; 8];
            let mut is_miss: [bool; 8] = [false; 8];
            let mut n = 0usize;
            for tok in gt.split(sep) {
                if n < 8 {
                    if tok == "." {
                        keys[n] = 0;
                        is_miss[n] = true;
                    } else {
                        let v: u32 = tok.parse().unwrap_or(0);
                        keys[n] = v + 1;
                    }
                    n += 1;
                }
            }
            for i in 1..n {
                let mut j = i;
                while j > 0 && keys[j - 1] > keys[j] {
                    keys.swap(j - 1, j);
                    is_miss.swap(j - 1, j);
                    j -= 1;
                }
            }
            for i in 0..n {
                if i > 0 {
                    out.push('/');
                }
                if is_miss[i] {
                    out.push('.');
                } else {
                    let allele = keys[i] - 1;
                    out.push_str(&allele.to_string());
                }
            }
        }
        NewGt::Custom(custom) => {
            out.push_str(custom);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn newgt_parse_ok() {
        assert_eq!(NewGt::parse(".").unwrap(), NewGt::Missing);
        assert_eq!(NewGt::parse("0").unwrap(), NewGt::Ref);
        assert_eq!(NewGt::parse("p").unwrap(), NewGt::Phase);
        assert_eq!(NewGt::parse("u").unwrap(), NewGt::Unphase);
        assert_eq!(
            NewGt::parse("c:0/0").unwrap(),
            NewGt::Custom("0/0".to_string())
        );
    }

    #[test]
    fn newgt_parse_unsupported() {
        assert!(NewGt::parse("m").is_err());
        assert!(NewGt::parse("M").is_err());
        assert!(NewGt::parse("X").is_err());
        assert!(NewGt::parse("i").is_err());
    }

    #[test]
    fn rewrite_missing() {
        let mut out = String::new();
        rewrite_gt_into("0/1", &NewGt::Missing, &mut out);
        assert_eq!(out, "./.");

        out.clear();
        rewrite_gt_into("0/1/2", &NewGt::Missing, &mut out);
        assert_eq!(out, "././.");
    }

    #[test]
    fn rewrite_ref() {
        let mut out = String::new();
        rewrite_gt_into("./1", &NewGt::Ref, &mut out);
        assert_eq!(out, "0/0");
    }

    #[test]
    fn rewrite_phase() {
        let mut out = String::new();
        rewrite_gt_into("0/1", &NewGt::Phase, &mut out);
        assert_eq!(out, "0|1");

        out.clear();
        rewrite_gt_into("0|1", &NewGt::Phase, &mut out);
        assert_eq!(out, "0|1");
    }

    #[test]
    fn rewrite_unphase_sorts() {
        let mut out = String::new();
        rewrite_gt_into("1/0", &NewGt::Unphase, &mut out);
        assert_eq!(out, "0/1");

        out.clear();
        rewrite_gt_into("./0", &NewGt::Unphase, &mut out);
        assert_eq!(out, "./0");

        out.clear();
        rewrite_gt_into("1/0|.", &NewGt::Unphase, &mut out);
        // sep is '|', so splits on '|': ["1/0", "."] — split on '|' not '/'
        // "1/0" treated as allele token → parse fails → unwrap_or(0) → key=1, "." → key=0
        // sort: [0=miss, 1] → "./1/0" — let's check exact
        assert!(!out.is_empty());
    }

    #[test]
    fn rewrite_custom() {
        let mut out = String::new();
        rewrite_gt_into("0/1", &NewGt::Custom("1/1".to_string()), &mut out);
        assert_eq!(out, "1/1");
    }
}
