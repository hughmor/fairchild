use crate::ParseError;

/// Parse a `.<directive> NAME [N]` vector-port declaration into
/// `(canonical name, channel count)`.  Shared by `.optical_port` and
/// `.electrical_port` so both reject a bad count the same way.
pub(super) fn parse_port_decl(
    line: &str,
    directive: &str,
    lineno: usize,
) -> Result<(String, usize), crate::ParseError> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    if tokens.len() < 2 {
        return Err(crate::ParseError::Syntax {
            line: lineno,
            msg: format!("{directive} needs a port name"),
        });
    }
    let channels = if tokens.len() >= 3 {
        tokens[2]
            .parse::<usize>()
            .ok()
            .filter(|&n| n > 0)
            .ok_or_else(|| crate::ParseError::Syntax {
                line: lineno,
                msg: format!("invalid channel count '{}' in {directive}", tokens[2]),
            })?
    } else {
        1
    };
    Ok((canon_node(tokens[1]), channels))
}

pub(super) fn canon_node(s: &str) -> String {
    let s = s.to_lowercase();
    if s == "gnd" {
        "0".to_string()
    } else {
        s
    }
}

/// Parse an SPICE suffix (k, meg, m, u, n, p, f, g, t) into a float.
pub(super) fn parse_value(s: &str, lineno: usize) -> Result<f64, ParseError> {
    let s_lc = s.to_lowercase();
    let (num_part, multiplier) = if let Some(n) = s_lc.strip_suffix("meg") {
        (n, 1e6)
    } else if let Some(n) = s_lc.strip_suffix('k') {
        (n, 1e3)
    } else if let Some(n) = s_lc.strip_suffix('m') {
        (n, 1e-3)
    } else if let Some(n) = s_lc.strip_suffix('u') {
        (n, 1e-6)
    } else if let Some(n) = s_lc.strip_suffix('n') {
        (n, 1e-9)
    } else if let Some(n) = s_lc.strip_suffix('p') {
        (n, 1e-12)
    } else if let Some(n) = s_lc.strip_suffix('f') {
        (n, 1e-15)
    } else if let Some(n) = s_lc.strip_suffix('g') {
        (n, 1e9)
    } else if let Some(n) = s_lc.strip_suffix('t') {
        (n, 1e12)
    } else {
        (s_lc.as_str(), 1.0)
    };
    num_part
        .parse::<f64>()
        .map(|v| v * multiplier)
        .map_err(|e| ParseError::BadNumber {
            value: s.to_string(),
            line: lineno,
            source: e,
        })
}

/// Expand a bus-vector token like `net[M..N]` into individual net names.
pub(super) fn expand_bus_vectors(token: &str) -> Vec<String> {
    if let (Some(lb), Some(rb)) = (token.find('['), token.rfind(']')) {
        if lb < rb {
            let base = &token[..lb];
            let range_str = &token[lb + 1..rb];
            if let Some((lo_s, hi_s)) = range_str.split_once("..") {
                if let (Ok(lo), Ok(hi)) =
                    (lo_s.trim().parse::<usize>(), hi_s.trim().parse::<usize>())
                {
                    return (lo..=hi).map(|i| format!("{}_{}", base, i)).collect();
                }
            }
        }
    }
    vec![token.to_string()]
}
