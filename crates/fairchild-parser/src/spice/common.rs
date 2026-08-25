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

/// SPICE's engineering multipliers, longest spelling first so `meg` is matched
/// before `m` and before `g`.
///
/// One table, read by both the trailing form (`4.7k`) and the RKM infix form
/// (`4k7`). Two tables would be two chances for the same letter to mean
/// different things depending on where in the token it appeared.
const MULTIPLIERS: &[(&str, f64)] = &[
    ("meg", 1e6),
    ("k", 1e3),
    ("g", 1e9),
    ("t", 1e12),
    ("m", 1e-3),
    ("u", 1e-6),
    ("n", 1e-9),
    ("p", 1e-12),
    ("f", 1e-15),
];

/// A SPICE value, with an engineering suffix or in RKM form.
///
/// ## RKM (`4k7`, `2n2`, `1meg5`)
///
/// The multiplier letter stands in for the decimal point, which is the notation
/// silkscreens and BOMs use because a period does not survive photocopying.
/// ngspice does **not** support it: it reads `4k7` as 4000 and drops the `7`
/// without a word. That is why this is safe to add rather than a compatibility
/// hazard — a deck validated against ngspice cannot contain RKM and mean
/// anything by it, so the only decks that reach this path are ones written
/// deliberately for fairchild.
///
/// Recognised only when digits sit on *both* sides of a known multiplier, so
/// `10nF` and `1e3` are untouched and every token that parses today parses to
/// the same number.
///
/// ## `m` is refused in RKM position, on purpose
///
/// SPICE reads `m` as milli and `meg` as mega. RKM, coming from component
/// marking, reads `M` as mega. So `4M7` means 4.7 MΩ to the person who wrote it
/// and 4.7 mΩ to a SPICE parser — a factor of 10⁹, in a token whose whole point
/// is being read at a glance. Neither reading is safe to pick, so both are
/// refused with an error naming the two unambiguous spellings.
pub(crate) fn parse_value(s: &str, lineno: usize) -> Result<f64, ParseError> {
    let s_lc = s.to_lowercase();

    // RKM first: it is the only form with a letter in the middle, so it cannot
    // shadow the trailing-suffix form below.
    if let Some(alpha) = s_lc.find(|c: char| c.is_ascii_alphabetic()) {
        let (head, rest) = s_lc.split_at(alpha);
        let rkm = MULTIPLIERS
            .iter()
            .find(|(name, _)| rest.starts_with(name))
            .filter(|(name, _)| {
                let tail = &rest[name.len()..];
                !head.is_empty()
                    && !tail.is_empty()
                    && tail.bytes().all(|b| b.is_ascii_digit())
                    && head
                        .bytes()
                        .all(|b| b.is_ascii_digit() || b == b'-' || b == b'+')
            });
        if let Some((name, mult)) = rkm {
            // The letter put back as the decimal point it stands in for.
            let joined = format!("{head}.{}", &rest[name.len()..]);
            if *name == "m" {
                return Err(ParseError::Syntax {
                    line: lineno,
                    msg: format!(
                        "'{s}' is ambiguous: SPICE reads 'm' as milli, but the RKM \
                         notation this looks like reads 'M' as mega — a factor of 1e9 \
                         apart. Write '{joined}m' for milli or '{joined}meg' for mega"
                    ),
                });
            }
            return joined
                .parse::<f64>()
                .map(|v| v * mult)
                .map_err(|e| ParseError::BadNumber {
                    value: s.to_string(),
                    line: lineno,
                    source: e,
                });
        }
    }

    let (num_part, multiplier) = MULTIPLIERS
        .iter()
        .find_map(|(name, mult)| s_lc.strip_suffix(name).map(|n| (n, *mult)))
        .unwrap_or((s_lc.as_str(), 1.0));
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
