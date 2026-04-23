# ADR 0005 — Licensing audit: root LICENSE vs upstream GPL v2

Status: accepted (2026-04-23, during rescue Phase 2); remediation committed in the same task as this ADR

## Context

Phase 0 baseline verification flagged an inconsistency:

- Repository-root `LICENSE` was the **GNU Lesser General Public License v2.1**. 179 lines; first line referenced `lgpl-2.1.en.html`.
- Upstream SWAP 4.2.0 (preserved on `legacy/swap-4.2.0` at `source_swp_4.2.0/LICENSE`) declares the **GNU General Public License v2**. TTUTIL 4.27 is separately LGPL v2.1.

The modernization tree in `src/` is a derivative work of the upstream GPL v2 SWAP code. Under GPL v2 §2(b) a derivative work must be licensed "as a whole at no charge to all third parties under the terms of this License" — i.e., GPL v2 (or v2+). It cannot be redistributed under LGPL v2.1 without permission from every copyright holder of the original SWAP work. For SWAP that means Wageningen University & Research (WUR), not the current modernization maintainer.

The mismatch is almost certainly accidental inheritance — TTUTIL's LGPL text was copied into the root `LICENSE` slot at some point, or the LGPL v2.1 text was dropped in as a placeholder. There is no evidence of an explicit relicensing agreement with WUR.

## Decision

The repository's license is corrected to GPL v2 to match upstream SWAP 4.2.0's license. The TTUTIL subset, when included verbatim or in derivative form, retains its LGPL v2.1 notice in its own file (preserved at build time via `subprojects/ttutil/` and at rest in `legacy/swap-4.2.0`).

Remediation committed in the same task as this ADR:

1. Root `LICENSE` replaced with the full GPL v2 text from https://www.gnu.org/licenses/old-licenses/gpl-2.0.txt.
2. `NOTICE.md` added clarifying that SWAP portions are GPL v2 and bundled TTUTIL portions are LGPL v2.1.

## Consequences

- **Positive**: the modernization repo complies with GPL v2 §2(b) going forward. No legal ambiguity for downstream users.
- **Positive**: TTUTIL's LGPL v2.1 is preserved at the per-file level, where it was always intended to live.
- **Negative**: anyone who previously relied on the root-level LGPL v2.1 notice to justify linking SWAP into proprietary code was mistaken. This ADR + fix makes the actual licensing explicit and forecloses that interpretation.
- **Negative**: the correction is a historical fix; downstream consumers who pulled the repo during the LGPL-v2.1-at-root period received it under that notice. We cannot retroactively change their copy's terms. Going forward, every new clone sees GPL v2.

## When to revisit

Closed by the remediation commit. If upstream SWAP ever relicenses (later GPL version, dual-license), that triggers ADR 0006.
