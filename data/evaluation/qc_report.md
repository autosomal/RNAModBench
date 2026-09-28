# harmonisation QC report

Generated: 2026-09-20 02:34:29

## Registry

- (sample, tool) pairs: 540
- status: {'ok': 407, 'missing': 125, 'ok_zero': 5, 'empty': 3}

## Pending / gaps (2 rows)

| sample | tool | reason |
|---|---|---|
| Arabidopsis_fip37_rep3 | EpiNano_Error | discovered but empty/unreadable |
| Arabidopsis_WT_rep3 | EpiNano_Error | discovered but empty/unreadable |

## Callsets

- extracted callsets: 407 (4534460 rows)
- non-ok: {'missing': 125, 'ok_zero': 5, 'empty_parsed': 3}

- `empty_parsed` (3) = the RNA004 Dorado `m6A_guitar` model splits: those pileup files hold no call row after the modkit `name`-code split (the guitar model has no m6A code), so they parse to zero rows. Not a missing sample.

- additionally filled by R2Dtool liftover: 1 callsets (replicates the legacy pipeline never converted)

### Availability matrix (per sample == per replicate)

| tool | Arabidopsis_WT_rep1 | Arabidopsis_WT_rep2 | Arabidopsis_WT_rep3 | Arabidopsis_fip37_rep1 | Arabidopsis_fip37_rep2 | Arabidopsis_fip37_rep3 | Curlcake_IVT_rep1 | Curlcake_IVT_rep2_partial | Curlcake_IVT_rep3 | Curlcake_RNA004_IVT | Curlcake_m6A_rep1 | Curlcake_m6A_rep2 | E_IVT_neg1 | E_IVT_neg2 | E_ss_rd_RNA1 | E_ss_rd_RNA2 | HeLa_IVT_rep1 | HeLa_IVT_rep2 | HeLa_IVT_rep3 | HeLa_RNA004_IVT | HeLa_RNA004_WT | HeLa_WT1 | HeLa_WT2 | HeLa_WT3 | mESCs_Mettl3_KO | mESCs_Mettl3_WT | mES_KO | mES_WT |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| CHEUI_m5C | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 1 | 0 | 0 | 1 | 1 | 1 | 0 | 0 | 0 | 0 |
| CHEUI_m6A | 1 | 1 | 1 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 0 | 0 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| DENA | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 0 | 1 | 1 | 0 | 0 | 0 | 0 | 1 | 1 | 1 | 0 | 0 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| DRUMMER | 1 | 1 | 1 | 1 | 0 | 1 | 0 | 0 | 0 | 0 | 1 | 1 | 0 | 0 | 1 | 1 | 1 | 1 | 1 | 0 | 1 | 1 | 1 | 1 | 0 | 1 | 1 | 1 |
| Dorado_hac@v5.0.0_m6A@v1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Dorado_hac@v5.0.0_m6A_DRACH | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Dorado_hac@v5.0.0_m6A_DRACH@v1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Dorado_hac@v5.0.0_pseU@v1_otherMod | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Dorado_hac@v5.0.0_pseU_m6A | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Dorado_hac@v5.0.0_pseU_m6A_Psi | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Dorado_hac@v5.1.0_all | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Dorado_hac@v5.1.0_all_Psi | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Dorado_hac@v5.1.0_all_m5C | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Dorado_hac@v5.1.0_inosine_m6A | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Dorado_hac@v5.1.0_inosine_m6A_inosine | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Dorado_hac@v5.1.0_inosine_m6A_otherMod | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Dorado_hac@v5.1.0_m5C@v1_otherMod | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Dorado_hac@v5.1.0_m6A_DRACH@v1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Dorado_hac@v5.1.0_pseU@v1_otherMod | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Dorado_sup@v5.0.0_m6A@v1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Dorado_sup@v5.0.0_m6A_DRACH | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Dorado_sup@v5.0.0_m6A_DRACH@v1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Dorado_sup@v5.0.0_pseU@v1_otherMod | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Dorado_sup@v5.0.0_pseU_m6A | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Dorado_sup@v5.0.0_pseU_m6A_Psi | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Dorado_sup@v5.1.0_all | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Dorado_sup@v5.1.0_all_Psi | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Dorado_sup@v5.1.0_all_m5C | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Dorado_sup@v5.1.0_inosine_m6A | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Dorado_sup@v5.1.0_inosine_m6A_inosine | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Dorado_sup@v5.1.0_inosine_m6A_otherMod | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Dorado_sup@v5.1.0_m5C@v1_otherMod | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Dorado_sup@v5.1.0_m6A_DRACH@v1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Dorado_sup@v5.1.0_pseU@v1_otherMod | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ELIGOS2_diff | 1 | 1 | 1 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 | 1 | 0 | 0 | 1 | 1 | 1 | 1 | 1 | 0 | 1 | 1 | 1 | 1 | 0 | 1 | 1 | 1 |
| ELIGOS2_solo | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 0 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| EpiNano_Error | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 0 | 0 | 1 | 1 | 0 | 1 | 1 | 1 | 1 | 1 | 1 | 0 | 0 | 1 | 1 | 1 | 1 | 1 | 0 | 1 |
| MINES | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 0 | 1 | 1 | 0 | 0 | 0 | 0 | 1 | 1 | 1 | 0 | 0 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| NanoMUD_m1psi | 0 | 0 | 1 | 0 | 0 | 1 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 1 | 0 | 0 | 1 | 1 | 1 | 0 | 0 | 0 | 0 |
| NanoMUD_psi | 0 | 0 | 1 | 0 | 0 | 1 | 0 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 1 | 0 | 0 | 1 | 1 | 1 | 0 | 0 | 0 | 0 |
| NanoNm | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 1 | 0 | 0 | 1 | 1 | 1 | 0 | 0 | 0 | 0 |
| NanoPsu | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 0 | 0 | 0 | 0 |
| NanoSPA_m6A | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| NanoSPA_psU | 0 | 0 | 1 | 0 | 0 | 1 | 0 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 0 | 0 | 0 | 0 |
| Nanocompore | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 0 | 1 | 1 | 0 | 1 | 1 | 1 | 1 | 1 | 1 | 0 | 0 | 1 | 1 | 1 | 1 | 1 | 0 | 1 |
| Nanom6A | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 0 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 0 | 0 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| m6Anet | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| xPore | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 0 | 1 | 1 | 0 | 1 | 1 | 1 | 1 | 1 | 1 | 0 | 0 | 1 | 1 | 1 | 0 | 1 | 1 | 1 |
| yanocomp | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 0 | 0 | 1 | 1 | 0 | 1 | 1 | 1 | 1 | 1 | 1 | 0 | 0 | 1 | 1 | 1 | 0 | 1 | 1 | 1 |

## Liftover reproduction check

- cases compared against the legacy `*_liftover.txt` files: 93
- verdicts: {'identical': 88, 'different': 5}

## Universe

| sample | universe rows (cov>=5) | cov>=10 | DRACH |
|---|---|---|---|
| HeLa_WT1 | 3535909 | 2053505 | 150428 |
| HeLa_WT1 | 3173856 | 1877331 | 0 |
| HeLa_WT1 | 3468987 | 2014722 | 0 |
| HeLa_WT1 | 3468987 | 2014722 | 0 |
| HeLa_WT1 | 3535909 | 2053505 | 150428 |
| HeLa_WT1 | 13463884 | 7889565 | 150428 |
| HeLa_WT2 | 5268781 | 3478914 | 252602 |
| HeLa_WT2 | 4628949 | 3075227 | 0 |
| HeLa_WT2 | 5189761 | 3439400 | 0 |
| HeLa_WT2 | 5189761 | 3439400 | 0 |
| HeLa_WT2 | 5268781 | 3478914 | 252602 |
| HeLa_WT2 | 19876709 | 13179289 | 252602 |
| HeLa_WT3 | 2696539 | 1537920 | 114649 |
| HeLa_WT3 | 2563990 | 1514783 | 0 |
| HeLa_WT3 | 2660288 | 1514431 | 0 |
| HeLa_WT3 | 2660288 | 1514431 | 0 |
| HeLa_WT3 | 2696539 | 1537920 | 114649 |
| HeLa_WT3 | 10546402 | 6109730 | 114649 |
| HeLa_IVT_rep1 | 3634955 | 2314317 | 170399 |
| HeLa_IVT_rep1 | 3156486 | 2018934 | 0 |
| HeLa_IVT_rep1 | 3481656 | 2209896 | 0 |
| HeLa_IVT_rep1 | 3481656 | 2209896 | 0 |
| HeLa_IVT_rep1 | 3634955 | 2314317 | 170399 |
| HeLa_IVT_rep1 | 13633986 | 8705772 | 170399 |
| HeLa_IVT_rep2 | 4789231 | 3232443 | 239239 |
| HeLa_IVT_rep2 | 4292712 | 2903531 | 0 |
| HeLa_IVT_rep2 | 4558758 | 3057069 | 0 |
| HeLa_IVT_rep2 | 4558758 | 3057069 | 0 |
| HeLa_IVT_rep2 | 4789231 | 3232443 | 239239 |
| HeLa_IVT_rep2 | 18169213 | 12280552 | 239239 |
| HeLa_IVT_rep3 | 2775464 | 1695666 | 125307 |
| HeLa_IVT_rep3 | 2448393 | 1495111 | 0 |
| HeLa_IVT_rep3 | 2627434 | 1595917 | 0 |
| HeLa_IVT_rep3 | 2627434 | 1595917 | 0 |
| HeLa_IVT_rep3 | 2775464 | 1695666 | 125307 |
| HeLa_IVT_rep3 | 10467120 | 6396625 | 125307 |
| HeLa_RNA004_WT | 8680595 | 6832804 | 491709 |
| HeLa_RNA004_WT | 7658621 | 5987241 | 0 |
| HeLa_RNA004_WT | 8524924 | 6730709 | 0 |
| HeLa_RNA004_WT | 8524924 | 6730709 | 0 |
| HeLa_RNA004_WT | 8680595 | 6832804 | 491709 |
| HeLa_RNA004_WT | 32813998 | 25780925 | 491709 |
| HeLa_RNA004_IVT | 10207535 | 8040365 | 578170 |
| HeLa_RNA004_IVT | 8802293 | 6882313 | 0 |
| HeLa_RNA004_IVT | 9897322 | 7777287 | 0 |
| HeLa_RNA004_IVT | 9897322 | 7777287 | 0 |
| HeLa_RNA004_IVT | 10207535 | 8040365 | 578170 |
| HeLa_RNA004_IVT | 38160225 | 29973234 | 578170 |

## Annotation QC

- callsets annotated: 13
- callsets where <80 % of calls sit on the expected base (possible coordinate/strand issues or transcript-space output): 0

## Known caveats

- Mouse WT/KO come from two independent studies; they are evaluated separately and never merged (see metadata/replicate_structure.csv).
- Curlcake `Curlcake_IVT_rep2_partial` is a depth-matched subset of `Curlcake_IVT_rep3` (same run, SRR8767348) and is NOT an independent replicate.
- E. coli WT (`E_ss_rd_RNA1/2`) are two runs of the same condition; the IVT control `E_IVT_neg1/2` is ONE sample split into two halves (SRR27228854, 882k/506k local reads) and is only used as a null-vs-null false-positive control.
- Comparison outputs are filed under the tool's own *test side*; the two pairs whose BOTH sides are special samples carry the full comparison in the directory name instead: `E_IVT_neg1_vs_E_IVT_neg2` (the single null-vs-null comparison; owned by `E_IVT_neg1` for DRUMMER/ELIGOS2_diff and by `E_IVT_neg2` for EpiNano_DiffErr/xPore/yanocomp/Nanocompore) and `mESCs_Mettl3_KO_vs_mES_KO` (cross-study KO-vs-KO; owned by `mES_KO` for ELIGOS2_diff/DRUMMER/xPore/yanocomp and by `mESCs_Mettl3_KO` for EpiNano_DiffErr/Nanocompore). There is deliberately no comparison-index numbering.
- The Nanom6A fill samples were produced with the f5c-mode pipeline; their provenance is recorded in the callsets manifest.
- purified sites are WT-called/KO-absent on the common universe and carry the R3-7 circularity caveat.
