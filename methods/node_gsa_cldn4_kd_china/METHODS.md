# Methods

Searches were run on 2026-09-21 against the public HTTP APIs. No login, no controlled-data request, no FASTQ or matrix download.

## NODE

`POST https://www.biosino.org/node/api/app/browse/search` with `queryWord` set to `CLDN4`, `claudin-4`, `claudin`, `CLDN`, `knockdown`, `shRNA`, `siRNA`, `PRJCA023797`, `HRA006761`, `CRAD`, and `ELFN1`. Pages of 30 were followed until `totalElements` was covered. A record was counted as a CLDN4 mention when the returned JSON contained `CLDN4`, `claudin-4`, `claudin 4`, or `claudin4`.

## CNCB

`GET https://ngdc.cncb.ac.cn/search/api/specific` with `db` in `gsa`, `hra` (GSA-Human), `omix`, and `bioproject`. `start`/`length` paging was used until `recordsTotal` was covered. GSA-Human is the database id `hra`; `gsa-human` is not the search id.

HRA006761 file rows came from `GET https://ngdc.cncb.ac.cn/gsa-human/ajaxb/runinstudy?accession=HRA006761&pageNo=&pageSize=10`. That endpoint returns the public catalog (run, file name, size, experiment title). It does not return sequence.

## Papers

PubMed esearch/efetch for CLDN4 or claudin-4 plus knockdown, knockout, shRNA, or siRNA, plus lung, ovarian, or breast, with a China affiliation filter. Abstracts were used to separate a true CLDN4 knockdown from a paper that only measures CLDN4 after silencing another gene. GEO series headers (`acc.cgi?targ=self&form=text`) supplied the country of GSE22493, GSE207704, and GSE50927.

## Reproduce

```bash
python3 methods/node_gsa_cldn4_kd_china/inventory.py
```

Outputs: `results/gate.json` and `results/tables/*.tsv`.
