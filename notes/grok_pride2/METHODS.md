# grok_pride2 methods (short)

- Date: 2026-08-16
- APIs: PRIDE Archive v2 search + project + files; ProteomeXchange GetDataset JSON; OmicsDI dataset search; iProX PROXI `/proxi/datasets/{PXD}`
- Excluded: PXD042091, PXD059688, PXD019061
- Download rule: protein tables only, 80 MB cap, no raw MS
- iProX HTTPS used `ssl` unverified context because `download.iprox.org` presented a hostname-mismatched certificate
- Auto-keep (88) was not trusted; curated list is the analysis set
