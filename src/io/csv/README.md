# CSV layer

`csv_reader.f90` — generic header-validated CSV parser (untyped `(:,:)` output).

`meteo_csv.f90` — typed table records for meteo CSVs (pilot for ADR 0044).
Each family wraps `csv_reader` and exposes typed `rows(:)` with named fields.
