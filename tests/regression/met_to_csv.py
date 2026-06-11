"""Convert legacy SWAP .met meteo files to the modern CSV format (ADR 0014).

Legacy .met:  comment lines starting '*', a header
  Station,DD,MM,YYYY,Rad,Tmin,Tmax,Hum,Wind,Rain,ETref,Wet
  then rows: '260',01,01,1980,2530.0,-0.8,2.3,0.594468,2.6,5.800,0.3,0.092
Modern .csv:  date,rad,tmin,tmax,hum,wind,rain,etref,wet
  2002-01-01,3810.0,-3.2,-0.1,0.523764,4.9,0.000,0.4,0.000

Usage: python met_to_csv.py in.met out.csv
"""
import csv, sys

def convert(src, dst):
    rows = []
    with open(src) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('*'):
                continue
            if s.lower().startswith('station'):
                continue  # header
            parts = [p.strip().strip("'\"") for p in s.split(',')]
            if len(parts) < 12:
                continue
            _, dd, mm, yyyy, rad, tmin, tmax, hum, wind, rain, etref, wet = parts[:12]
            date = f"{int(yyyy):04d}-{int(mm):02d}-{int(dd):02d}"
            rows.append([date, rad, tmin, tmax, hum, wind, rain, etref, wet])
    with open(dst, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(["date", "rad", "tmin", "tmax", "hum", "wind", "rain", "etref", "wet"])
        w.writerows(rows)
    return len(rows)

if __name__ == "__main__":
    n = convert(sys.argv[1], sys.argv[2])
    print(f"wrote {n} rows -> {sys.argv[2]}")
