# Pluto PS1 NIR Finder

`pluto_ps1_nir_finder.py` creates a **printable finder chart** for Pluto using:

- Pluto position at a chosen date/time (Astropy ephemerides)
- Pan-STARRS DR1 stars (`II/349/ps1`) for RA/Dec and brightness
- **z-band magnitudes only** for star size scaling
- Black-on-white atlas style with RA/Dec grid, N/E arrows, and Pluto marker

---

## 1) Get the code from GitHub

Repository:

- https://github.com/blackhaz/pluto_ps1_nir_finder

Clone with Git:

```bash
git clone https://github.com/blackhaz/pluto_ps1_nir_finder.git
cd pluto_ps1_nir_finder
```

If you do not use Git, you can also download ZIP from GitHub:

1. Open the repository URL above
2. Click **Code** → **Download ZIP**
3. Extract ZIP
4. Open terminal in extracted folder

---

## 2) Requirements

- Python **3.10+** (3.11 recommended)
- Internet access (queries Pan-STARRS through CDS VizieR)
- Python packages:
  - `numpy`
  - `matplotlib`
  - `astropy`
  - `astroquery`

---

## 3) Install on Windows

### 3.1 Install system tools

1. Install Python from https://www.python.org/downloads/windows/
   - During install, check **Add Python to PATH**.
2. (Optional) Install Git: https://git-scm.com/download/win

### 3.2 Create virtual environment and install dependencies

Open **PowerShell** in the project folder and run:

```powershell
py -3 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install numpy matplotlib astropy astroquery
```

> If PowerShell blocks activation, run once:
>
> ```powershell
> Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
> ```

### 3.3 Test installation

```powershell
python pluto_ps1_nir_finder.py --help
```

---

## 4) Install on macOS

### 4.1 Install system tools

If needed:

```bash
xcode-select --install
```

Install Homebrew (optional but recommended):

- https://brew.sh

Install Python + Git:

```bash
brew install python git
```

### 4.2 Create virtual environment and install dependencies

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install numpy matplotlib astropy astroquery
```

### 4.3 Test installation

```bash
python pluto_ps1_nir_finder.py --help
```

---

## 5) Install on Linux

### 5.1 Install system tools

#### Debian/Ubuntu

```bash
sudo apt update
sudo apt install -y python3 python3-venv python3-pip git
```

#### Fedora

```bash
sudo dnf install -y python3 python3-pip git
```

#### Arch

```bash
sudo pacman -S --needed python python-pip git
```

### 5.2 Create virtual environment and install dependencies

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install numpy matplotlib astropy astroquery
```

### 5.3 Test installation

```bash
python pluto_ps1_nir_finder.py --help
```

---

## 6) Basic usage

Required argument:

- `--datetime` observation date/time in ISO-8601 format

Example:

```bash
python pluto_ps1_nir_finder.py \
  --datetime 2026-04-10T00:00:00+00:00 \
  --fov 3.0 \
  --zmag-limit 15.0 \
  --dpi 240 \
  --output pluto_finder.png
```

This creates `pluto_finder.png` in the current folder.

Topocentric (observer-specific) example:

```bash
python pluto_ps1_nir_finder.py \
  --datetime 2026-04-10T00:00:00+00:00 \
  --observer-lat 52.5200 --observer-lon 13.4050 --observer-elevation-m 35 \
  --fov 1.0 --zmag-limit 17.5 \
  --output pluto_finder_berlin.png
```

---

## 7) CLI options

```text
--datetime         Required. ISO datetime, e.g. 2026-07-15T22:30:00+00:00
--fov              Field width in degrees (default: 1.2)
--zmag-limit       Faintest z magnitude to include (default: 18.0)
--dpi              Output DPI (default: 240)
--output           Output image filename (default: pluto_finder_ps1_z_printable.png)
--center-ra        Optional fixed chart center RA in deg
--center-dec       Optional fixed chart center Dec in deg
--center-datetime  Optional fixed center = Pluto position at this datetime
--observer-lat     Optional observer latitude in deg (for topocentric Pluto)
--observer-lon     Optional observer longitude in deg, east positive
--observer-elevation-m  Optional observer elevation in meters (default: 0)
```

### Centering behavior

- **Default**: chart center follows Pluto at `--datetime`
- For motion studies (same star field, Pluto moves), use:
  - `--center-datetime` (recommended), or
  - `--center-ra` and `--center-dec`

---

## 8) Show Pluto motion between dates (recommended workflow)

Generate two charts with the **same center**:

```bash
python pluto_ps1_nir_finder.py \
  --datetime 2026-04-10T00:00:00+00:00 \
  --center-datetime 2026-04-10T00:00:00+00:00 \
  --fov 0.8 --zmag-limit 17.5 --output chart_2026-04-10.png

python pluto_ps1_nir_finder.py \
  --datetime 2026-05-10T00:00:00+00:00 \
  --center-datetime 2026-04-10T00:00:00+00:00 \
  --fov 0.8 --zmag-limit 17.5 --output chart_2026-05-10.png
```

Now the stars stay fixed and Pluto position shifts.

---

## 9) Star size model used

Star diameters are atlas-style, based on z magnitude:

- mag 0 → 20 mm
- mag 1 → 10 mm
- mag 3 → 5 mm
- mag 4 → 2 mm

Then all diameters are scaled by the script constant:

- `STAR_DIAMETER_SCALE = 2.0`

Faintest stars at `--zmag-limit` are forced to at least **2 output pixels diameter** (because of the 2x scale).

---

## 10) Troubleshooting

### No chart or “No Pan-STARRS DR1 sources returned”

- Increase `--fov` (e.g., `1.2` → `3.0`)
- Use a fainter `--zmag-limit` (e.g., `15` → `18`)
- Check internet/firewall/proxy

### Pluto seems not moving when date changes

This is expected if chart center follows Pluto (default). Use fixed center:

- `--center-datetime <baseline_date>`
  or
- `--center-ra ... --center-dec ...`

### Pluto position does not match Stellarium

Use the same time, location, and coordinate mode when comparing:

- Set `--datetime` exactly (including timezone offset)
- If needed, pass observer site:
  - `--observer-lat`
  - `--observer-lon`
  - `--observer-elevation-m`
- In Stellarium, compare using **J2000/ICRF-like equatorial coordinates** (not apparent-of-date readout)

The script computes Pluto as an Earth-observed apparent position (GCRS) and uses Pan-STARRS J2000 star coordinates.
Small arcsecond-level differences may still occur vs. Stellarium due to ephemeris/model differences.

### `ModuleNotFoundError`

Make sure your virtual environment is activated and dependencies are installed:

```bash
pip install numpy matplotlib astropy astroquery
```

### PowerShell activation blocked

Use:

```powershell
Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
```

---

## 11) Quick verification command

```bash
python pluto_ps1_nir_finder.py \
  --datetime 2026-04-10T00:00:00+00:00 \
  --fov 0.8 --zmag-limit 17.5 \
  --output test_chart.png
```

If successful, the script prints:

```text
Saved finder chart: test_chart.png
```
