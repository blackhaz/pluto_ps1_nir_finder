#!/usr/bin/env python3
"""
Create a printable Pluto finder chart using Pan-STARRS DR1 only:
  - Pluto apparent geocentric/topocentric position at a user-specified date/time
  - Pan-STARRS DR1 stellar positions (RA/Dec)
  - Pan-STARRS DR1 z-band magnitudes for star-dot scaling

This version does NOT query APASS.

By default the chart is centered on Pluto at the requested time.
Use --center-ra/--center-dec (or --center-datetime) to keep a fixed field center
across times so Pluto motion becomes obvious between charts.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone

import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord, EarthLocation, get_body, solar_system_ephemeris
from astropy.time import Time
from matplotlib.patches import Ellipse

PS1_CATALOG = "II/349/ps1"
PLUTO_RADIUS = 1188.3 * u.km
STAR_DIAMETER_SCALE = 2.0


def parse_time(dt_text: str) -> Time:
    """Parse ISO-8601 datetime. Naive values are treated as UTC."""
    if dt_text.endswith("Z"):
        dt_text = dt_text[:-1] + "+00:00"
    dt = datetime.fromisoformat(dt_text)
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    return Time(dt)


def pluto_apparent(obstime: Time, location: EarthLocation | None = None) -> SkyCoord:
    """
    Get Pluto apparent coordinates as seen from Earth (GCRS).

    This is the coordinate type to compare with planetarium software views.
    """
    last_error = None
    for eph in ("de440s", "de440", "de432s"):
        try:
            with solar_system_ephemeris.set(eph):
                return get_body("pluto", obstime, location=location)
        except Exception as exc:
            last_error = exc
            continue

    raise RuntimeError(f"Could not compute Pluto position with available ephemerides: {last_error}")


def pluto_apparent_diameter_arcsec(pluto: SkyCoord) -> float:
    """Compute Pluto apparent angular diameter in arcseconds."""
    if not hasattr(pluto, "distance"):
        return float("nan")

    distance_km = pluto.distance.to(u.km)
    ratio = (PLUTO_RADIUS / distance_km).decompose().value
    diameter_rad = 2.0 * np.arctan(ratio) * u.rad
    return float(diameter_rad.to_value(u.arcsec))


def col_as_float(col) -> np.ndarray:
    """Convert an astropy table column (possibly masked) to float ndarray with NaN."""
    arr = np.array(col)
    if np.ma.isMaskedArray(arr):
        arr = arr.filled(np.nan)
    return np.array(arr, dtype=float)


def query_ps1_sources(center: SkyCoord, radius: u.Quantity, zmag_limit: float):
    """Query Pan-STARRS DR1 positions and z magnitudes around the field center."""
    v = Vizier(
        columns=["RAJ2000", "DEJ2000", "zmag"],
        column_filters={"zmag": f"<={zmag_limit:.3f}"},
        row_limit=-1,
    )
    tables = v.query_region(center, radius=radius, catalog=PS1_CATALOG)
    if len(tables) == 0:
        raise RuntimeError("No Pan-STARRS DR1 sources returned for this field.")

    t = tables[0]
    ra = col_as_float(t["RAJ2000"])
    dec = col_as_float(t["DEJ2000"])
    zmag = col_as_float(t["zmag"])

    good = np.isfinite(ra) & np.isfinite(dec) & np.isfinite(zmag)
    if not np.any(good):
        raise RuntimeError("Pan-STARRS query returned no usable (RA, Dec, zmag) rows.")

    coords = SkyCoord(ra=ra[good] * u.deg, dec=dec[good] * u.deg)
    return coords, zmag[good]


def magnitude_to_marker_size(
    zmag: np.ndarray,
    zmag_limit: float,
    output_dpi: int,
) -> np.ndarray:
    """
    Map z-band magnitudes to marker area (points^2) using atlas-like diameters.

    Diameter anchors (in mm):
      mag 0 -> 20 mm
      mag 1 -> 10 mm
      mag 3 -> 5 mm
      mag 4 -> 2 mm

    All diameters are multiplied by STAR_DIAMETER_SCALE.

    For magnitudes fainter than 4, diameter decreases smoothly (log-space interpolation)
    down to STAR_DIAMETER_SCALE output pixels at zmag_limit.
    """
    zmag = np.asarray(zmag, dtype=float)

    scale = STAR_DIAMETER_SCALE

    # Scaled output-pixel floor converted to mm at the chosen save DPI.
    min_diam_mm = scale * 25.4 / float(output_dpi)

    if zmag_limit <= 4.0:
        mag_nodes = np.array([0.0, 1.0, 3.0, 4.0], dtype=float)
        diam_nodes_mm = np.array([20.0 * scale, 10.0 * scale, 5.0 * scale, 2.0 * scale], dtype=float)
    else:
        mag_nodes = np.array([0.0, 1.0, 3.0, 4.0, zmag_limit], dtype=float)
        diam_nodes_mm = np.array([20.0 * scale, 10.0 * scale, 5.0 * scale, 2.0 * scale, min_diam_mm], dtype=float)

    # Interpolate diameters in log-space so size changes feel atlas-like and smooth.
    mags = np.clip(zmag, mag_nodes[0], mag_nodes[-1])
    log_d = np.interp(mags, mag_nodes, np.log10(diam_nodes_mm))
    diam_mm = 10.0 ** log_d

    # Convert diameter in mm -> area in points^2 for matplotlib scatter(s=...)
    diam_pt = diam_mm * 72.0 / 25.4
    area_pt2 = np.pi * (diam_pt / 2.0) ** 2

    # Enforce at least STAR_DIAMETER_SCALE pixels diameter in output.
    min_diam_pt = STAR_DIAMETER_SCALE * 72.0 / float(output_dpi)
    min_area_pt2 = np.pi * (min_diam_pt / 2.0) ** 2
    return np.maximum(area_pt2, min_area_pt2)


def unwrap_ra_around(ra_deg: np.ndarray, center_deg: float) -> np.ndarray:
    """Wrap RA values into a continuous range around a center RA."""
    ra_deg = np.asarray(ra_deg, dtype=float)
    return ((ra_deg - center_deg + 180.0) % 360.0) - 180.0 + center_deg


def format_ra_hhmm(ra_deg: float) -> str:
    """Format RA in degrees to hh:mm (wrapped to 0..24h)."""
    ra_deg = ra_deg % 360.0
    total_minutes = ra_deg / 15.0 * 60.0
    hh = int(total_minutes // 60) % 24
    mm = int(round(total_minutes % 60))
    if mm == 60:
        hh = (hh + 1) % 24
        mm = 0
    return f"{hh:02d}:{mm:02d}"


def format_dec_ddmm(dec_deg: float) -> str:
    """Format declination in degrees to ±dd°mm'."""
    sign = "+" if dec_deg >= 0 else "-"
    val = abs(dec_deg)
    dd = int(val)
    mm = int(round((val - dd) * 60.0))
    if mm == 60:
        dd += 1
        mm = 0
    return f"{sign}{dd:02d}°{mm:02d}'"


def add_orientation_arrows(ax):
    """Add North and East arrows directly on the chart."""
    x0, y0 = 0.12, 0.12

    ax.annotate(
        "",
        xy=(x0, y0 + 0.10),
        xytext=(x0, y0),
        xycoords=ax.transAxes,
        arrowprops=dict(arrowstyle="-|>", color="black", lw=1.2),
    )
    ax.text(
        x0,
        y0 + 0.11,
        "N",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        color="black",
        fontsize=9,
    )

    # East points left on sky charts with RA increasing to the left.
    ax.annotate(
        "",
        xy=(x0 - 0.10, y0),
        xytext=(x0, y0),
        xycoords=ax.transAxes,
        arrowprops=dict(arrowstyle="-|>", color="black", lw=1.2),
    )
    ax.text(
        x0 - 0.11,
        y0,
        "E",
        transform=ax.transAxes,
        ha="right",
        va="center",
        color="black",
        fontsize=9,
    )


def make_plot(
    pluto: SkyCoord,
    chart_center: SkyCoord,
    obstime: Time,
    star_coords: SkyCoord,
    zmag: np.ndarray,
    zmag_limit: float,
    fov_deg: float,
    pluto_diameter_arcsec: float,
    output: str,
    output_dpi: int,
):
    """Render and save the finder chart (black-on-white printable style)."""
    ra_center = chart_center.ra.deg
    dec_center = chart_center.dec.deg

    # Convert RA to continuous values around Pluto so axis is not broken at RA=0/24.
    ra_plot = unwrap_ra_around(star_coords.ra.deg, ra_center)
    dec_plot = star_coords.dec.deg
    pluto_ra_plot = unwrap_ra_around(np.array([pluto.ra.deg]), ra_center)[0]
    pluto_dec_plot = pluto.dec.deg

    # Set field limits in absolute RA/Dec.
    dec_half = fov_deg / 2.0
    cos_dec = max(np.cos(np.deg2rad(dec_center)), 0.12)
    ra_half = dec_half / cos_dec

    ra_min = ra_center - ra_half
    ra_max = ra_center + ra_half
    dec_min = dec_center - dec_half
    dec_max = dec_center + dec_half

    inside = (
        (ra_plot >= ra_min)
        & (ra_plot <= ra_max)
        & (dec_plot >= dec_min)
        & (dec_plot <= dec_max)
        & np.isfinite(zmag)
    )
    pluto_in_field = (ra_min <= pluto_ra_plot <= ra_max) and (dec_min <= pluto_dec_plot <= dec_max)

    ra_plot = ra_plot[inside]
    dec_plot = dec_plot[inside]
    zmag = zmag[inside]

    if len(zmag) == 0:
        raise RuntimeError("No stars remain inside the requested chart bounds.")

    sizes = magnitude_to_marker_size(zmag, zmag_limit=zmag_limit, output_dpi=output_dpi)

    fig, ax = plt.subplots(figsize=(8.2, 8.2), facecolor="white")
    ax.set_facecolor("white")

    # Stars: black opaque dots, size from z-band magnitude.
    ax.scatter(ra_plot, dec_plot, s=sizes, c="black", alpha=1.0, linewidths=0, zorder=2)

    # Pluto true apparent-size marker (ellipse in RA/Dec coordinates).
    if np.isfinite(pluto_diameter_arcsec) and pluto_diameter_arcsec > 0:
        pluto_diam_dec_deg = pluto_diameter_arcsec / 3600.0
        pluto_diam_ra_deg = pluto_diam_dec_deg / cos_dec
        ax.add_patch(
            Ellipse(
                (pluto_ra_plot, pluto_dec_plot),
                width=pluto_diam_ra_deg,
                height=pluto_diam_dec_deg,
                facecolor="#c62828",
                edgecolor="black",
                linewidth=0.5,
                zorder=6,
            )
        )

    # Small solid marker for visibility at normal chart scales.
    ax.scatter(
        [pluto_ra_plot],
        [pluto_dec_plot],
        s=20,
        c="#c62828",
        edgecolors="black",
        linewidths=0.5,
        zorder=7,
    )

    # Astronomical convention: East to the left => RA increases to the left.
    ax.set_xlim(ra_max, ra_min)
    ax.set_ylim(dec_min, dec_max)

    # Real RA/Dec grid and coordinate labels.
    ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=7))
    ax.yaxis.set_major_locator(mticker.MaxNLocator(nbins=7))
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, pos: format_ra_hhmm(x)))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda y, pos: format_dec_ddmm(y)))
    ax.grid(which="major", color="0.65", alpha=0.9, linestyle="-", linewidth=0.6, zorder=1)

    ax.set_xlabel("Right Ascension (J2000)", color="black")
    ax.set_ylabel("Declination (J2000)", color="black")
    ax.tick_params(colors="black")
    for spine in ax.spines.values():
        spine.set_color("black")

    add_orientation_arrows(ax)

    title = (
        f"Pluto Finder Chart  |  {obstime.utc.isot} UTC\n"
        f"Pan-STARRS DR1 positions and z-band magnitudes"
    )
    ax.set_title(title, color="black", fontsize=11)

    info = (
        f"Center RA: {format_ra_hhmm(chart_center.ra.deg)}\n"
        f"Center Dec: {format_dec_ddmm(chart_center.dec.deg)}\n"
        f"Pluto RA:  {format_ra_hhmm(pluto.ra.deg)}\n"
        f"Pluto Dec: {format_dec_ddmm(pluto.dec.deg)}\n"
        f"Pluto diameter: {pluto_diameter_arcsec:.3f}\"\n"
        f"Pluto in field: {'yes' if pluto_in_field else 'no'}\n"
        f"Stars plotted: {len(zmag)}"
    )
    ax.text(
        0.99,
        0.01,
        info,
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        color="black",
        fontsize=8,
        bbox=dict(facecolor="white", edgecolor="black", alpha=0.85, boxstyle="round,pad=0.3"),
    )

    fig.tight_layout()
    fig.savefig(output, dpi=output_dpi, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Plot a printable Pluto finder chart with Pan-STARRS DR1 stars (z-band sizing)."
    )
    parser.add_argument(
        "--datetime",
        required=True,
        help="Observation time (ISO-8601), e.g. 2026-07-15T22:30:00+00:00",
    )
    parser.add_argument(
        "--fov",
        type=float,
        default=1.2,
        help="Field of view width in degrees (default: 1.2)",
    )
    parser.add_argument(
        "--zmag-limit",
        type=float,
        default=18.0,
        help="Faint z-band magnitude limit (default: 18.0)",
    )
    parser.add_argument(
        "--output",
        default="pluto_finder_ps1_z_printable.png",
        help="Output image path",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=240,
        help="Output DPI (default: 240). Used for faint-star minimum diameter scaling.",
    )
    parser.add_argument(
        "--center-ra",
        type=float,
        default=None,
        help="Optional fixed chart center RA in degrees. Use with --center-dec.",
    )
    parser.add_argument(
        "--center-dec",
        type=float,
        default=None,
        help="Optional fixed chart center Dec in degrees. Use with --center-ra.",
    )
    parser.add_argument(
        "--center-datetime",
        default=None,
        help="Optional fixed center from Pluto position at this ISO datetime. Ignored if --center-ra/--center-dec are set.",
    )
    parser.add_argument(
        "--observer-lat",
        type=float,
        default=None,
        help="Observer latitude in degrees (optional; use with --observer-lon for topocentric Pluto).",
    )
    parser.add_argument(
        "--observer-lon",
        type=float,
        default=None,
        help="Observer longitude in degrees (east positive; optional; use with --observer-lat).",
    )
    parser.add_argument(
        "--observer-elevation-m",
        type=float,
        default=0.0,
        help="Observer elevation in meters (default: 0).",
    )
    args = parser.parse_args()

    if (args.observer_lat is None) != (args.observer_lon is None):
        raise ValueError("Provide both --observer-lat and --observer-lon, or neither.")

    location = None
    if args.observer_lat is not None:
        location = EarthLocation.from_geodetic(
            lon=args.observer_lon * u.deg,
            lat=args.observer_lat * u.deg,
            height=args.observer_elevation_m * u.m,
        )

    obstime = parse_time(args.datetime)
    pluto = pluto_apparent(obstime, location=location)
    pluto_diam_arcsec = pluto_apparent_diameter_arcsec(pluto)

    if (args.center_ra is None) != (args.center_dec is None):
        raise ValueError("Provide both --center-ra and --center-dec, or neither.")

    if args.center_ra is not None:
        chart_center = SkyCoord(ra=args.center_ra * u.deg, dec=args.center_dec * u.deg, frame="icrs")
    elif args.center_datetime is not None:
        center_time = parse_time(args.center_datetime)
        chart_center = pluto_apparent(center_time, location=location)
    else:
        # Default behavior: chart follows Pluto (so Pluto stays near center)
        chart_center = pluto

    # Circumscribed radius so corners of the square field are covered.
    query_radius = np.hypot(args.fov / 2.0, args.fov / 2.0) * u.deg

    # Query stars using the chart center's angular position only (no distance/origin translation).
    query_center = SkyCoord(ra=chart_center.ra.deg * u.deg, dec=chart_center.dec.deg * u.deg, frame="icrs")
    ps1_coords, zmag = query_ps1_sources(query_center, query_radius, args.zmag_limit)

    make_plot(
        pluto=pluto,
        chart_center=chart_center,
        obstime=obstime,
        star_coords=ps1_coords,
        zmag=zmag,
        zmag_limit=args.zmag_limit,
        fov_deg=args.fov,
        pluto_diameter_arcsec=pluto_diam_arcsec,
        output=args.output,
        output_dpi=args.dpi,
    )

    print(f"Saved finder chart: {args.output}")


if __name__ == "__main__":
    main()
