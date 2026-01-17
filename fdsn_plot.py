#!/usr/bin/env python3

import argparse
import os
import sys
import subprocess
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
import matplotlib.pyplot as plt


DEFAULT_SNCL = "IU.ANMO.00.BHZ"
DEFAULT_HOURS = 1
DEFAULT_CLIENT = "IRIS"


def open_image(path):
    """Open an image using the default OS viewer (Windows, macOS, Linux, WSL)."""
    try:
        # WSL detection
        if "microsoft" in os.uname().release.lower():
            # Convert WSL path to Windows path
            win_path = subprocess.check_output(
                ["wslpath", "-w", path],
                text=True
            ).strip()
            subprocess.Popen(["explorer.exe", win_path])
        elif sys.platform.startswith("darwin"):
            subprocess.Popen(["open", path])
        elif os.name == "nt":
            os.startfile(path)
        else:
            subprocess.Popen(["xdg-open", path])
    except Exception as e:
        print(f"Could not open image automatically: {e}")



def parse_sncl(sncl):
    """Split SNCL string into components."""
    try:
        net, sta, loc, cha = sncl.split(".")
        return net, sta, loc, cha
    except ValueError:
        raise ValueError("SNCL must be in NET.STA.LOC.CHA format")


def main():
    parser = argparse.ArgumentParser(
        description="Query FDSN data and plot waveform to PNG"
    )
    parser.add_argument(
        "--sncl",
        help="SNCL code (NET.STA.LOC.CHA)",
        default=DEFAULT_SNCL,
    )
    parser.add_argument(
        "--hours",
        type=float,
        default=DEFAULT_HOURS,
        help="Hours of data to fetch (default: 1)",
    )
    parser.add_argument(
        "--interactive",
        action="store_true",
        help="Prompt for SNCL and time window",
    )
    parser.add_argument(
        "--client",
        default=DEFAULT_CLIENT,
        help="FDSN client name (default: IRIS)",
    )

    args = parser.parse_args()

    # Interactive overrides
    if args.interactive:
        sncl = input(f"SNCL [{DEFAULT_SNCL}]: ").strip() or DEFAULT_SNCL
        hours = input(f"Hours [{DEFAULT_HOURS}]: ").strip()
        hours = float(hours) if hours else DEFAULT_HOURS
    else:
        sncl = args.sncl
        hours = args.hours

    net, sta, loc, cha = parse_sncl(sncl)

    # Time window
    endtime = UTCDateTime.now()
    starttime = endtime - hours * 3600

    print(f"Querying {sncl}")
    print(f"Time window: {starttime} → {endtime}")

    client = Client(args.client)

    try:
        st = client.get_waveforms(
            network=net,
            station=sta,
            location=loc,
            channel=cha,
            starttime=starttime,
            endtime=endtime,
        )
    except Exception as e:
        print(f"Error fetching data: {e}")
        sys.exit(1)

    if len(st) == 0:
        print("No data returned.")
        sys.exit(1)

    # Prepare data
    st.merge()
    st.detrend("demean")
    st.detrend("linear")

    tr = st[0]

    times = tr.times("matplotlib")
    data = tr.data

    # Plot
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot_date(times, data, "-", linewidth=0.8)

    ax.set_title(tr.id)
    ax.set_xlabel("Time (UTC)")
    ax.set_ylabel("Counts")

    fig.autofmt_xdate()

    out_png = f"{sncl.replace('.', '_')}.png"
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)


    print(f"Saved plot to {out_png}")
    open_image(out_png)


if __name__ == "__main__":
    main()
