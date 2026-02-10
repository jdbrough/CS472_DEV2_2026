#!/usr/bin/env python3

import argparse
import os
import sys
import subprocess
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy import Stream
from scipy.signal import coherence
import matplotlib.pyplot as plt
import numpy as np


DEFAULT_ID = "12060740"
DEFAULT_MINUTES = 10
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

def get_event_info(client, event_id):
    catalog = client.get_events(eventid = event_id)
    event = catalog[0]
    origin = event.preferred_origin() or event.origins[0]

    return origin.time, origin.latitude, origin.longitude

def sliding_coherence(x, y, fs, win_len, step_len, seg_len, fmin, fmax):
    """
    Time-dependent, band-averaged coherence using Welch averaging.
    """
    nwin  = int(win_len * fs)
    nstep = int(step_len * fs)
    nseg  = int(seg_len * fs)

    times = []
    coh_vals = []

    for start in range(0, len(x) - nwin, nstep):
        xs = x[start:start + nwin]
        ys = y[start:start + nwin]

        f, Cxy = coherence(xs, ys, fs=fs, nperseg=nseg)

        band = (f >= fmin) & (f <= fmax)
        coh_vals.append(Cxy[band].mean())

        times.append((start + nwin / 2) / fs)

    return np.array(times), np.array(coh_vals)

def main():
    parser = argparse.ArgumentParser(
        description="Query FDSN data and plot waveform to PNG"
    )
    parser.add_argument(
        "--eventid",
        help = "FDSN Event ID (e.g. 11843205)",
        default = DEFAULT_ID
    )
    parser.add_argument(
        "--radius",
        type = float,
        default = 1.0,
        help = "Search radius in degrees around epicenter"
    )
    parser.add_argument(
        "--direction",
        default = "N",
        help = "Direction to fetch (N, E, Z; default: N)"
    )
    parser.add_argument(
        "--client",
        default=DEFAULT_CLIENT,
        help="FDSN client name (default: IRIS)",
    )
    parser.add_argument(
        "--minutes",
        type=float,
        default=DEFAULT_MINUTES,
        help="Minutes of data to fetch (default: 10)",
    )
    parser.add_argument(
        "--interactive",
        action="store_true",
        help="Prompt for SNCL and time window",
    )

    args = parser.parse_args()

    # Interactive overrides
    if args.interactive:
        args.eventid = input(f"Enter an Event ID default: [{args.eventid}]: ").strip() or args.eventid
        
        radius_input = input(f"Enter Search Radius default: [{args.radius}]: ").strip()
        args.radius = float(radius_input) if radius_input else args.radius
        
        args.channel = input(f"Enter Channel default: [{args.channel}]: ").strip() or args.channel
        
        minutes_input = input(f"Enter Minutes default: [{args.minutes}]: ").strip()
        args.minutes = float(minutes_input) if minutes_input else args.minutes

    client = Client(args.client)

    print(f"--- Fetching Event {args.eventid} ---")
    ev_time, ev_lat, ev_long = get_event_info(client, args.eventid)

    #Time starts five minutes before event
    starttime = ev_time - 300
    endtime = ev_time + (args.minutes * 60)

    print(f"Searching stations within {args.radius} degrees...")

    inventory = client.get_stations(
        latitude=ev_lat, 
        longitude=ev_long, 
        maxradius=args.radius, 
        level="station"
    )

    all_data = Stream()
    all_coh_ts = []

    included_networks = ['AK']

    for network in inventory:
        # Certain networks will be identified but the data not publicly accessible so just skip them
        # You may see requests to networks not included, if they don't show on the graph they can probably be excluded
        if network.code not in included_networks:
            continue
        for station in network:
            try:
                st = client.get_waveforms(
                    network=network.code,
                    station=station.code,
                    location="*",
                    channel="*H"+args.direction,
                    starttime=starttime,
                    endtime=endtime,
                    attach_response=True
                )
                temp = st.copy()
                tr_h = st[0]
                st = client.get_waveforms(
                    network=network.code,
                    station=station.code,
                    location="*",
                    channel="*N"+args.direction,
                    starttime=starttime,
                    endtime=endtime,
                    attach_response=True
                )
                temp += st.copy()
                tr_n = st[0]

                if len(temp) >= 2:
                    tr_h.remove_response(output="VEL", water_level=60)
                    tr_n.remove_response(output="VEL", water_level=60)

                    for tr in (tr_h, tr_n):
                        tr.detrend("demean")
                        tr.detrend("linear")
                        tr.taper(0.05)
                    
                    fs = min(tr_h.stats.sampling_rate, tr_n.stats.sampling_rate)

                    tr_h.resample(fs)
                    tr_n.resample(fs)

                    start = max(tr_h.stats.starttime, tr_n.stats.starttime)
                    end = min(tr_h.stats.endtime, tr_n.stats.endtime)

                    tr_h.trim(start, end)
                    tr_n.trim(start, end)

                    times, coh_ts = sliding_coherence(
                        tr_h.data,
                        tr_n.data,
                        fs=fs,
                        win_len=60.0,      # Coherence window in seconds
                        step_len=5.0,      # update every 5 seconds
                        seg_len=10.0,      # FFT segment length
                        fmin=0.5,
                        fmax=5.0
                    )

                    label = f"{tr_h.id} vs {tr_n.id}"
                    all_coh_ts.append((times, coh_ts, label))

                    print(f"Requested {network.code}.{station.code}...")
                    all_data += temp
            except Exception:
                continue

    if not all_data:
        print("No waveform data found for the given parameters.")
        return

    # Processing
    all_data.merge(method=1, fill_value='interpolate')
    all_data.detrend("demean")
    all_data.detrend("linear")
    all_data.taper(max_percentage=0.05)

    # Plotting

    plt.style.use("ggplot")

    fig, (ax_ts, ax_coh) = plt.subplots(
        2, 1,
        figsize=(12, 9),
        sharex=True,
        gridspec_kw={"height_ratios": [2, 1]}
    )

    # ---- Top panel: time series ----
    for tr in all_data:
        ax_ts.plot(
            tr.times(),
            tr.data,
            lw=0.8,
            label=tr.id
        )

    ax_ts.set_ylabel("Velocity")
    ax_ts.set_title(f"Event {args.eventid}: Broadband vs Strong Motion")
    ax_ts.grid(True)

    # Keep legend readable
    if len(all_data) <= 8:
        ax_ts.legend(fontsize="x-small", ncol=2)

    # ---- Bottom panel: coherence vs time ----
    for times, coh_ts, label in all_coh_ts:
        ax_coh.plot(
            times,
            coh_ts,
            lw=1.2,
            label=label
        )

    ax_coh.set_xlabel("Time since start (s)")
    ax_coh.set_ylabel("Coherence")
    ax_coh.set_ylim(0, 1.05)
    ax_coh.grid(True)

    if len(all_coh_ts) <= 6:
        ax_coh.legend(fontsize="x-small", ncol=2)

    # ---- Save and show ----
    out_png = f"event_{args.eventid}_timeseries_coherence.png"
    plt.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved presentation figure to {out_png}")
    open_image(out_png)


if __name__ == "__main__":
    main()
