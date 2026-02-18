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
    mag_obj = event.preferred_magnitude() or event.magnitudes[0]
    magnitude = mag_obj.mag

    return origin.time, origin.latitude, origin.longitude, magnitude

def sliding_coherence(x, y, fs, win_len, step_len, seg_len, fmin, fmax):
    """
    Time-dependent, band-averaged coherence using Welch averaging.
    """
    nwin  = int(win_len * fs)
    nstep = int(step_len * fs)
    nseg  = int(seg_len * fs)

    times = []
    coh_vals = []
    top_3 = []

    for start in range(0, len(x) - nwin, nstep):
        xs = x[start:start + nwin]
        ys = y[start:start + nwin]

        f, Cxy = coherence(xs, ys, fs=fs, nperseg=nseg)

        band = (f >= fmin) & (f <= fmax)
        coh_val = Cxy[band].mean()
        coh_vals.append(coh_val)
        top_3.append(coh_val)
        top_3.sort(reverse=True)
        if len(top_3) > 3:
            top_3.pop()

        times.append((start + nwin / 2) / fs)

    return np.array(times), np.array(coh_vals), top_3

def mag_to_range(magnitude):
    if magnitude >= 6:
        dist_km = 500
    elif magnitude >= 5:
        dist_km = 300
    elif magnitude >= 4:
        dist_km = 200
    elif magnitude >= 3:
        dist_km = 100
    else:
        dist_km = 100

    radius_deg = round((dist_km / 111.19), 2)

    return radius_deg

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
    ev_time, ev_lat, ev_long, ev_mag = get_event_info(client, args.eventid)

    #Time starts five minutes before event
    starttime = ev_time - 300
    endtime = ev_time + (args.minutes * 60)

    args.radius = mag_to_range(ev_mag)
    
    print(f"Searching stations within {args.radius} degrees...")

    inventory = client.get_stations(
        latitude=ev_lat, 
        longitude=ev_long, 
        maxradius=args.radius, 
        level="station",
        channel="BN?,HN?,BH?,HH?"
    )

    per_station_data = []

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
                    channel="BNN,BNE,BNZ,BHN,BHE,BHZ,HNN,HNE,HNZ,HHN,HHE,HHZ",
                    starttime=starttime,
                    endtime=endtime,
                    attach_response=True
                )
                temp = st.copy()

                # Group traces by channel prefix (channel minus last character), e.g., 'BH' from 'BHN'
                prefix_map = {}
                for tr in temp:
                    ch = tr.stats.channel
                    if not ch:
                        continue
                    comp = ch[-1].upper()
                    prefix = ch[:-1]
                    if prefix not in prefix_map:
                        prefix_map[prefix] = {}
                    prefix_map[prefix][comp] = tr

                coherence_entries = []

                # Build map of component -> list of traces (across prefixes)
                comp_map = {}
                for prefix, comps in prefix_map.items():
                    for comp, tr in comps.items():
                        comp_map.setdefault(comp, []).append((prefix, tr))

                # For each component (N/E/Z), compute coherence between traces of the same component
                for comp, tr_list in comp_map.items():
                    # need at least two traces to compute coherence
                    if len(tr_list) < 2:
                        continue
                    # compute coherence for all unique pairs among traces with same component
                    for i in range(len(tr_list)):
                        for j in range(i + 1, len(tr_list)):
                            pref_a, tr_a = tr_list[i]
                            pref_b, tr_b = tr_list[j]
                            tr_a = tr_a.copy()
                            tr_b = tr_b.copy()
                            try:
                                tr_a.remove_response(output="VEL", water_level=60)
                            except Exception:
                                pass
                            try:
                                tr_b.remove_response(output="VEL", water_level=60)
                            except Exception:
                                pass

                            for tr in (tr_a, tr_b):
                                tr.detrend("demean")
                                tr.detrend("linear")
                                tr.taper(0.05)

                            fs = min(tr_a.stats.sampling_rate, tr_b.stats.sampling_rate)

                            try:
                                tr_a.resample(fs)
                                tr_b.resample(fs)
                            except Exception:
                                pass

                            start = max(tr_a.stats.starttime, tr_b.stats.starttime)
                            end = min(tr_a.stats.endtime, tr_b.stats.endtime)
                            if end <= start:
                                continue

                            tr_a.trim(start, end)
                            tr_b.trim(start, end)

                            times, coh_ts, top_3 = sliding_coherence(
                                tr_a.data,
                                tr_b.data,
                                fs=fs,
                                win_len=60.0,
                                step_len=5.0,
                                seg_len=10.0,
                                fmin=0.5,
                                fmax=5.0
                            )

                            avg_coh = sum(top_3) / len(top_3) if len(top_3) > 0 else 0

                            label = f"{tr_a.id} vs {tr_b.id}"
                            coherence_entries.append({
                                "times": times,
                                "values": coh_ts,
                                "label": label,
                                "avg_coh": avg_coh
                            })

                if not coherence_entries:
                    # no usable pairs for this station
                    continue

                print(f"Requested {network.code}.{station.code}...")

                station_package = {
                    "station_id": f"{network.code}.{station.code}",
                    "waveforms": temp,
                    "coherence_entries": coherence_entries,
                }

                per_station_data.append(station_package)

            except Exception:
                continue

    if not per_station_data:
        print("No waveform data found for the given parameters.")
        return

    # Processing

    for entry in per_station_data:
        st = entry["waveforms"]
        st.merge(method=1, fill_value='interpolate')
        st.detrend("demean")
        st.detrend("linear")
        st.taper(max_percentage=0.05)

    # Plotting

    output_dir = f"event_{args.eventid}_plots"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    plt.style.use("ggplot")
    for entry in per_station_data:
        station_name = entry["station_id"]
        out_png = os.path.join(output_dir, f"station_{station_name}_coherence.png")

        fig, ax_coh = plt.subplots(
            figsize=(12, 6)
        )

        # ---- Coherence vs time ----
        for coh in entry.get("coherence_entries", []):
            ax_coh.plot(
                coh["times"],
                coh["values"],
                lw=1.2,
                label=f"{coh.get('label', 'coherence')} (avg={coh.get('avg_coh', 0):.2f})"
            )

        ax_coh.set_xlabel("Time since start (s)")
        ax_coh.set_ylabel("Coherence")
        ax_coh.set_ylim(0, 1.05)
        ax_coh.set_title(f"Event {args.eventid}, Station {entry['station_id']}: Coherence")
        ax_coh.grid(True)
        ax_coh.legend(fontsize="x-small")

        # ---- Save and show ----
        plt.savefig(out_png, dpi=200, bbox_inches="tight")
        plt.close(fig)

        print(f"Saved coherence plot to {out_png}")

    # Compilation plot - all coherence data together
    fig, ax_comp = plt.subplots(figsize=(14, 8))
    
    for entry in per_station_data:
        station_id = entry["station_id"]
        for coh in entry.get("coherence_entries", []):
            ax_comp.plot(
                coh["times"],
                coh["values"],
                lw=1.0,
                label=f"{station_id}: {coh.get('label', 'coherence')} (avg={coh.get('avg_coh', 0):.2f})"
            )
    
    ax_comp.set_xlabel("Time since start (s)")
    ax_comp.set_ylabel("Coherence")
    ax_comp.set_ylim(0, 1.05)
    ax_comp.set_title(f"Event {args.eventid}: All Station Coherence Compilation")
    ax_comp.grid(True)
    ax_comp.legend(fontsize="x-small", loc="best")
    
    compilation_png = os.path.join(output_dir, f"event_{args.eventid}_coherence_compilation.png")
    plt.savefig(compilation_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    
    print(f"Saved compilation plot to {compilation_png}")


if __name__ == "__main__":
    main()
