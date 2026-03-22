#!/usr/bin/env python3

import argparse
import os
import sys
import numpy as np

DEFAULT_MINUTES = 5
DEFAULT_NETWORK = "AK"
DEFAULT_CLIENT = "IRIS"
DEFAULT_STATION = "HDA"
DEFAULT_STATIONS = False

def get_station_info(client, station):
    # Check if station exists and has BN*, HN*, BH*, HH* channels
    try:        
        catalog = client.get_stations(
            station=station,
            channel="BN?,HN?,BH?,HH?"
        )
    except Exception as e:
        print(f"Error fetching station info: {e}")
        sys.exit(1)

    for net in catalog:
        for sta in net:
            latitude = sta.latitude
            longitude = sta.longitude

    return latitude, longitude

def get_lat_and_long_bounds(lat, long):
    lat_N = lat + 4.5
    lat_S = lat - 4.5
    long_E = long + 3.5
    long_W = long - 4.5

    return lat_N, lat_S, long_E, long_W

def calculate_data(temp, eventID):
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

    print(f"Requested {eventID}")

    station_package = {
        "event_id": f"{eventID}",
        "waveforms": temp,
        "coherence_entries": coherence_entries,
    }

    return station_package

def sliding_coherence(x, y, fs, win_len, step_len, seg_len, fmin, fmax):
    from scipy.signal import coherence
    #Time-dependent, band-averaged coherence using Welch averaging.
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
        add_help=False,
        description="Query FDSN data and plot waveform to PNG"
    )
    parser.add_argument(
        "-h", "--help",
        action="store_true",
        help="Show short help message and exit"
    )
    parser.add_argument(
        "--client",
        help = "FDSN Client (e.g. IRIS, AK, NCEDC)",
        default = DEFAULT_CLIENT
    )
    parser.add_argument(
        "--network",
        help = "FDSN Network Code (e.g. AK, IRIS)",
        default = DEFAULT_NETWORK
    )
    parser.add_argument(
        "--station",
        help = "FDSN Station Code (e.g. ANMO)",
        default = DEFAULT_STATION
    )
    parser.add_argument(
        "--stations",
        action="store_true",
        help = "Comma-separated list of FDSN Station Codes (e.g. ANMO,AFI)"
    )

    args = parser.parse_args()

    if args.help:
        print("station_first.py: query FDSN station data and generate coherence plots for a station.")
        print("Options:")
        print("  -h, --help      Show this short help message and exit")
        print("  --client        FDSN Client (default IRIS)")
        print("  --network       FDSN Network Code (default AK)")
        print("  --station       Skip the prompt and use a FDSN Station Code")
        print("  --stations      List available stations with both broadband and strong motion channels for the specified client/network")
        return
    
    from obspy.clients.fdsn import Client
    if args.stations:
        print("Available stations with both broadband and strong motion channels:")
        catalog = Client(args.client).get_stations(network=args.network, station="*", channel="BN?,HN?,BH?,HH?")
        stations = []
        for net in catalog:
            for sta in net:
                stations.append(sta.code)
        print(", ".join(stations))
        return

    client = Client(args.client)

    print(f"--- Fetching Data for Station {args.station} ---")
    st_lat, st_long = get_station_info(client, args.station)

    lat_N, lat_S, long_E, long_W = get_lat_and_long_bounds(st_lat, st_long)
    
    print(f"Searching events within 4.5 degrees...")

    from obspy import UTCDateTime

    inventory = client.get_events(
        minlatitude=lat_S,
        maxlatitude=lat_N,
        minlongitude=long_W,
        maxlongitude=long_E,
        minmagnitude=2.0,
        maxmagnitude=7.0,
        starttime=UTCDateTime.now() - 90 * 24 * 3600,  # last 3 months
    )

    per_event_data = []

    for event in inventory:

        #Skip events that are too far away to have data available using magnitude as a proxy for distance
        ev_lat = event.origins[0].latitude
        ev_long = event.origins[0].longitude
        magnitude = event.magnitudes[0].mag
        radius_deg = mag_to_range(magnitude)
        if abs(ev_lat - st_lat) > radius_deg or abs(ev_long - st_long) > radius_deg:
            continue

        # Certain networks will be identified but the data not publicly accessible so just skip them
        # You may see requests to networks not included, if they don't show on the graph they can probably be excluded
        try:
            id = event.resource_id.id.split('=')[-1]
        except AttributeError:
            id = 0

        try:
            st = client.get_waveforms(
                network=args.network,
                station=args.station,
                location="*",
                channel="BNN,BNE,BNZ,BHN,BHE,BHZ,HNN,HNE,HNZ,HHN,HHE,HHZ",
                attach_response=True,
                starttime=event.origins[0].time - (DEFAULT_MINUTES * 20),
                endtime=event.origins[0].time + (DEFAULT_MINUTES * 40)
            )
            temp = st.copy()

            per_event_data.append(calculate_data(temp, id))
        except Exception as e:
            continue

    if not per_event_data:
        print("No waveform data found for the given parameters.")
        return
    
    for entry in per_event_data:
        print(entry["event_id"])

    # Processing

    for entry in per_event_data:
        st = entry["waveforms"]
        st.merge(method=1, fill_value='interpolate')
        st.detrend("demean")
        st.detrend("linear")
        st.taper(max_percentage=0.05)

    # Plotting
    import matplotlib.pyplot as plt

    output_dir = f"station_{args.station}_plots"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    plt.style.use("ggplot")
    for entry in per_event_data:
        event_name = entry["event_id"]
        out_png = os.path.join(output_dir, f"station_{event_name}_coherence.png")

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
        ax_coh.set_title(f"Station {args.station}, Event {entry['event_id']}: Coherence")
        ax_coh.grid(True)
        ax_coh.legend(fontsize="x-small")

        # ---- Save and show ----
        plt.savefig(out_png, dpi=200, bbox_inches="tight")
        plt.close(fig)

        print(f"Saved coherence plot to {out_png}")

    # Compilation plot - all coherence data together
    fig, ax_comp = plt.subplots(figsize=(14, 8))
    
    for entry in per_event_data:
        event_id = entry["event_id"]
        for coh in entry.get("coherence_entries", []):
            ax_comp.plot(
                coh["times"],
                coh["values"],
                lw=1.0,
                label=f"{event_id}: {coh.get('label', 'coherence')} (avg={coh.get('avg_coh', 0):.2f})"
            )
    
    ax_comp.set_xlabel("Time since start (s)")
    ax_comp.set_ylabel("Coherence")
    ax_comp.set_ylim(0, 1.05)
    ax_comp.set_title(f"Station {args.station}: All Station Coherence Compilation")
    ax_comp.grid(True)
    ax_comp.legend(fontsize="x-small", loc="best")
    
    compilation_png = os.path.join(output_dir, f"station_{args.station}_coherence_compilation.png")
    plt.savefig(compilation_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    
    print(f"Saved compilation plot to {compilation_png}")

if __name__ == "__main__":
    main()
