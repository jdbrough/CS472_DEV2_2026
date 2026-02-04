#!/usr/bin/env python3

import argparse
import os
import sys
import subprocess
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy import Stream
import matplotlib.pyplot as plt


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
        "--channel",
        default = "BHZ",
        help = "Seismic channel to fetch (default: BHZ)"
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
                    channel=args.channel,
                    starttime=starttime,
                    endtime=endtime,
                )
                temp_data = st.copy()
                st = client.get_waveforms(
                    network=network.code,
                    station=station.code,
                    location="*",
                    channel="BNZ",
                    starttime=starttime,
                    endtime=endtime,
                )
                temp_data += st
                if len(temp_data) > 1:
                    print(f"Requesting {network.code}.{station.code}...")
                    all_data += temp_data
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
    plt.style.use('ggplot')
    fig, ax = plt.subplots(figsize=(12, 6))
    
    for tr in all_data:
        ax.plot(tr.times("matplotlib"), tr.data, label=tr.id, lw=1)

    ax.set_title(f"Event {args.eventid} Comparison | Channel: {args.channel}")
    ax.set_ylabel("Counts")
    ax.xaxis_date()
    fig.autofmt_xdate()
    ax.legend(loc='upper right', fontsize='small', ncol=2)

    out_png = f"event_{args.eventid}_comparison.png"
    plt.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"Success! Saved to {out_png}")
    open_image(out_png)


if __name__ == "__main__":
    main()
