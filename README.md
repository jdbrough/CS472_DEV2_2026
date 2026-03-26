# AECQAQC Seismic Coherence Plot Tool

This project is a Python-based QA/QC tool for comparing seismic waveform data from multiple stations associated with a specific earthquake event. The script queries FDSN web services, retrieves waveform data, computes time-dependent coherence between channels, and generates visualization plots.

## Table of Contents

**[Install Instructions](#install-instructions)**

**[Station First Instructions](#usage-instructions-for-station_firstpy)**

**[Event First Instructions](#usage-instructions-for-fdsn_scriptpy)**

### Features:
* Queries seismic event information using FDSN web services
* Automatically selects stations within a radius based on event magnitude
* Fetches broadband and strong-motion waveform data
* Generates PNG visualization plots including:
    *  Time series waveform comparison
    *	Band-averaged coherence vs time
---
## Install Instructions
0. Install python and pip onto your machine.
   Detailed instructions can be found [>here<](https://github.com/PackeTsar/Install-Python?tab=readme-ov-file#windows-) 

1 **(Optional)** Create and activate a venv environment to run script in
   - Create venv environment with
     ```
     python -m venv scipt_env
     ```
   - Activate venv environment
		- Windows
    	```
    	script_env/Scripts/activate
    	```
        - macOS/Linux
    	```
		source script_env/bin/activate
     	``` 
2. Clone Script form Github Repo:
	 ```
	 git clone https://github.com/jdbrough/CS472_DEV2_2026.git
  	 ```
  
3. Install required python libraries from file
	 ```
  	 pip install -r requirements.txt
  	 ```

---
## Usage Instructions for station_first.py

### Basic usage : 
```
python station_first.py
```

You should then be prompted with the following: 
```
Enter a FDSN Station Code (e.g. HDA):
```
Enter an FDSN Station Code to query for events within a relavent range around the station

### The script will:
1. Query All events in the last 6 months within a radius around the station
2. Check the magnitude of each event to see if it's within a valid radius of then
3. Makes a coherency plot for each event that effects the station
4. Makes an aggregation plot showing all of the events
5. Saves all of the plots into an output folder

### Optional flags :
 -h, --help      Show a short help message and exit
  --client        FDSN Client (default IRIS)
  --network       FDSN Network Code (default AK)
  --station       Skip the prompt and use a FDSN Station Code
  --stations      List available stations with both broadband and strong motion channels for the specified client/network


---
## Usage Instructions for fdsn_script.py

### Basic usage :
```
python fdsn_plot.py
```

### The script will:
1.	Query earthquake event 10976411
2.	Find nearby stations
3.	Download waveform data 
4.	Compute coherence
5.	Save plots into an output folder


### Selecting the Earthquake Event

```
python fdsn_plot.py --eventid <EVENT_ID>
```

### The script will:
1.	Query earthquake event `<EVENT_ID>`
2.	Find nearby stations
3.	Download waveform data 
4.	Compute coherence
5.	Save plots into an output folder

---
## fdsn_plot.py Outputs 

Generated plots are saved placed in a generated directory named:
event_<EVENT_ID>_plots/


Each station will produce a PNG file:
station_<NETWORK>.<STATION>_coherence.png

Each figure contains:
•	Top panel: velocity time series
•	Bottom panel: sliding-window coherence

How It Works

1.	Fetch event metadata (time, location, magnitude)
2.	Convert magnitude to search radius
3.	Query station inventory
4.	Download waveform channels
5.	Preprocess signals
6.	Compute band-averaged coherence
7.	Generate and save plots

   
