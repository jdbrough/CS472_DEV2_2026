# AECQAQC Seismic Coherence Plot Tool

This project is a Python-based QA/QC tool for comparing seismic waveform data from multiple stations associated with a specific earthquake event. The script queries FDSN web services, retrieves waveform data, computes time-dependent coherence between channels, and generates visualization plots.

### Features:
* Queries seismic event information using FDSN web services
* Automatically selects stations within a radius based on event magnitude
* Fetches broadband and strong-motion waveform data
* Generates PNG visualization plots including:
    *  Time series waveform comparison
    *	Band-averaged coherence vs time
---
## Install Instructions
0. **(Optional)** Create and activate a venv environment to run script in
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
1. Clone Script form Github Repo:
	 ```
	 git clone https://github.com/jdbrough/CS472_DEV2_2026.git
  	 ```
  
1. Install required python libraries from file
	 ```
  	 pip install -r requirements.txt
  	 ```

---
## Running the Script

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
## Outputs 

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
