import sys
import os
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                             QHBoxLayout, QLineEdit, QPushButton, QLabel, 
                             QTextEdit, QFileDialog)
from PyQt6.QtCore import QThread, pyqtSignal

class SeismicWorker(QThread):
    """Worker thread to handle FDSN data fetching without freezing the GUI."""
    finished = pyqtSignal(list)
    log = pyqtSignal(str)

    def __init__(self, client_name, network, station):
        super().__init__()
        self.client_name = client_name
        self.network = network
        self.station = station

    def run(self):
        from obspy.clients.fdsn import Client
        from obspy import UTCDateTime
        from concurrent.futures import ThreadPoolExecutor, as_completed

        try:
            self.log.emit(f"Connecting to {self.client_name}...")
            client = Client(self.client_name)
            
            self.log.emit(f"Fetching station {self.station} info...")
            st_lat, st_long = get_station_info(client, self.station)
            
            lat_N, lat_S, long_E, long_W = get_lat_and_long_bounds(st_lat, st_long)
            
            self.log.emit("Searching for events in the last 90 days...")
            inventory = client.get_events(
                minlatitude=lat_S, maxlatitude=lat_N,
                minlongitude=long_W, maxlongitude=long_E,
                minmagnitude=2.0, starttime=UTCDateTime.now() - 90 * 24 * 3600
            )

            per_event_data = []
            with ThreadPoolExecutor(max_workers=8) as executor:
                class Args: pass
                args = Args()
                args.network = self.network
                args.station = self.station

                futures = {executor.submit(fetch_event, ev, client, args, st_lat, st_long): ev for ev in inventory}
                for future in as_completed(futures):
                    res = future.result()
                    if res:
                        per_event_data.append(res)
                        self.log.emit(f"Data retrieved for event: {res['event_id']}")

            self.finished.emit(per_event_data)
        except Exception as e:
            self.log.emit(f"Critical Error: {str(e)}")

class SeismicGui(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Seismic Coherence Analyzer")
        self.setMinimumSize(600, 400)
        self.init_ui()

    def init_ui(self):
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)

        # Input Row
        input_layout = QHBoxLayout()
        self.client_input = QLineEdit("IRIS")
        self.net_input = QLineEdit("AK")
        self.sta_input = QLineEdit("HDA")
        
        input_layout.addWidget(QLabel("Client:"))
        input_layout.addWidget(self.client_input)
        input_layout.addWidget(QLabel("Network:"))
        input_layout.addWidget(self.net_input)
        input_layout.addWidget(QLabel("Station:"))
        input_layout.addWidget(self.sta_input)
        
        # Action Buttons
        self.run_btn = QPushButton("Run Analysis")
        self.run_btn.clicked.connect(self.start_analysis)
        
        # Log Output
        self.log_output = QTextEdit()
        self.log_output.setReadOnly(True)
        self.log_output.setStyleSheet("background-color: #1e1e1e; color: #00ff00; font-family: Courier;")

        layout.addLayout(input_layout)
        layout.addWidget(self.run_btn)
        layout.addWidget(QLabel("Process Log:"))
        layout.addWidget(self.log_output)

    def start_analysis(self):
        self.log_output.clear()
        self.run_btn.setEnabled(False)
        
        # Initialize Worker
        self.worker = SeismicWorker(
            self.client_input.text(),
            self.net_input.text(),
            self.sta_input.text()
        )
        self.worker.log.connect(self.update_log)
        self.worker.finished.connect(self.on_finished)
        self.worker.start()

    def update_log(self, text):
        self.log_output.append(text)

    def on_finished(self, data):
        self.run_btn.setEnabled(True)
        if not data:
            self.log_output.append("Analysis complete. No data found.")
        else:
            self.log_output.append(f"\nSUCCESS: Processed {len(data)} events.")
            self.log_output.append(f"Plots saved to folder: station_{self.sta_input.text()}_plots")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = SeismicGui()
    window.show()
    sys.exit(app.exec())