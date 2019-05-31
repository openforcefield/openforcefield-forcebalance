#!/usr/bin/env python
from monitor_torsiondrives import TorsionMonitor

monitor = TorsionMonitor('torsion_submit_checkpoint.json', client_conf_file='public_client.yaml')

monitor.sync_from_dataset("OpenFF Group1 Torsions", "default")

monitor.get_update()

monitor.download_complete()
