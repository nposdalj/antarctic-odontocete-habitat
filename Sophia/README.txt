This dataset contains daily marine mammal acoustic detections and associated environmental variables collected near Antarctica. Each row represents one day of observations at a specific recording site.

General Columns
Column	Description
date	Date of the observation
Site	Acoustic recording site location (e.g., EI = Elephant Island, KGI = King George Island, CI = Clarence Island)

Marine Mammal Detection Columns

These columns represent daily acoustic detections or presence of different marine mammal species/call types.

Column	Description
BW29	Presence/count of BW29 beaked whale click type (associated with Southern bottlenose whales)
BW37	Presence/count of BW37 beaked whale click type
BW58	Presence/count of BW58 beaked whale click type
Gm	Long-finned pilot whale detections (Globicephala melas)
Pm	Sperm whale detections (Physeter macrocephalus)
Oo	Killer whale detections (Orcinus orca)
Bp	Fin whale acoustic detections/index (Balaenoptera physalus)
Bm	Blue whale acoustic detections/index (Balaenoptera musculus)
Mn	Humpback whale detections (Megaptera novaeangliae)

Environmental Variable Columns

These columns describe oceanographic and environmental conditions corresponding to the same day and location.

Column	Description
temperature_0	Sea surface temperature
salinity_0	Sea surface salinity
EKE_0	Eddy kinetic energy (measure of ocean turbulence and circulation variability)
chla_0	Chlorophyll-a concentration (proxy for phytoplankton productivity)
productivity_0	Primary productivity estimate
o2_0	Dissolved oxygen concentration
EKE_mad_0	Variability in eddy kinetic energy
SSH	Sea surface height
mixed_layer	Mixed layer depth
ice_conc	Sea ice concentration
ice_thickness	Sea ice thickness
FSLE	Finite-Size Lyapunov Exponent (used to identify ocean fronts and mixing structures)